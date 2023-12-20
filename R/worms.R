#' Get geometry by MRGID.
#'
#' @export
fetch_geom <- function(mrgid) {
  message(glue("Fetching geom for mrgid {mrgid}"))
  url <- glue::glue("https://marineregions.org/rest/getGazetteerGeometries.jsonld/{mrgid}/")
  geom <- jsonlite::fromJSON(url)$"mr:hasGeometry"$"gsp:asWKT"
  if (!is.null(geom)) {
    geom_fixed <- stringr::str_replace(geom, "<.*>\\s+", "") %>%
    return(geom_fixed)
  }
}

fetch_geom_cached <- memoise(fetch_geom, cache = speedy_disk_cache)

wm_distribution_possibly <- possibly(wm_distribution, otherwise = NULL)

worms_geoms_by_aphiaid <- function(aphiaid) {
  message(glue("Fetching geoms for aphiaid {aphiaid}"))
  # TODO: buffer non polygon types
  dist <- wm_distribution_possibly(aphiaid)
  if (is.null(dist)) {
    return(data.frame())
  }
  dist <- dist %>%
    mutate(mrgid = str_replace_all(locationID, "http://marineregions.org/mrgid/", ""))
  res <- map(1:nrow(dist), function(i) {
    geoms <- fetch_geom_cached(dist$mrgid[i])
    if (!is.null(geoms)) {
      split_geoms <- c()
      for (j in 1:length(geoms)) {
        if (str_detect(geoms[j], "GEOMETRYCOLLECTION")) {
          content <- str_match(geoms[j], "GEOMETRYCOLLECTION\\s?\\((.*)\\)")[,2]
          parts <- str_extract_all(content, "(?:POLYGON|POINT|LINESTRING)\\s?\\([^A-Z]+\\)")[[1]]
          split_geoms <- c(split_geoms, parts)
        } else {
          split_geoms <- c(split_geoms, geoms[j])
        }
      }
      return(data.frame(wkt = split_geoms, locality = dist$locality[i], establishmentMeans = dist$establishmentMeans[i]))
    }
  })
  result <- bind_rows(res)
  return(result)
}

#' Get all synonymized AphiaIDs for an AphiaID.
get_all_aphiaids_by_aphiaid <- function(aphiaid) {
  message(glue("Fetching aphiaids for aphiaid {aphiaid}"))
  accepted <- possibly(worrms::wm_record, otherwise = data.frame(valid_AphiaID = character(0)))(aphiaid) %>% pull(valid_AphiaID)
  synonyms <- possibly(worrms::wm_synonyms, otherwise = data.frame(AphiaID = character(0)))(accepted) %>% pull(AphiaID)
  aphiaids <- unique(c(aphiaid, accepted, synonyms))
  return(aphiaids)
}

#' Get WRiMS distribution.
#'
#' @export
get_worms_dist <- function(scientificname = NULL, aphiaid = NULL, taxonkey = NULL) {
  taxonomy <- insistent_resolve_taxonomy(scientificname, aphiaid, taxonkey)
  stopifnot(is.numeric(taxonomy$aphiaid))
  message(glue("Generating dist for aphiaid {aphiaid}"))

  aphiaids <- get_all_aphiaids_by_aphiaid(taxonomy$aphiaid)
  message("Fetching geoms for all linked aphiaids")
  dist_map <- map(aphiaids, worms_geoms_by_aphiaid)
  if (sum(sapply(dist_map, nrow)) == 0) {
    return(vect())
  }
  message("Processing geometries")
  dist_all <- do.call(rbind, dist_map)
  dist_all <- dist_all[!is.na(dist_all$wkt) & dist_all$wkt != "",]
  wkt_vect <- map(dist_all$wkt, vect)
  which_polygons <- which(map(wkt_vect, geomtype) == "polygons")
  wkt_vect <- wkt_vect[which_polygons]
  if (length(wkt_vect) == 0) {
    return(vect())
  }
  data <- dist_all[which_polygons, c("locality", "establishmentMeans")] %>% mutate(
    establishmentMeans = case_match(
      as.character(establishmentMeans),
      "Alien" ~ "Alien",
      "Native" ~ "Native",
      .default = "Other"
    )
  )
  combined <- vect(wkt_vect)
  values(combined) <- data
  simplified <- simplifyGeom(combined, tolerance = 0.002, preserveTopology = FALSE, makeValid = TRUE) # TODO: check
  simplified_aggregated <- aggregate(disagg(simplified), by = "establishmentMeans", dissolve = TRUE)
  simplified_aggregated_buffered <- buffer(simplified_aggregated, 0.01, quadsegs = 4)
  bbox <- vect("POLYGON ((-180 -90, -180 90, 180 90, 180 -90, -180 -90))")
  simplified_aggregated_buffered_cropped <- crop(simplified_aggregated_buffered, bbox)
  simplified_aggregated_buffered_cropped_noholes <- fillHoles(simplified_aggregated_buffered_cropped)
  simplified_aggregated_buffered_cropped_noholes_simplified <- simplifyGeom(simplified_aggregated_buffered_cropped_noholes, 0.002)
  crs(simplified_aggregated_buffered_cropped_noholes_simplified) <- "epsg:4326"
  layers <- list()
  for (establishment in simplified_aggregated_buffered_cropped_noholes_simplified$establishmentMeans) {
    layers[[establishment]] <- simplified_aggregated_buffered_cropped_noholes_simplified[simplified_aggregated_buffered_cropped_noholes_simplified$establishmentMeans == establishment]
  }
  layers <- layers[order(names(layers), decreasing = TRUE)]
  if (length(layers) > 1) {
    for (i in 2:length(layers)) {
      layers[[1]] <- cover(layers[[1]], layers[[i]])
    }
  }

  return(layers[[1]])
}

get_worms_dist_cached <- memoise(get_worms_dist, cache = speedy_disk_cache_short)

#' Get name for AphiaID.
#'
#' @export
aphiaid_for_name <- function(scientificname) {
  records <- wm_records_name(scientificname, fuzzy = FALSE) %>%
    filter(match_type == "exact")
  return(records$AphiaID[1])
}

#' Get name for AphiaID.
#'
#' @export
name_for_aphiaid <- function(aphiaid) {
  records <- wm_record(aphiaid)
  return(records$scientificname[1])
}
