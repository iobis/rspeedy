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
    geoms <- fetch_geom(dist$mrgid[i])
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
  synonyms <- tryCatch({
    worrms::wm_synonyms(aphiaid) %>%
      pull(AphiaID)
  }, error = function(err) {
    NULL
  })
  aphiaids <- unique(c(aphiaid, synonyms))
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
  dist_map <- map(aphiaids, worms_geoms_by_aphiaid)
  if (sum(sapply(dist_map, nrow)) == 0) {
    return(vect())
  }
  dist_all <- do.call(rbind, dist_map)
  dist_all <- dist_all[!is.na(dist_all$wkt) & dist_all$wkt != "",]
  wkt_vect <- map(dist_all$wkt, vect)
  which_polygons <- which(map(wkt_vect, geomtype) == "polygons")
  wkt_vect <- wkt_vect[which_polygons]
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
  aggregated <- aggregate(combined, by = "establishmentMeans", dissolve = TRUE)
  aggregated_buffered <- buffer(aggregated, 0.01)
  bbox <- vect("POLYGON ((-180 -90, -180 90, 180 90, 180 -90, -180 -90))")
  aggregated_buffered_cropped <- crop(aggregated_buffered, bbox)
  aggregated_buffered_cropped_noholes <- fillHoles(aggregated_buffered_cropped)
  aggregated_buffered_cropped_noholes_simplified <- simplifyGeom(aggregated_buffered_cropped_noholes, 0.001)
  crs(aggregated_buffered_cropped_noholes_simplified) <- "epsg:4326"

  layers <- split(aggregated_buffered_cropped_noholes_simplified, "establishmentMeans")
  names(layers) <- sapply(layers, function(x) { x$establishmentMeans })
  layers <- layers[order(names(layers), decreasing = TRUE)]

  if (length(layers) > 1) {
    for (i in 2:length(layers)) {
      layers[[1]] <- cover(layers[[1]], layers[[i]])
    }
  }

  return(layers[[1]])
}

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
