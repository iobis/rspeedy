#' Get geometries by MRGID.
#'
#' @export
mr_geometries <- function(mrgids) {
  geoms <- sapply(mrgids, function(mrgid) {
    message(glue("Fetching geom for mrgid {mrgid}"))
    url <- glue::glue("https://marineregions.org/rest/getGazetteerGeometries.jsonld/{mrgid}/")
    geom <- jsonlite::fromJSON(url)$"mr:hasGeometry"$"gsp:asWKT"
    if (is.null(geom)) {
      return(st_sfc(st_polygon()))
    } else {
      geom_fixed <- stringr::str_replace(geom, "<.*>\\s+", "")
      geom_fixed <- geom_fixed[geom_fixed != ""]
      geom_union <- st_union(st_as_sfc(geom_fixed))
      return(geom_union)
    }
  })
  geoms
}

worms_geoms_by_aphiaid <- function(aphiaid) {
  message(glue("Fetching geoms for aphiaid {aphiaid}"))
  # TODO: buffer non polygon types
  tryCatch({
    dist <- wm_distribution(aphiaid) %>%
      mutate(mrgid = str_replace_all(locationID, "http://marineregions.org/mrgid/", ""))
    mrgids <- unique(dist$mrgid)
    geoms <- mr_geometries(mrgids)
    spatial <- tibble(geom = geoms, mrgid = names(geom)) %>%
      as.data.frame()
    dist_all <- dist %>%
      left_join(spatial, by = "mrgid") %>%
      st_as_sf(crs = 4326)
    dist_collection <- dist_all %>% filter(st_is(., "GEOMETRYCOLLECTION")) %>%
      st_collection_extract()
    dist_not_collection <- dist_all %>% filter(!st_is(., "GEOMETRYCOLLECTION"))
    dist <- rbind(dist_collection, dist_not_collection) %>%
      filter(!st_is(., "MULTIPOINT") & !st_is(., "POINT") & !st_is(., "LINESTRING") & !st_is(., "MULTILINESTRING") & !st_is_empty(.)) %>%
      st_make_valid()
    return(dist)
  }, error = function(err) {
    message(err$message)
    empty_sf()
  })
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
  sf_use_s2(FALSE)

  aphiaids <- get_all_aphiaids_by_aphiaid(taxonomy$aphiaid)

  dist_map <- map(aphiaids, worms_geoms_by_aphiaid)
  if (sum(sapply(dist_map, nrow)) == 0) {
    return(empty_sf())
  }
  dist_all <- do.call(rbind, dist_map)

  dist_grouped <- dist_all %>%
    nngeo::st_remove_holes(max_area = 500000000) %>%
    mutate(
      establishmentMeans = case_match(
        as.character(establishmentMeans),
        "Alien" ~ "Alien",
        "Native" ~ "Native",
        .default = "Other"
      )
    ) %>%
    group_by(establishmentMeans) %>%
    summarize() %>%
    nngeo::st_remove_holes(max_area = 500000000)

  dist_simplified <- dist_grouped %>%
    ms_simplify(keep_shapes = TRUE, snap = FALSE)

  dist_alien <- dist_simplified %>%
    filter(establishmentMeans == "Alien")
  dist_native <- dist_simplified %>%
    filter(establishmentMeans == "Native")
  dist_other <- dist_simplified %>%
    filter(establishmentMeans == "Other")
  if (nrow(dist_alien) > 0) {
    dist_other <- dist_other %>% st_difference(st_geometry(dist_alien))
    dist_native <- dist_native %>% st_difference(st_geometry(dist_alien))
  }
  if (nrow(dist_native) > 0) {
    dist_other <- dist_other %>% st_difference(st_geometry(dist_native))
  }

  return(bind_rows(dist_alien, dist_native, dist_other))
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
