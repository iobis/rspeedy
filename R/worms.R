#' Get geometries by MRGID.
#'
#' @export
mr_geometries <- function(mrgids) {
  geoms <- sapply(mrgids, function(mrgid) {
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

#' Get WRiMS distribution.
#'
#' @export
get_worms_dist <- function(scientificname = NULL, aphiaid = NULL, taxonkey = NULL) {
  taxonomy <- resolve_taxonomy(scientificname, aphiaid, taxonkey)
  stopifnot(is.numeric(taxonomy$aphiaid))
  sf_use_s2(FALSE)
  dist <- wm_distribution(taxonomy$aphiaid) %>%
    mutate(mrgid = str_replace_all(locationID, "http://marineregions.org/mrgid/", ""))
  mrgids <- unique(dist$mrgid)
  geoms <- mr_geometries(mrgids)
  spatial <- tibble(geom = geoms, mrgid = names(geom)) %>%
    as.data.frame()
  dist_spatial <- dist %>%
    left_join(spatial, by = "mrgid") %>%
    st_as_sf(crs = 4326) %>%
    st_collection_extract() %>%
    filter(!st_is(., "LINESTRING") & !st_is(., "MULTILINESTRING") & !st_is_empty(.)) %>%
    st_make_valid()
  dist_grouped <- dist_spatial %>%
    nngeo::st_remove_holes(max_area = 500000000) %>%
    mutate(
      establishmentMeans = case_match(
        establishmentMeans,
        "Alien" ~ "Alien",
        .default = "Other"
      )
    ) %>%
    group_by(establishmentMeans) %>%
    summarize() %>%
    nngeo::st_remove_holes(max_area = 500000000)
  dist_simplified <- dist_grouped %>%
    ms_simplify(keep_shapes = TRUE, snap = FALSE)
  dist_alien <- dist_simplified %>% filter(establishmentMeans == "Alien") %>%
    st_make_valid()
  dist_other <- dist_simplified %>% filter(establishmentMeans == "Other") %>%
    st_make_valid() %>%
    st_difference(dist_alien)
  return(bind_rows(dist_alien, dist_other))
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
