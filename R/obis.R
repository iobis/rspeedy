#' Get OBIS gridded distribution.
#'
#' @export
get_obis_dist <- function(scientificname = NULL, aphiaid = NULL, taxonkey = NULL, res = 3) {
  taxonomy <- insistent_resolve_taxonomy(scientificname, aphiaid, taxonkey)
  stopifnot(is.numeric(res))
  stopifnot(is.numeric(taxonomy$aphiaid))
  url <- glue("https://api.obis.org/occurrence/grid/{res}?taxonid={taxonomy$aphiaid}")
  dist <- GET(URLencode(url)) %>%
      content(as = "text") %>%
      geojson_sf() %>%
      vect()
  return(dist)
}

get_obis_dist_safely <- safely(get_obis_dist, otherwise = NULL)
get_obis_dist_possibly <- possibly(get_obis_dist, otherwise = NULL)
