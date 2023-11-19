#' Get GBIF gridded distribution.
#'
#' @export
get_gbif_dist <- function(scientificname = NULL, aphiaid = NULL, taxonkey = NULL, zoom = 0, square_size = 16) {
  taxonomy <- insistent_resolve_taxonomy(scientificname, aphiaid, taxonkey)
  stopifnot(is.numeric(zoom))
  stopifnot(is.numeric(taxonomy$taxonkey))

  tile_ids <- expand_grid(x = 0:(2^zoom-1), y = 0:(2^zoom-1))
  tiles <- map2(tile_ids$x, tile_ids$y, function(x, y) {
    tile <- mvt_fetch(taxonKey = taxonomy$taxonkey, source = "density", bin = "square", squareSize = square_size, x = x, y = y, z = zoom)
    if (is.null(tile)) {
      empty_sf()
    } else {
      tile %>%
        rename(c(n = total))
    }
  })

  return(vect(do.call(rbind, tiles)))
}

#' @export
name_for_taxonkey <- function(taxonkey) {
  records <- name_usage(key = taxonkey)$data
  return(records$canonicalName[1])
}

#' @export
taxonkey_for_name <- function(scientificname) {
  records <- rgbif::name_backbone(scientificname) %>%
    filter(matchType == "EXACT")
  if (nrow(records) > 0) {
    return(records$usageKey[1])
  } else {
    return(NULL)
  }
}
