#' @title Marine species distribution summaries using OBIS, GBIF, and WoRMS.
#'
#' @docType package
#' @name speedy
#' @import dplyr httr sf glue geojsonsf stringr worrms tibble rmapshaper mapview rgbif tidyr purrr stars sdmpredictors ks ggplot2 nngeo eks
#' @author Pieter Provoost, \email{p.provoost@unesco.org}
#' @author Silas Principe
NULL

#' @export
resolve_taxonomy <- function(scientificname = NULL, aphiaid = NULL, taxonkey = NULL) {
  if (is.null(scientificname)) {
    if (!is.null(aphiaid)) {
      scientificname <- name_for_aphiaid(aphiaid)
    } else if (!is.null(taxonkey)) {
      scientificname <- name_for_taxonkey(taxonkey)
    }
  }
  if (is.null(aphiaid)) {
    if (!is.null(scientificname)) {
      aphiaid <- aphiaid_for_name(scientificname)
    }
  }
  if (is.null(taxonkey)) {
    if (!is.null(scientificname)) {
      taxonkey <- taxonkey_for_name(scientificname)
    }
  }
  return(list(scientificname = scientificname, aphiaid = aphiaid, taxonkey = taxonkey))
}

#' @export
get_dist <- function(scientificname = NULL, aphiaid = NULL, taxonkey = NULL) {
  return(list(
    obis = get_obis_dist(scientificname, aphiaid, taxonkey),
    gbif = get_gbif_dist(scientificname, aphiaid, taxonkey),
    worms = get_worms_dist(scientificname, aphiaid, taxonkey),
    envelope = get_thermal_envelope(scientificname, aphiaid, taxonkey)
  ))
}

#' @export
plot_dist <- function(dist, title = "") {
  ggplot() +
    geom_sf(data = st_as_sfc(st_bbox(c(xmin = -180, ymin = -85, xmax = 180, ymax = 85), crs = 4326)), fill = "#E2E2E0", color = NA) +
    geom_sf(data = dist$envelope$envelope, aes(fill = "Thermal envelope"), color = NA, alpha = 0.3) +
    geom_sf(data = landr::get_land_polygons(simplified = 20000), fill = "#fafafc", color = NA) +
    geom_sf(data = dist$worms %>% filter(establishmentMeans == "Other"), aes(color = "WoRMS"), fill = NA) +
    geom_sf(data = dist$worms %>% filter(establishmentMeans == "Alien"), aes(color = "WoRMS (introduced)"), fill = NA) +
    geom_sf(data = dist$obis, aes(color = "OBIS"), fill = NA) +
    geom_sf(data = dist$gbif, aes(color = "GBIF"), fill = NA) +
    scale_color_manual(name = "", values = c("#8934eb", "#8934eb", "#eb3489", "#32a877")) +
    scale_fill_manual(name = "", values = "#C75386") +
    theme_void() +
    ggtitle(title) +
    theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5, vjust = 1)) +
    guides(color = guide_legend(override.aes = list(linewidth = 0.2))) +
    coord_sf(expand = 0.1)
}
