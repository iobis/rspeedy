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
      aphiaid <- tryCatch({
        aphiaid_for_name(scientificname)
      }, error = function(err) {
        NULL
      })
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
  taxonomy <- insistent_resolve_taxonomy(scientificname, aphiaid, taxonkey)
  gbif_dist <- tryCatch({
    get_gbif_dist(scientificname, aphiaid, taxonkey)
  }, error = function(err) {
    message(err)
    empty_sf()
  })
  obis_dist <- tryCatch({
    get_obis_dist(scientificname, aphiaid, taxonkey)
  }, error = function(err) {
    message(err)
    empty_sf()
  })
  return(list(
    taxonomy = taxonomy,
    obis = obis_dist,
    gbif = gbif_dist,
    worms = get_worms_dist(scientificname, aphiaid, taxonkey),
    envelope = get_thermal_envelope(scientificname, aphiaid, taxonkey)
  ))
}

#' @export
plot_dist <- function(dist) {
  p <- ggplot() +
    geom_sf(data = st_as_sfc(st_bbox(c(xmin = -180, ymin = -85, xmax = 180, ymax = 85), crs = 4326)), fill = "#E2E2E0", color = NA)
  if (!is.null(dist$envelope)) {
    if (nrow(dist$envelope$envelope) > 0) {
      p <- p + geom_sf(data = dist$envelope$envelope, aes(fill = "Thermal envelope"), color = NA, alpha = 0.3)
    }
  }
  p <- p + geom_sf(data = landr::get_land_polygons(simplified = 20000), fill = "#fafafc", color = NA)
  if (nrow(dist$worms) > 0) {
    p <- p +
      geom_sf(data = dist$worms %>% filter(establishmentMeans == "Other"), aes(color = "WoRMS"), fill = NA) +
      geom_sf(data = dist$worms %>% filter(establishmentMeans == "Native"), aes(color = "WoRMS (native)"), fill = NA) +
      geom_sf(data = dist$worms %>% filter(establishmentMeans == "Alien"), aes(color = "WoRMS (introduced)"), fill = NA)
  }
  p <- p + geom_sf(data = dist$obis, aes(color = "OBIS"), fill = NA) +
    geom_sf(data = dist$gbif, aes(color = "GBIF"), fill = NA) +
    scale_color_manual(name = "", values = c("GBIF" = "#8934eb", "OBIS" = "#8934eb", "WoRMS" = "#eb3489", "WoRMS (introduced)" = "#d96b04", "WoRMS (native)" = "#32a877")) +
    scale_fill_manual(name = "", values = "#C75386") +
    theme_void() +
    theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5, vjust = 1)) +
    guides(color = guide_legend(override.aes = list(linewidth = 0.2))) +
    coord_sf(expand = 0.1)
  p
}

#' @export
export_dist <- function(dist, file) {
  st_write(dist$obis, file, "obis", driver = "GPKG", delete_dsn = TRUE)
  st_write(dist$gbif, file, "gbif", append = TRUE, driver = "GPKG")
  st_write(dist$worms, file, "worms", append = TRUE, driver = "GPKG")
  if (!is.null(dist$envelope)) {
    st_write(dist$envelope$envelope, file, "envelope", append = TRUE, driver = "GPKG")
  }
}

#' @export
upload_dist_s3 <- function(dist) {
  library("aws.s3")

  aphiaid <- dist$taxonomy$aphiaid
  stopifnot(is.numeric(aphiaid))

  tmp <- tempdir()
  file_gpkg <- file.path(tmp, "speedy.gpkg")
  file_png <- file.path(tmp, "speedy.png")

  export_dist(dist, file_gpkg)
  plot_dist(dist) + ggtitle(dist$taxonomy$scientificname)
  ggsave(file_png, height = 7, width = 12, dpi = 300, bg = "white")

  s3_base <- glue("speedy/{aphiaid}")
  delete_object(s3_base, bucket = "obis-products")
  put_object(file = file_gpkg, object = glue("{s3_base}/speedy.gpkg"), bucket = "obis-products")
  put_object(file = file_png, object = glue("{s3_base}/speedy.png"), bucket = "obis-products")
}

#' @export
download_dist_s3 <- function(aphiaid) {
  library("aws.s3")
  stopifnot(is.numeric(aphiaid))
  taxonomy <- insistent_resolve_taxonomy(aphiaid = aphiaid)
  tmp <- tempdir()
  s3_base <- glue("speedy/{aphiaid}")
  file_gpkg <- file.path(tmp, "speedy.gpkg")
  save_object(object = glue("{s3_base}/speedy.gpkg"), bucket = "obis-products", file = file_gpkg)

  dist_obis <- tryCatch({
    st_read(file_gpkg, layer = "obis")
  }, error = function(err) {
    message(err)
    empty_sf()
  })
  dist_gbif <- tryCatch({
    st_read(file_gpkg, layer = "gbif")
  }, error = function(err) {
    message(err)
    empty_sf()
  })
  dist_worms <- tryCatch({
    st_read(file_gpkg, layer = "worms")
  }, error = function(err) {
    message(err)
    empty_sf()
  })
  envelope <- tryCatch({
    st_read(file_gpkg, layer = "envelope")
  }, error = function(err) {
    message(err)
    empty_sf()
  })
  dist <- list(
    taxonomy = taxonomy,
    obis = dist_obis,
    gbif = dist_gbif,
    worms = dist_worms,
    envelope = envelope
  )

  return(dist)
}

empty_sf <- function() {
  st_sf(st_sfc(), crs = 4326)
}

#' @export
insistent_resolve_taxonomy <- purrr::insistently(resolve_taxonomy, rate = rate_backoff(pause_base = 10, pause_cap = 300, max_times = 3), quiet = FALSE)

