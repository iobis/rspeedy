# library(stars)
# library(ks)
# library(rgbif)
# library(mapview)
# library(tidyr)
# library(ggplot2)
# library(eks)

#' @export
view_envelope <- function(envelope, points = FALSE) {
  if (points) {
    mapview(st_as_sf(envelope$envelope)) + mapview(st_as_sf(envelope$points))
  } else {
    mapview(st_as_sf(envelope$envelope))
  }
}

#' @export
plot_envelope <- function(envelope) {
  ggplot() +
    geom_sf(data = st_as_sfc(st_bbox(c(xmin = -180, ymin = -85, xmax = 180, ymax = 85), crs = 4326)), fill = "#E2E2E0", color = NA) +
    geom_sf(data = envelope$envelope, aes(fill = BO22_tempmean_ss), fill = "#C75386", color = NA) +
    geom_sf(data = landr::get_land_polygons(simplified = 20000), fill = "#fafafc", color = NA) +
    theme_void()
}

#' @export
plot_kd <- function(envelope) {
  ggplot() +
    geom_line(data = data.frame(eval.points = envelope$kd$eval.points, estimate = envelope$kd$estimate), aes(eval.points, estimate), col = 4) +
    geom_jitter(data = data.frame(y = -0.03, x = envelope$kd$x), aes(x, y), height = 0.01, size = 0.1) +
    geom_vline(xintercept = envelope$percentiles, col = 4) +
    geom_vline(xintercept = quantile(envelope$temperatures, probs = c(0.01, 0.99)), col = "black", linetype = "dashed") +
    xlab("temperature") +
    ylab("density") +
    ggtitle("Temperature kernel density estimate", subtitle = "Percentiles 0.01 and 0.99")
}

#' @export
get_thermal_envelope <- function(scientificname = NULL, aphiaid = NULL, taxonkey = NULL) {
  taxonomy <- insistent_resolve_taxonomy(scientificname, aphiaid, taxonkey)
  stopifnot(is.numeric(taxonomy$aphiaid) | is.numeric(taxonomy$taxonkey))

  sf_use_s2(FALSE)

  points_per_m2 <- 1 / 300000000
  max_points <- 10000
  kde_bandwidth <- 0.8
  min_area <- 2000000000
  min_temperatures <- 3

  if (!is.null(taxonomy$aphiaid)) {
    obis_dist <- get_obis_dist_possibly(aphiaid = taxonomy$aphiaid, res = 4)
  } else {
    obis_dist <- NULL
  }
  if (!is.null(taxonomy$taxonkey)) {
    gbif_dist <- get_gbif_dist_possibly(taxonkey = taxonomy$taxonkey, zoom = 2, square_size = 8)
  } else {
    gbif_dist <- NULL
  }
  if (is.null(obis_dist) & is.null(gbif_dist)) {
    return(NULL)
  }
  # TODO: subtract environmental layer from dist
  dist <- aggregate(vect(c(obis_dist, gbif_dist)))
  if (nrow(dist) == 0) {
    return(NULL)
  }

  surface_area <- expanse(dist, "m")
  if (surface_area < min_area) {
    return(NULL)
  }
  n_points <- min(max_points, round(surface_area * points_per_m2))
  # TODO: spatSample not returning expected number of points
  points <- st_sample(st_as_sf(dist), n_points)

  # TODO: convert to terra
  env <- download_layers(dataset_id = "thetao_baseline_2000_2019_depthsurf", variables = "thetao_mean", constraints = list(time = c("2010-01-01", "2010-01-01")), fmt = "raster", directory = tempdir(), verbose = T) %>%
    st_as_stars()
  temperatures <- st_extract(env, points) %>% as.data.frame() %>% pull(1) %>% na.omit() %>% as.numeric()
  if (length(temperatures) < min_temperatures) {
    return(NULL)
  }
  kd <- ks::kde(temperatures, h = kde_bandwidth)
  percentiles <- qkde(c(0.01, 0.99), kd)

  masked <- env
  masked[env >= percentiles[1] & env <= percentiles[2]] <- 1
  masked[env < percentiles[1] | env > percentiles[2]] <- NA

  return(list(
    envelope = st_as_sf(masked, as_points = FALSE, merge = TRUE),
    kd = kd,
    points = vect(points),
    percentiles = percentiles,
    temperatures = temperatures
  ))

}
