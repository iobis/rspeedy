library(sdmpredictors)
library(stars)
library(ks)
library(rgbif)
library(mapview)
library(tidyr)
library(ggplot2)
library(eks)

#' @export
view_envelope <- function(envelope, points = FALSE) {
  if (points) {
    mapview(envelope$envelope) + mapview(envelope$points)
  } else {
    mapview(envelope$envelope)
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
  taxonomy <- resolve_taxonomy(scientificname, aphiaid, taxonkey)
  stopifnot(is.numeric(taxonomy$aphiaid) & is.numeric(taxonomy$taxonkey))

  variable <- "BO22_tempmean_ss"
  points_per_m2 <- 1 / 300000000
  max_points <- 10000
  kde_bandwidth <- 0.8

  obis_dist <- get_obis_dist(aphiaid = taxonomy$aphiaid, res = 4)
  gbif_dist <- get_gbif_dist(taxonkey = taxonomy$taxonkey, zoom = 2, square_size = 8)

  dist <- bind_rows(obis_dist, gbif_dist) %>%
    st_union() %>%
    st_make_valid()

  n <- min(max_points, round(st_area(dist) * points_per_m2))
  points <- st_sample(dist, n)

  env <- load_layers(variable) %>%
    st_as_stars()

  temperatures <- st_extract(env, points) %>% pull(variable) %>% na.omit() %>% as.vector()
  kd <- ks::kde(temperatures, h = kde_bandwidth)
  percentiles <- qkde(c(0.01, 0.99), kd)

  masked <- env
  masked[env >= percentiles[1] & env <= percentiles[2]] <- 1
  masked[env < percentiles[1] | env > percentiles[2]] <- NA

  return(list(
    envelope = st_as_sf(masked, as_points = FALSE, merge = TRUE),
    kd = kd,
    points = points,
    percentiles = percentiles,
    temperatures = temperatures
  ))

}
