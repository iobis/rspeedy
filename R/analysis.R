#' @export
calculate_plausibility <- function(dist, buffer = 1, obis = TRUE, gbif = TRUE, worms = TRUE) {
  sf_use_s2(FALSE)

  ranges <- data.frame(geometry = st_sfc(crs = 4326))
  if (obis & nrow(dist$obis) > 0) {
    ranges <- bind_rows(
      ranges,
      data.frame(geometry = st_union(dist$obis))
    )
  }
  if (gbif & nrow(dist$gbif) > 0) {
    ranges <- bind_rows(
      ranges,
      data.frame(geometry = st_union(dist$gbif))
    )
  }
  if (worms & nrow(dist$worms) > 0) {
    ranges <- bind_rows(
      ranges,
      data.frame(geometry = st_union(dist$worms))
    )
  }

  geom_range <- ranges %>%
    st_as_sf() %>%
    st_union() %>%
    st_buffer(buffer) %>%
    st_crop(st_as_sfc(st_bbox(c(xmin = -180, xmax = 180, ymin = -90, ymax = 90), crs = 4326)))

  geom_env <- dist$envelope$envelope

  # clip geom_env with geom_range
  if (length(geom_range) > 0 & length(geom_env) > 0) {
    geom_env <- geom_env %>%
      st_difference(geom_range)
  }
  if (length(geom_env) > 0) {
    geom_env <- geom_env %>%
      st_union()
  }

  # combine

  plausibility <- data.frame(plausibility = numeric(), geometry = st_sfc(crs = 4326))
  if (length(geom_range) > 0) {
    plausibility <- bind_rows(plausibility, data.frame(plausibility = 1, geometry = geom_range))
  }
  if (length(geom_env) > 0) {
    plausibility <- bind_rows(plausibility, data.frame(plausibility = 0.5, geometry = geom_env))
  }

  plausibility <- plausibility %>%
    st_sf()

  return(plausibility)
}

#' @export
plot_plausibility <- function(pl) {
  ggplot() +
    geom_sf(data = st_as_sfc(st_bbox(c(xmin = -180, ymin = -85, xmax = 180, ymax = 85), crs = 4326)), fill = "#E2E2E0", color = NA) +
    geom_sf(data = landr::get_land_polygons(simplified = 20000), fill = "#fafafc", color = NA) +
    geom_sf(data = pl, aes(fill = plausibility), color = NA) +
    scale_color_brewer(palette = "Spectral") +
    theme_void() +
    theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5, vjust = 1)) +
    coord_sf(expand = 0.1)
}
