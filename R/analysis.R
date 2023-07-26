#' @export
calculate_plausibility <- function(dist, buffer = 1) {
  sf_use_s2(FALSE)

  geom_range <- c(st_geometry(dist$obis), st_geometry(dist$gbif), st_geometry(dist$worms)) %>%
    st_union() %>%
    st_buffer(1) %>%
    st_crop(st_as_sfc(st_bbox(c(xmin = -180, xmax = 180, ymin = -90, ymax = 90), crs = 4326)))
  geom_env <- dist$envelope$envelope %>%
    st_difference(geom_range) %>%
    st_union()

  plausibility <- data.frame(
    plausibility = c(2, 1),
    geometry = c(geom_range, geom_env)
  ) %>%
    st_as_sf()

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
