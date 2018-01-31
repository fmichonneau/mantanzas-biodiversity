## seqtab <- readRDS("data/seqtab_final.rds")

## stations <- dimnames(seqtab)[[1]]
## stations <- gsub(".+(MCPK-.+$)", "\\1", stations[grepl("MCPK", stations)])

extract_metabarcoding_stations <- function(station_data, sample_metacommunity,
                                           extraction_metacommunity, pcr_metacommunity,
                                           library_id = "BCS-003") {

    pcr_metacommunity %>%
        dplyr::filter(library_id == library_id) %>%
        dplyr::left_join(dplyr::select(extraction_metacommunity, .data$extraction_id, .data$metacommunity_sample_id),
                         by = "extraction_id") %>%
        dplyr::left_join(sample_metacommunity, by = "metacommunity_sample_id") %>%
        dplyr::left_join(station_data, by = "station_number") %>%
        dplyr::select(.data$extraction_id, .data$extraction_id_unique,
                      .data$metacommunity_sample_id, .data$station_number,
                      .data$fraction, .data$year, .data$month, .data$day,
                      starts_with("longitude"), starts_with("latitude"))
}


map_coord_limits <- function() {
    ## limits of plotting area for area sampled in Matanzas river
    list(
        xlim = c(-81.25, -81.21),
        ylim = c(29.685, 29.725)
    )
}

fortify_shp_file <- function(path, layer) {
    shp <- rgdal::readOGR(path, layer)
    ggplot2::fortify(shp)
}


map_metabarcoding_stations <- function(station_data, map_data) {
    ## limit to Whitney/Matanzas samples
    station_data <- station_data %>%
        dplyr::mutate(date = lubridate::ymd(paste(year, month, day, sep = "-"))) %>%
        dplyr::filter(longitude_start > -81.5 & latitude_start > 29.5)

    ggplot() +
        geom_polygon(data = map_data, aes(x= long, y = lat, group = group), fill="gray20") +
        geom_point(data = station_data,
                   aes(x = longitude_start, y = latitude_start, color = as.factor(date)),
                   position = position_jitter(.0005, .0005), inherit.aes = FALSE) +
        coord_quickmap(xlim = map_coord_limits()$xlim,
                       ylim = map_coord_limits()$ylim) +
        xlab("Longitude") + ylab("Latitude") +
        scale_color_discrete(name = "Collection date") +
        theme_bw()
}

map_florida <- function(map_data) {
    ## extract & convert limits for red square on Florida map
    rect_limits <- data_frame(
        xmin = min(map_coord_limits()$xlim) - .05,
        xmax = max(map_coord_limits()$xlim) + .05,
        ymin = min(map_coord_limits()$ylim) - .05,
        ymax = max(map_coord_limits()$ylim) + .05

    )

    ggplot() +
        geom_polygon(data = map_data, aes(x = long, y = lat, group = group), fill = "gray20") +
        geom_rect(data = rect_limits, aes(xmin = xmin, xmax = xmax,
                                          ymin = ymin, ymax = ymax),
                  color = "red", fill = NA) +
        coord_quickmap() +
        xlab("Longitude") + ylab("Latitude") +
        theme_bw()
}

draw_map_sampling <- function(fl_map, mat_map) {
    ## use patchwork magic to combine and annotate the plots
    fl_map + mat_map +
        plot_annotation(tag_levels = "A",
                        tag_suffix = ".")
}
