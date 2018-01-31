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
