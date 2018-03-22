read_matrix_data <- function(seqtab_rds) {
    mat <- readRDS(seqtab_rds)
    dimnames(mat)[[2]] <- map_chr(dimnames(mat)[[2]], digest::digest)
    dt <- as.data.frame(mat)
    dt$stations <- dimnames(mat)[[1]]
    res <- tidyr::gather(dt, seq_hash, n_seqs, -stations) %>%
        dplyr::filter(n_seqs > 0) %>%
        dplyr::mutate(stations = gsub("(.+)_(.+)", "\\2", stations))
}

plot_sequence_abundance <- function(seq_abdce) {
    seq_abdce %>%
        ggplot() +
        geom_col(aes(x = seq_hash, y = n_seqs, fill = stations)) +
        facet_wrap(~ stations)
}


if (FALSE) {
    mat <- readRDS("data/seqtab_final.rds")
    bcs_d <- dist(mat)
    library(MASS)
    bcs_nmds <- isoMDS(bcs_d)
    lbls <- gsub("(.+)_(.+)", "\\2", rownames(bcs_nmds$points))
    plot(bcs_nmds$points)
    text(bcs_nmds$points, label = lbls)

}

plot_species_richness <- function(seq_abdce, sp_assign, project_file = "data/station_project.csv") {

    sp_assign_d <- sp_assign %>%
        dplyr::filter(query_seqid >=  90, !is.na(phylum)) %>%
        dplyr::group_by(query_seqid) %>%
        dplyr::arrange(desc(pident)) %>%
        dplyr::slice(1) %>%
        dplyr::ungroup()

    res <- seq_abdce %>%
        dplyr::left_join(sp_assign_d, by = c("seq_hash" = "query_seqid"))

    proj <- readr::read_csv(project_file, col_types = "cc")

    res <- res %>%
        dplyr::left_join(proj, by = "stations")

    res %>%
        dplyr::group_by(project, stations, phylum, class, order) %>%
        dplyr::summarize(
                   n_seqs = sum(n_seqs),
                   ) %>%
        dplyr::group_by(stations) %>%
        dplyr::mutate(p_seqs = n_seqs/sum(n_seqs)) %>%
        dplyr::ungroup() %>%
        ggplot() +
        geom_col(aes(x = interaction(stations, project), y = p_seqs, fill = interaction(order, class, phylum)),
                 position = "stack") +
        theme(axis.text.x = element_text(angle = 90, hjust = 0))
}
