get_blast_data <- function(path, pattern = "blasted_.+\\.txt",
                           database = c("nt", "midori")) {

    lst_files <- list.files(path = path, pattern = pattern,
                            full.names = TRUE)

    res <- purrr::map_dfr(lst_files, readr::read_delim,
                          delim = "\t", col_names = FALSE,
                          .id = "identity", col_types = c("ccdiiiiiiidic"))

    names(res)[2:ncol(res)] <- c("query_seqid", "all_subject_seqid",
                                 "pident", "length", "mismatch", "gapopen",
                                 "qstart", "qend", "start_alignment",
                                 "end_alignment", "evalue", "bitscore",
                                 "taxonomy")
    id_level <- gsub(".+_i([0-9]+)\\.txt", "\\1", lst_files)
    res$identity <- id_level[as.numeric(res$identity)]

    ## the midori database has the taxonomy ID and string in the taxononmy
    ## column but the nt database doesn't
    if (identical(database, "midori")) {
        res %>%
            tidyr::separate(taxonomy, into = c("taxid", "taxstring"),
                            sep = "\\s{3}") %>%
            dplyr::mutate(classif = purrr::map(taxid, ~ ncbi_classification_from_genbank(.))) %>%
            return()

    } else if (identical(database, "nt")) {
        tt <- res %>%
            dplyr::mutate(genbank_id = stringr::str_extract(all_subject_seqid, "[A-Z]{1,2}[0-9]{5,6}\\.[1-9]{1}"),
                          genbank_id = gsub("\\.[1-9]{1}$", "", genbank_id)) %>%
            dplyr::mutate(classif = purrr::map(genbank_id, ~ ncbi_classification_from_genbank(.))) %>%
            return()

    } else stop("invalid")
}


add_classification <- function(blast_res) {
    blast_res %>%
        dplyr::mutate(tmp_classif = purrr::map(classif, function(x, ...) {

            if (!inherits(x, "data.frame"))
                return(
                    data_frame(
                        phylum = NA,
                        class = NA,
                        order = NA,
                        family = NA,
                        genus = NA,
                        species = NA
                    ))
            x %>%
                dplyr::filter(rank %in% c("phylum", "class", "order", "family", "genus", "species")) %>%
                dplyr::select(-id) %>%
                tidyr::spread(rank, name)
        })) %>%
        dplyr::select(-classif) %>%
        tidyr::unnest()
}


fetch_hook_ncbi_id <- function(key, namespace) {
    message("ncbi sequence ID: ", key)
    taxize::genbank2uid(key)
}

store_ncbi_id <- function(path = "data-store/ncbi-id") {
    storr::storr_external(storr::driver_rds(path),
                          fetch_hook_ncbi_id)
}

fetch_hook_ncbi_classification <- function(key, namespace) {
    message("ncbi taxonomy ID: ", key)
    taxize::classification(key, db = "ncbi")
}

store_ncbi_classification <- function(path = "data-store/ncbi-classification") {
    storr::storr_external(storr::driver_rds(path),
                          fetch_hook_ncbi_classification)
}

ncbi_classification_from_genbank <- function(id) {
    message("ID: ", id)
    ncbi_id    <- store_ncbi_id()$get(id)
    if (is.na(ncbi_id)) return(NA)
    classif_id <- store_ncbi_classification()$get(ncbi_id[[1]])
    select_classif <- classif_id[[1]][classif_id[[1]]$rank %in%
                                      c("superkingdom",
                                        "kingdom",
                                        "phylum",
                                        "class",
                                        "order",
                                        "genus",
                                        "species"), ]
    select_classif
}
