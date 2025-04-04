#' Get all human orthologs for an organism
#'
#' @param organism name of the organism: first letter of the genus immediately
#' followed by the specific name
#'
#' @return returns a list with three elements. First: `orthologs`, a data frame
#' with 9 columns (EnsemblID, humanEnsemblID, Type, humanInQuery, queryInHuman,
#' GOC, WGA, confidence, dataset), Second: `pkginfo`, a data frame with the
#' package version information, Third: `accessdate`: date of data accession
#' (UTC time zone).
#' @export
#'
get_all_orthologs <- function(organism = "mmusculus", target = "hsapiens") {
  local_locale(c("LC_TIME" = "C"))
  local_timezone("UTC")
  time <- format(Sys.time(), "%Y_%m_%d_%H_%M")
  ensembl <- biomaRt::useEnsembl(biomart = "genes")
  dataset <- paste(tolower(organism), "_gene_ensembl", sep = "")
  ensembl <- biomaRt::useDataset(dataset = dataset, mart = ensembl)
  orthogs <- paste(target, c("homolog_ensembl_gene",
                             "homolog_orthology_type",
                             "homolog_perc_id",
                             "homolog_perc_id_r1",
                             "homolog_goc_score",
                             "homolog_wga_coverage",
                             "homolog_orthology_confidence"), sep = "_")
  idmap <- biomaRt::getBM(
    attributes = c(
      "ensembl_gene_id",
      orthogs
    ),
    mart = ensembl
  ) %>%
    mutate(dataset = organism)
  pkginfo <- data.frame(package = c("biomaRt", "tidyverse")) %>%
    rowwise() %>%
    mutate(version = packageVersion(package))
  idmap <- list(orthologs = idmap, pkginfo = pkginfo, accessdate = time)
  return(idmap)
}

#' Get all IDs and human orthologs for an organism
#'
#' @param organism name of the organism: first letter of the genus immediately
#' followed by the specific name
#'
#' @return returns a list with three elements. First: `orthologs`, a data frame
#' with 13 columns (Gene, EnsemblID, EntrezID, UniprotID, geneType, dataset,
#' humanEnsemblID, Type, humanInQuery, queryInHuman, GOC, WGA, confidence),
#' Second: `pkginfo`, a data frame with the package version information,
#' Third: `accessdate`: date of data accession (UTC time zone).
#' @export
#'
tidy_orthologs <- function(organism = 'mmusculus', target = 'hsapiens') {
  local_locale(c("LC_TIME" = "C"))
  local_timezone("UTC")
  time <- format(Sys.time(), "%Y_%m_%d_%H_%M")
  idmap <- get_all_geneIDs(organism = organism)
  ortholog <- get_all_orthologs(organism = organism, target = target)
  pkginfo <- data.frame(package = c("biomaRt", "tidyverse")) %>%
    rowwise() %>%
    mutate(version = packageVersion(package))
  ortholog <- left_join(idmap$geneIDs,ortholog$orthologs)
  ortholog <- list(orthologs = ortholog, pkginfo = pkginfo, accessdate = time)
  return(ortholog)
}


