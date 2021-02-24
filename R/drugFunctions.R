#' Download drug gene interactions from DGIdb
#'
#' @param intfile link to the interactionfile from DGIdb
#'
#' @return returns a dataframe of 7 columns
#' @export
#'
#' @examples
download_drug_interactions <- function(intfile = "https://dgidb.org/data/monthly_tsvs/2020-Nov/interactions.tsv") {
  suppressMessages(require(tidyverse))
  interactions <- readr::read_tsv(intfile, n_max = 100) %>%
    select(-2, -6, -7, -8) %>%
    set_names(c("Gene", "EntrezID", "interactionSource", "interactionType",
                "ChEMBLID", "interactionScore", "PMIDs")) %>%
    mutate(ChEMBLID = gsub("chembl:", "", ChEMBLID))
  return(interactions)
}

