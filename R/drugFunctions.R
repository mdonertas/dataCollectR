#' Download drug gene interactions from DGIdb
#'
#' @param intfile link to the interactionfile from DGIdb
#'
#' @return returns a dataframe of 7 columns
#' @export
#'
download_drug_interactions <- function(intfile = "https://dgidb.org/data/monthly_tsvs/2020-Nov/interactions.tsv") {
  interactions <- read_tsv(intfile, n_max = 100) %>%
    select(-2, -6, -7, -8) %>%
    set_names(c("Gene", "EntrezID", "interactionSource", "interactionType",
                "ChEMBLID", "interactionScore", "PMIDs")) %>%
    mutate(ChEMBLID = gsub("chembl:", "", ChEMBLID))
  return(interactions)
}

#' Get molecule information for a ChEMBL ID
#'
#' @param ChEMBLID ChEMBLID
#'
#' @return returns a data.frame with 19 columns: "ChEMBLID", "name", "type",
#' "therapeutic", "first_approval", "indication_class", "max_phase",
#' "atc_classifications", "black_box_warning", "oral", "topical", "natural",
#' "availability", "chirality", "inorganic", "smiles", "inchi", "inchi_key",
#' "synonyms"
#'
get_info4ChEMBLID <- function(ChEMBLID) {
  nm <- URLencode(toupper(ChEMBLID), reserved = T)
  dat <- getURL(paste("https://www.ebi.ac.uk/chembl/api/data/molecule/",
                      nm, ".json",
                      sep = ""
  ))
  if(dat==''){
    dat = NA
  } else{
    dat <- fromJSON(dat)
    dat <- data.frame(
      ChEMBLID = ChEMBLID,
      name = dat$pref_name,
      type = dat$molecule_type,
      therapeutic = dat$therapeutic_flag,
      first_approval = dat$first_approval,
      indication_class = paste(sort(unique(dat$indication_class)),
                               collapse = ", "),
      max_phase = dat$max_phase,
      atc_classifications = paste(sort(unique(dat$atc_classifications)),
                                  collapse = ", "),
      black_box_warning = dat$black_box_warning,
      oral = dat$oral,
      topical = dat$topical,
      natural = dat$natural_product,
      availability = dat$availability_type,
      chirality = dat$chirality,
      inorganic = dat$inorganic_flag,
      smiles = dat$molecule_structures$canonical_smiles,
      inchi = dat$molecule_structures$standard_inchi,
      inchi_key = dat$molecule_structures$standard_inchi_key,
      synonyms = paste(sort(unique(c(dat$molecule_synonyms$molecule_synonym,
                                     dat$molecule_synonyms$synonyms))),
                       collapse = ", "))
  }
  return(dat)
}

#' Get molecule information for a vector of ChEMBL IDs
#'
#' @param ChEMBLIDs a character vector of ChEMBL IDs
#'
#' @return returns a list with three elements. First: `ChEMBL_info`, a data
#' frame with 19 columns ( "ChEMBLID", "name", "type", "therapeutic",
#' "first_approval", "indication_class", "max_phase", "atc_classifications",
#' "black_box_warning", "oral", "topical", "natural", "availability",
#' "chirality", "inorganic", "smiles", "inchi", "inchi_key", "synonyms"),
#' Second: `pkginfo`, a data frame with the package version information,
#' Third: `accessdate`: date of data accession (UTC time zone).
#' @export
#'
#' @examples
#' get_info4ChEMBLIDs(c('CHEMBL413','CHEMBL1431'))
#'
get_info4ChEMBLIDs <- function(ChEMBLIDs){
  local_locale(c("LC_TIME" = "C"))
  local_timezone("UTC")
  time <- format(Sys.time(), "%Y_%m_%d_%H_%M")
  info <- lapply(ChEMBLIDs, get_info4ChEMBLID)
  info <- melt(info, id.vars = colnames(info[[1]])) %>%
    select(-L1)
  pkginfo <- data.frame(package = c("RCurl", "jsonlite", "tidyverse",
                                    "reshape2","utils")) %>%
    rowwise() %>%
    mutate(version = packageVersion(package))
  chembl <- list(ChEMBL_info = info, pkginfo = pkginfo, accessdate = time)
  return(chembl)
}

