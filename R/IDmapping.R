#' Get gene name, Ensembl IDs, Entrez IDs, Uniprot IDs, and gene biotype for all
#' genes
#'
#' @param organism name of the organism: first letter of the genus immediately
#' followed by the specific name
#'
#' @return returns a list with three elements. First: `geneIDs`, a data frame
#' with 5 columns (Gene, EnsemblID, EntrezID, UniprotID, geneType), Second:
#' `pkginfo`, a data frame with the package version information, Third:
#' `accessdate`: date of data accession (UTC time zone).
#'
#' @export
get_all_geneIDs <- function(organism = "hsapiens") {
  withr::local_locale(c("LC_TIME" = "C"))
  withr::local_timezone("UTC")
  time <- format(Sys.time(), "%Y_%m_%d_%H_%M")
  suppressMessages(require(tidyverse))
  ensembl <- biomaRt::useEnsembl(biomart = "genes")
  dataset <- paste(tolower(organism), "_gene_ensembl", sep = "")
  ensembl <- biomaRt::useDataset(dataset = dataset, mart = ensembl)
  idmap <- biomaRt::getBM(
    attributes = c(
      "external_gene_name",
      "ensembl_gene_id",
      "entrezgene_id",
      "uniprotswissprot",
      "gene_biotype"
    ),
    mart = ensembl
  ) %>%
    set_names(c("Gene", "EnsemblID", "EntrezID", "UniprotID", "geneType"))
  pkginfo <- data.frame(package = c("biomaRt", "tidyverse")) %>%
    rowwise() %>%
    mutate(version = utils::packageVersion(package))
  idmap <- list(geneIDs = idmap, pkginfo = pkginfo, accessdate = time)
  return(idmap)
}

#' Chemical Name to PubChem CID
#'
#' @param chemname Chemical name
#'
#' @return returns pubchem CID
map_chemName2CID <- function(chemname) {
  nm <- utils::URLencode(chemname, reserved = T)
  name2cid <- RCurl::getURL(paste("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/",
    nm, "/cids/JSON",
    sep = ""
  ))
  name2cid <- jsonlite::fromJSON(name2cid)
  cid <- as.character(name2cid$IdentifierList$CID)
  return(cid)
}

#' Map chemicals to PubChem CIDs
#'
#' @param chem_names a character vector including names of the drugs / chemicals
#'
#' @return returns a list with three elements. First: `CIDs`, a data frame with
#' 2 columns (chemName and CID), Second: `pkginfo`, a data frame with the
#' package version information, Third: `accessdate`: date of data accession
#' (UTC time zone).
#' @import tidyverse
#' @export
#'
#' @examples
#' map_chemName2CIDs(c("metformin", "rapamycin"))
map_chemName2CIDs <- function(chem_names) {
  withr::local_locale(c("LC_TIME" = "C"))
  withr::local_timezone("UTC")
  time <- format(Sys.time(), "%Y_%m_%d_%H_%M")
  suppressMessages(require(tidyverse))
  cids <- lapply(chem_names, map_chemName2CID)
  names(cids) <- chem_names
  cids <- reshape2::melt(cids) %>%
    rename(CID = value, chemName = L1) %>%
    full_join(data.frame(chemName = as.character(chem_names))) %>%
    select(chemName, CID)
  pkginfo <- data.frame(package = c("RCurl", "jsonlite", "tidyverse",
                                    "reshape2")) %>%
    rowwise() %>%
    mutate(version = utils::packageVersion(package))
  cids <- list(CIDs = cids, pkginfo = pkginfo, accessdate = time)
  return(cids)
}

#' PubChem CID to ChEMBL ID
#'
#' @param cid PubChem CID
#'
#' @return returns ChEMBLID
map_CID2ChEMBLID <- function(cid) {
  nm <- utils::URLencode(cid, reserved = T)
  cid2chembl <- RCurl::getURL(paste("https://www.ebi.ac.uk/unichem/rest/src_compound_id/",
    nm, "/22/1",
    sep = ""
  ))
  cid2chembl <- jsonlite::fromJSON(cid2chembl)
  chembl <- as.character(cid2chembl$src_compound_id)
  return(chembl)
}

#' Map PubChem CIDs to ChEMBL IDs
#'
#' @param cids a character vector including PubChem CIDs
#'
#' @return returns a list with three elements. First: `ChEMBLIDs`, a data frame
#' with 2 columns (CID and ChEMBLID), Second: `pkginfo`, a data frame with the
#' package version information, Third: `accessdate`: date of data accession
#' (UTC time zone).
#' @import tidyverse
#' @export
#'
#' @examples
#' map_CID2ChEMBLIDs(c("4091", "5284616"))
map_CID2ChEMBLIDs <- function(cids) {
  withr::local_locale(c("LC_TIME" = "C"))
  withr::local_timezone("UTC")
  time <- format(Sys.time(), "%Y_%m_%d_%H_%M")
  suppressMessages(require(tidyverse))
  chembl <- lapply(cids, map_CID2ChEMBLID)
  names(chembl) <- cids
  chembl <- reshape2::melt(chembl) %>%
    rename(ChEMBLID = value, CID = L1) %>%
    full_join(data.frame(CID = as.character(cids))) %>%
    select(CID, ChEMBLID)
  pkginfo <- data.frame(package = c("RCurl", "jsonlite", "tidyverse",
                                    "reshape2","utils")) %>%
    rowwise() %>%
    mutate(version = utils::packageVersion(package))
  chembl <- list(ChEMBLIDs = chembl, pkginfo = pkginfo, accessdate = time)
  return(chembl)
}

#' Map chemicals to ChEMBL IDs
#'
#' @param chem_names a character vector including drug / chemical names
#'
#' @return returns a list with three elements. First: `ChEMBLIDs`, a data frame
#' with 3 columns (chemName, CID and ChEMBLID), Second: `pkginfo`, a data frame
#' with the package version information, Third: `accessdate`: date of data
#' accession (UTC time zone).
#' @import tidyverse
#' @export
#'
#' @examples
#' map_chemName2ChEMBLIDs(c("metformin", "rapamycin"))
map_chemName2ChEMBLIDs <- function(chem_names) {
  withr::local_locale(c("LC_TIME" = "C"))
  withr::local_timezone("UTC")
  time <- format(Sys.time(), "%Y_%m_%d_%H_%M")
  suppressMessages(require(tidyverse))
  cids <- lapply(chem_names, map_chemName2CID)
  names(cids) <- chem_names
  cids <- reshape2::melt(cids) %>%
    rename(CID = value, chemName = L1) %>%
    full_join(data.frame(chemName = as.character(chem_names))) %>%
    select(chemName, CID)
  chembl <- lapply(unique(cids$CID), map_CID2ChEMBLID)
  names(chembl) <- unique(cids$CID)
  chembl <- reshape2::melt(chembl) %>%
    rename(ChEMBLID = value, CID = L1) %>%
    full_join(data.frame(CID = as.character(cids))) %>%
    select(CID, ChEMBLID)
  chembl <- left_join(cids, chembl)
  pkginfo <- data.frame(package = c("RCurl", "jsonlite", "tidyverse",
                                    "reshape2","utils")) %>%
    rowwise() %>%
    mutate(version = utils::packageVersion(package))
  chembl <- list(ChEMBLIDs = chembl, pkginfo = pkginfo, accessdate = time)
  return(chembl)
}

#' Get molecule information for a ChEMBL ID
#'
#' @param ChEMBLID
#'
#' @return returns a data.frame with 19 columns: "ChEMBLID", "name", "type",
#' "therapeutic", "first_approval", "indication_class", "max_phase",
#' "atc_classifications", "black_box_warning", "oral", "topical", "natural",
#' "availability", "chirality", "inorganic", "smiles", "inchi", "inchi_key",
#' "synonyms"
#'
get_info4ChEMBLID <- function(ChEMBLID) {
  nm <- utils::URLencode(toupper(ChEMBLID), reserved = T)
  dat <- RCurl::getURL(paste("https://www.ebi.ac.uk/chembl/api/data/molecule/",
    nm, ".json",
    sep = ""
  ))
  if(dat==''){
    dat = NA
  } else{
    dat <- jsonlite::fromJSON(dat)
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
  withr::local_locale(c("LC_TIME" = "C"))
  withr::local_timezone("UTC")
  time <- format(Sys.time(), "%Y_%m_%d_%H_%M")
  suppressMessages(require(tidyverse))
  info <- lapply(ChEMBLIDs, get_info4ChEMBLID)
  info <- reshape2::melt(info, id.vars = colnames(info[[1]])) %>%
    select(-L1)
  pkginfo <- data.frame(package = c("RCurl", "jsonlite", "tidyverse",
                                    "reshape2","utils")) %>%
    rowwise() %>%
    mutate(version = utils::packageVersion(package))
  chembl <- list(ChEMBL_info = info, pkginfo = pkginfo, accessdate = time)
  return(chembl)
}

