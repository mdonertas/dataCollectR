#' Get gene name, Ensembl IDs, Entrez IDs, Uniprot IDs, and gene biotype for all
#' genes
#'
#' @param organism name of the organism: first letter of the genus immediately
#' followed by the specific name
#'
#' @return returns a list with three elements. First: `geneIDs`, a data frame
#' with 9 columns (external_gene_name, ensembl_gene_id, entrezgene_id,
#' uniprotswissprot, uniprot_gn_id, external_synonym, gene_biotype, and
#' description, organism). Second: `pkginfo`, a data frame with the package
#' version information, Third: `accessdate`: date of data accession (UTC time 
#' zone).
#'
#' @importFrom utils packageVersion
#' @export
get_all_geneIDs <- function(organism = "hsapiens") {
  local_locale(c("LC_TIME" = "C"))
  local_timezone("UTC")
  time <- format(Sys.time(), "%Y_%m_%d_%H_%M")
  ensembl <- biomaRt::useEnsembl(biomart = "genes")
  dataset <- paste(tolower(organism), "_gene_ensembl", sep = "")
  ensembl <- biomaRt::useDataset(dataset = dataset, mart = ensembl)
  attrlist = intersect(c(
    "external_gene_name",
    "ensembl_gene_id",
    "entrezgene_id",
    "uniprotswissprot",
    "uniprot_gn_id",
    "uniprotsptrembl",
    # "external_synonym",
    "gene_biotype",
    "description"
  ), biomaRt::listAttributes(ensembl)$name)
  idmap <- biomaRt::getBM(
    attributes = attrlist,
    mart = ensembl
  ) %>%
    mutate(dataset = organism)
  pkginfo <- data.frame(package = c("biomaRt", "tidyverse")) %>%
    rowwise() %>%
    mutate(version = packageVersion(package))
  idmap <- list(geneIDs = idmap, pkginfo = pkginfo, accessdate = time)
  return(idmap)
}

#' Chemical Name to PubChem CID
#'
#' @param chemname Chemical name
#' @importFrom utils URLencode
#' @importFrom RCurl getURL
#' @importFrom jsonlite fromJSON
#' @return returns pubchem CID
map_chemName2CID <- function(chemname) {
  nm <- URLencode(chemname, reserved = T)
  name2cid <- getURL(paste("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/",
    nm, "/cids/JSON",
    sep = ""
  ))
  name2cid <- fromJSON(name2cid)
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
#' @importFrom reshape2 melt
#' @export
#'
#' @examples
#' map_chemName2CIDs(c("metformin", "rapamycin"))
map_chemName2CIDs <- function(chem_names) {
  local_locale(c("LC_TIME" = "C"))
  local_timezone("UTC")
  time <- format(Sys.time(), "%Y_%m_%d_%H_%M")
  cids <- lapply(chem_names, map_chemName2CID)
  names(cids) <- chem_names
  cids <- melt(cids) %>%
    rename(CID = value, chemName = L1) %>%
    full_join(data.frame(chemName = as.character(chem_names))) %>%
    select(chemName, CID)
  pkginfo <- data.frame(package = c("RCurl", "jsonlite", "tidyverse",
                                    "reshape2")) %>%
    rowwise() %>%
    mutate(version = packageVersion(package))
  cids <- list(CIDs = cids, pkginfo = pkginfo, accessdate = time)
  return(cids)
}

#' PubChem CID to ChEMBL ID
#'
#' @param cid PubChem CID
#'
#' @return returns ChEMBLID
map_CID2ChEMBLID <- function(cid) {
  nm <- URLencode(cid, reserved = T)
  cid2chembl <- getURL(paste("https://www.ebi.ac.uk/unichem/rest/src_compound_id/",
    nm, "/22/1",
    sep = ""
  ))
  cid2chembl <- fromJSON(cid2chembl)
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
#' @export
#'
#' @examples
#' map_CID2ChEMBLIDs(c("4091", "5284616"))
map_CID2ChEMBLIDs <- function(cids) {
  local_locale(c("LC_TIME" = "C"))
  local_timezone("UTC")
  time <- format(Sys.time(), "%Y_%m_%d_%H_%M")
  chembl <- lapply(cids, map_CID2ChEMBLID)
  names(chembl) <- cids
  chembl <- melt(chembl) %>%
    rename(ChEMBLID = value, CID = L1) %>%
    full_join(data.frame(CID = as.character(cids))) %>%
    select(CID, ChEMBLID)
  pkginfo <- data.frame(package = c("RCurl", "jsonlite", "tidyverse",
                                    "reshape2","utils")) %>%
    rowwise() %>%
    mutate(version = packageVersion(package))
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
#' @export
#'
#' @examples
#' map_chemName2ChEMBLIDs(c("metformin", "rapamycin"))
map_chemName2ChEMBLIDs <- function(chem_names) {
  local_locale(c("LC_TIME" = "C"))
  local_timezone("UTC")
  time <- format(Sys.time(), "%Y_%m_%d_%H_%M")
  cids <- lapply(chem_names, map_chemName2CID)
  names(cids) <- chem_names
  cids <- melt(cids) %>%
    rename(CID = value, chemName = L1) %>%
    full_join(data.frame(chemName = as.character(chem_names))) %>%
    select(chemName, CID)
  chembl <- lapply(unique(cids$CID), map_CID2ChEMBLID)
  names(chembl) <- unique(cids$CID)
  chembl <- melt(chembl) %>%
    rename(ChEMBLID = value, CID = L1) %>%
    full_join(data.frame(CID = as.character(cids))) %>%
    select(CID, ChEMBLID)
  chembl <- left_join(cids, chembl)
  pkginfo <- data.frame(package = c("RCurl", "jsonlite", "tidyverse",
                                    "reshape2","utils")) %>%
    rowwise() %>%
    mutate(version = packageVersion(package))
  chembl <- list(ChEMBLIDs = chembl, pkginfo = pkginfo, accessdate = time)
  return(chembl)
}

