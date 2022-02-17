#' Download GenAge Human Data
#'
#' @param outprefix prefix for the name of output file, excluding the directory,
#' defaults to 'genage_human'
#' @param outdir directory to save the file, defaults to the current directory
#'
#' @return returns the file path
#' @import magrittr
#' @import dplyr
#' @importFrom withr local_locale local_timezone local_dir
#' @importFrom fs path dir_ls
#' @export
download_genage_human <- function(outprefix = "genage_human", outdir = "./") {
  local_locale(c("LC_TIME" = "C"))
  local_timezone("UTC")
  local_dir(outdir)
  time <- format(Sys.time(), "%Y%m%d%H%M%S")
  system("wget https://genomics.senescence.info/genes/human_genes.zip")
  system("unzip human_genes")
  system("rm human_genes.zip")
  outfile <- paste(outprefix, "_", time, ".csv", sep = "")
  fname <- grep("human", grep(".csv", dir_ls(), value = T), ignore.case = T,
                value = T)
  system(paste("mv", fname, outfile, sep = " "))
  system(paste("mv release.html ", outprefix, "_release_", time, ".html",
    sep = ""
  ))
  outfile <- path(outdir, outfile)
  return(outfile)
}

#' Tidy GenAge Human Data
#'
#' @param path2file path to genage human csv file
#'
#' @return returns a data frame with two columns: Gene and Why
#' @export
#'
tidy_genage_human <- function(path2file){
  dat <- read_csv(path2file) %>%
    select(2,6) %>%
    set_names(c('Gene','Why')) %>%
    unique()
  return(dat)
}

#' Download GenAge Model Organism Data
#'
#' @param outprefix prefix for the name of output file, excluding the directory,
#' defaults to 'genage_models'
#' @param outdir directory to save the file, defaults to the current directory
#'
#' @return returns the file path
#' @export
download_genage_model <- function(outprefix = "genage_models", outdir = "./") {
  withr::local_locale(c("LC_TIME" = "C"))
  withr::local_timezone("UTC")
  withr::local_dir(outdir)
  time <- format(Sys.time(), "%Y%m%d%H%M%S")
  system("wget https://genomics.senescence.info/genes/models_genes.zip")
  system("unzip models_genes")
  system("rm models_genes.zip")
  outfile <- paste(outprefix, "_", time, ".csv", sep = "")
  fname <- grep("model", grep(".csv", fs::dir_ls(), value = T), ignore.case = T,
                value = T)
  system(paste("mv", fname, outfile, sep = " "))
  system(paste("mv release.html ", outprefix, "_release_", time, ".html",
    sep = ""
  ))
  outfile <- fs::path(outdir, outfile)
  return(outfile)
}


#' Tidy GenAge Model Organism Data
#'
#' @param path2file path to genage model csv file
#'
#' @return returns a data frame with 6 columns: 'Gene','avgChange','Direction',
#' 'Type','organism','dataset'
#' @export
#'
tidy_genage_model <- function(path2file){
  dat <- read_csv(path2file)
  dat <- dat %>% select(organism) %>% unique() %>%
    separate(organism, into = c('genus','specific'), remove = F) %>%
    mutate(dataset = tolower(paste(substr(genus,1,1),specific,sep=''))) %>%
    select(organism, dataset) %>%
    right_join(dat) %>%
    rename(Gene = symbol,
           EntrezID = `entrez gene id`) %>%
    select(4, 7, 8, 9, 1, 2) %>%
    set_names(c('Gene','avgChange','Direction','Type','organism','dataset'))
  return(dat)
}

#' Download DrugAge Data
#'
#' @param outprefix prefix for the name of output file, excluding the directory,
#' defaults to 'drugage'
#' @param outdir directory to save the file, defaults to the current directory
#'
#' @return returns the file path
#' @export
download_drugage <- function(outprefix = "drugage", outdir = "./") {
  local_locale(c("LC_TIME" = "C"))
  local_timezone("UTC")
  local_dir(outdir)
  time <- format(Sys.time(), "%Y%m%d%H%M%S")
  system("wget https://genomics.senescence.info/drugs/dataset.zip")
  system("unzip dataset.zip")
  system("rm dataset.zip")
  outfile <- paste(outprefix, "_", time, ".csv", sep = "")
  fname <- grep("drug", grep(".csv", dir_ls(), value = T), ignore.case = T,
                value = T)
  system(paste("mv", fname, outfile, sep = " "))
  system(paste("mv release.html ", outprefix, "_release_", time, ".html",
    sep = ""
  ))
  outfile <- path(outdir, outfile)
  return(outfile)
}



#' Tidy DrugAge Data
#'
#' @param path2file path to drugage csv file
#'
#' @return returns a data frame with 10 columns: 'organism','dataset','compound_name',
#' 'strain','dosage','avg_lifespan_change','max_lifespan_change','gender','significance','pubmed_id'
#' @export
#'
tidy_drugage <- function(path2file){
  gendermap <- setNames(c('Male','Female','Male','Female','Both','Hermaphrodite','Pooled','Mixed',NA),
                        c('Male','Female','MALE','FEMALE','BOTH','Hermaphrodite','Pooled','Mixed','Not Sp'))
  signifmap <- setNames(c('ns','s','ns'),
                        c('NS','S','0'))
  dat <- read_csv(path2file)
  dat <- dat %>% select(species) %>% unique() %>%
    separate(species, into = c('genus','specific'), remove = F) %>%
    mutate(dataset = ifelse(!is.na(specific),
                            tolower(paste(substr(genus,1,1),specific,sep='')),
                            NA)) %>%
    select(species, dataset) %>%
    right_join(dat) %>%
    mutate(gender = gendermap[as.character(gender)],
           significance = signifmap[as.character(significance)]) %>%
    rename(organism = species) %>%
    rename(chemName = compound_name)
  return(dat)
}

#' Download CellAge Data
#'
#' @param outprefix prefix for the name of output file, excluding the directory,
#' defaults to 'cellage'
#' @param outdir directory to save the file, defaults to the current directory
#'
#' @return returns the file path
#' @export
download_cellage <- function(outprefix = "cellage", outdir = "./") {
  local_locale(c("LC_TIME" = "C"))
  local_timezone("UTC")
  local_dir(outdir)
  time <- format(Sys.time(), "%Y%m%d%H%M%S")
  system("wget https://genomics.senescence.info/cells/cellAge.zip")
  system("unzip cellAge.zip")
  system("rm cellAge.zip")
  outfile <- paste(outprefix, "_", time, ".csv", sep = "")
  fname <- grep("cell", grep(".csv", dir_ls(), value = T), ignore.case = T,
                value = T)
  system(paste("mv", fname, outfile, sep = " "))
  system(paste("mv release.html ", outprefix, "_release_", time, ".html",
               sep = ""))
  outfile <- path(outdir, outfile)
  return(outfile)
}

#' Tidy CellAge Data
#'
#' @param path2file path to cellage csv file
#'
#' @return returns a data frame with 7 columns: organism, dataset, Gene,
#' cancer_type, senescence_effect, description, notes
#' @export
#'
tidy_cellage <- function(path2file){
  organismmap = setNames(c('Homo sapiens'),c('Human'))
  cancermap = setNames(c('-','Cancer','Cancer'),
                       c('No','yes','Yes'))
  dat <- read_delim(path2file, delim = ';') %>%
    select(2,5,6,7,8,9) %>%
    rename(Gene = gene_name) %>%
    mutate(organism = organismmap[organism])
  dat <- dat %>% select(organism) %>% unique() %>%
    separate(organism, into = c('genus','specific'), remove = F) %>%
    mutate(dataset = ifelse(!is.na(specific),
                            tolower(paste(substr(genus,1,1),specific,sep='')),
                            NA)) %>%
    select(organism, dataset) %>%
    right_join(dat) %>%
    mutate(cancer_type = cancermap[cancer_type])
  return(dat)
}

#' Download GenAge, DrugAge, and CellAge Data
#'
#' @param outdir directory to save the file, defaults to the current directory
#'
#' @return returns a data.frame with four columns including file paths to
#' GenAge Human (`human`), GenAge Model Organism (`model`), DrugAge (`drugage`),
#' and CellAge (`cellage`) files.
#'
#' @export
download_HAGR <- function(outdir = "./") {
  human <- download_genage_human(outdir = outdir)
  model <- download_genage_model(outdir = outdir)
  drugage <- download_drugage(outdir = outdir)
  cellage <- download_cellage(outdir = outdir)
  HAGR_fn <- data.frame(human, model, drugage, cellage)
  return(HAGR_fn)
}
