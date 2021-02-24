#' Download GenAge Human Data
#'
#' @param outprefix prefix for the name of output file, excluding the directory,
#' defaults to 'genage_human'
#' @param outdir directory to save the file, defaults to the current directory
#'
#' @return returns the file path
#' @export
download_genage_human <- function(outprefix = "genage_human", outdir = "./") {
  withr::local_locale(c("LC_TIME" = "C"))
  withr::local_timezone("UTC")
  withr::local_dir(outdir)
  time <- format(Sys.time(), "%Y%m%d%H%M%S")
  system("wget https://genomics.senescence.info/genes/human_genes.zip")
  system("unzip human_genes")
  system("rm human_genes.zip")
  outfile <- paste(outprefix, "_", time, ".csv", sep = "")
  fname <- grep("human", grep(".csv", fs::dir_ls(), value = T), ignore.case = T,
                value = T)
  system(paste("mv", fname, outfile, sep = " "))
  system(paste("mv release.html ", outprefix, "_release_", time, ".html",
    sep = ""
  ))
  outfile <- fs::path(outdir, outfile)
  return(outfile)
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

#' Download DrugAge Data
#'
#' @param outprefix prefix for the name of output file, excluding the directory,
#' defaults to 'drugage'
#' @param outdir directory to save the file, defaults to the current directory
#'
#' @return returns the file path
#' @export
download_drugage <- function(outprefix = "drugage", outdir = "./") {
  withr::local_locale(c("LC_TIME" = "C"))
  withr::local_timezone("UTC")
  withr::local_dir(outdir)
  time <- format(Sys.time(), "%Y%m%d%H%M%S")
  system("wget https://genomics.senescence.info/drugs/dataset.zip")
  system("unzip dataset.zip")
  system("rm dataset.zip")
  outfile <- paste(outprefix, "_", time, ".csv", sep = "")
  fname <- grep("drug", grep(".csv", fs::dir_ls(), value = T), ignore.case = T,
                value = T)
  system(paste("mv", fname, outfile, sep = " "))
  system(paste("mv release.html ", outprefix, "_release_", time, ".html",
    sep = ""
  ))
  outfile <- fs::path(outdir, outfile)
  return(outfile)
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
  withr::local_locale(c("LC_TIME" = "C"))
  withr::local_timezone("UTC")
  withr::local_dir(outdir)
  time <- format(Sys.time(), "%Y%m%d%H%M%S")
  system("wget https://genomics.senescence.info/cells/cellAge.zip")
  system("unzip cellAge.zip")
  system("rm cellAge.zip")
  outfile <- paste(outprefix, "_", time, ".csv", sep = "")
  fname <- grep("cell", grep(".csv", fs::dir_ls(), value = T), ignore.case = T,
                value = T)
  system(paste("mv", fname, outfile, sep = " "))
  system(paste("mv release.html ", outprefix, "_release_", time, ".html",
               sep = ""))
  outfile <- fs::path(outdir, outfile)
  return(outfile)
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
