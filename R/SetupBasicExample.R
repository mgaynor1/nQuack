#' @title Setup Basic Example
#' @description This function was made to download all files needed to run the Basic Example vignette.
#'  It downloads a zipped file from Zenodo and unzips it to the "inst/extdata/" directory  within the nQuack package.
#'
#' @param overwrite Logical. If TRUE, the function will overwrite any existing files
#' in the data directory. Default is FALSE.
#' @examples
#' if(exists("crazy")){
#'    SetupBasicExample(overwrite = FALSE)
#' }
#' @returns This function downloads files to the "inst/extdata/" directory within the nQuack package
#' and does not return any value. The function will
#' print messages as it downloads to keep you updated on the progress.
#'
#' @importFrom tools R_user_dir
#' @importFrom httr2 request req_perform resp_url req_timeout
#' @importFrom utils unzip
#' @import magrittr
#' @export

SetupBasicExample <- function(overwrite = FALSE){

  ## Setup directories
  dbdir <- tools::R_user_dir("nQuack")
  aimdir <- paste0(dbdir, "/inst/extdata")
  writedir <-   paste0(dbdir, "/inst/check")

  ## Check if downloaded prior
  if(dir.exists(paste0(aimdir, "/01_raw/"))){
      if(overwrite == FALSE){
        message("The basic example has already been setup. You are ready to rumble!")
        token <- 0
      }else{
        message("The basic example has already been setup, but you have chosen to overwrite it.")
        token <- 1
      }
  }else{
    message("Setting up the basic example for the first time. This may take a few minutes.")
    token <- 1
  }
  ## Download and unzip
  if(token == 1){
    message("Checking Zenodo DOI.")
    ### Zenodo DOI
    doi <- "10.5281/zenodo.21074264"

    landing_page <- tryCatch(
      httr2::request(paste0("https://doi.org/", doi)) |>
        httr2::req_perform() |>
        httr2::resp_url(),
      error = function(e) {
        stop("Could not resolve DOI ", doi, ".", call. = FALSE)
      }
    )

    # Extract record ID
    record_id <- sub(".*/records/([0-9]+).*", "\\1", landing_page)

    # Query Zenodo API
    url <- sprintf(
      "https://zenodo.org/records/%s/files/extdata.zip?download=1",
      record_id)

    message("Starting download of extdata.zip from Zenodo. This may take a few minutes.")
    zipfile <- tempfile(fileext = ".zip")

    tryCatch(
      httr2::request(url) |>
        httr2::req_timeout(seconds = 3600) |>
        httr2::req_perform(path = zipfile),
      error = function(e) {
        stop("Download from Zenodo failed:\n", conditionMessage(e),
             call. = FALSE)
      }
    )

    if (!file.exists(zipfile)) {
      stop("Download failed: zip file was not created.",
           call. = FALSE)
    }
    message("Zipfile Downloaded! Time to unzip!")
    utils::unzip(zipfile,
          exdir = writedir,
          overwrite = TRUE)

    if(!dir.exists(paste0(aimdir, "/01_raw/"))){
      stop("File did not unzip properly.",  call. = FALSE)
    }

      # Clean up
      unlink(zipfile)
      message("Done!")

    }else{
      message("Done!")
    }

}
