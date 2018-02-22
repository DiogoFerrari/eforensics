# required (by devtools) to link the cpp code 
#' @useDynLib eforensics
#' @importFrom Rcpp sourceCpp
NULL 

#' @importFrom magrittr %>%


.onAttach <- function(libname, pkgname){
    packageStartupMessage('

 ---------------------------------------
 Election Forensics Package (eforensics)
 ---------------------------------------

 Authors:

 Supported by NSF grant SES 1523355


 ')
}
