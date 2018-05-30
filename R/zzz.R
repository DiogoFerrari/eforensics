# required (by devtools) to link the cpp code 
#' @importFrom magrittr %>%
#' @importFrom foreach %dopar%

.onAttach <- function(libname, pkgname){
    packageStartupMessage('

 ---------------------------------------
 Election Forensics Package (eforensics)
 ---------------------------------------

 Authors:

 Supported by NSF grant SES 1523355


 ')
}

if(getRversion() >= "2.15.1")  utils::globalVariables(c(".", "contrasts",
                                                        "lower", "upper",
                                                        "Mean", "SD", "Parameter"))
