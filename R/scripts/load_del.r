

#Function to download, install and load the required libraries only when needed.

#' Package download function
#'
#' Function to download, install and load the required libraries only when needed
#' @param x The name of the package you wish to download and install? Defaults NULL.
#' @keywords cats
#' @export
#' @examples
#' packages()

packages <- function(x) {
  x <- as.character(match.call()[[2]])
  if (!require(x, character.only = TRUE)) {
    install.packages(pkgs = x, repos = "http://cran.r-project.org")
    require(x, character.only = TRUE)
  }
}



#Function to remove read count and index files created from main script.

#' File delete function
#'
#' Function to delete files generated in the process.
#' @param pattern The filename extension to be deleted.
#' @keywords pattern
#' @export
#' @examples
#' del_files()

del_files<- function(pattern){
  fl_rm<-list.files(".", pattern = pattern, all.files = F, full.names = F)
  for(i in fl_rm){
    i<-paste("rm", i , sep = " ")
    system(i)
  }
}


