
#Function to remove read count and index files created from main script.

#' File delete function
#'
#' @description Function to delete files generated in the process.
#' @param pattern The filename extension to be deleted.
#' @keywords pattern
#' @export del_files
#' @examples
#' del_files()

del_files<- function(pattern){
  fl_rm<-list.files(".", pattern = pattern, all.files = F, full.names = F)
  for(i in fl_rm){
    i<-paste("rm", i , sep = " ")
    system(i)
  }
}


