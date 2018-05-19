
#' RACEseqR dataframe function
#'
#' Function to generate and output the RACE seq dataframe of the specified start and end positions.
#' @param str The start position of the binding region.
#' @param end The end position of the binding region.
#' @param df_data The read count data generated from the alignment 
#' @param replicon_ref The reference sequence used in the alignment.
#' @keywords dataframe output
#' @export
#' @examples
#' out_csv()


out_csv<- function(str, end, replicon_ref) {
    
  #read the output  file
  
  #input reads in .txt format
  out_reads<- list.files(".", pattern ="read_count_", all.files = F, full.names = F)
  if ((length(out_reads))==0) {
    stop("No output reads file available")
  } 
  else if ((length(out_reads))>=2) {
    stop("More than one output reads file")
  } 
  else if ((length(out_reads))==1) {
    out_reads<- as.character(out_reads)
    reads<- read.delim(out_reads, header = F )
  }
  else { 
    stop("No read counts file available")}
  
  
  #reading and transforming reference sequence
  if(is.na(replicon_ref))
    #input the reference sequence in .fasta format
    refname<- list.files(".", pattern ="fasta", all.files = F, full.names = F)
  if ((length(refname))==0) {
    stop("No input .fasta reference file available")
  } 
  else if ((length(refname))>=2) {
    stop("More than one reference file")
  } 
  else if ((length(refname))==1) {
    replicon_ref<- as.character(refname)
    nt_reference <-strsplit((toString(readBStringSet(replicon_ref))), NULL , fixed = T)
    nt_reference<- data.frame(lapply(nt_reference, function(x) toupper(x)), stringsAsFactors = F)
  }
  else {
    nt_reference <-strsplit((toString(readBStringSet(replicon_ref))), NULL , fixed = T)
    nt_reference<- data.frame(lapply(nt_reference, function(x) toupper(x)), stringsAsFactors = F)}
  
  
  
  #create dataframe with reference and reads
  
  dataframe<- data.frame(reads, nt_reference , stringsAsFactors = F)
  
  #calculating the % and log10 columns
  dataframe[,5] <- (dataframe[,3]/sum(dataframe[,3])*100)
  dataframe[,6] <- (log10(dataframe[,3]))
  dataframe[dataframe== -Inf] <-0
  
  #focusing on target region can be ajusted acording to experiment
  binding_region <- dataframe[str:end,]
  
  outfilename<- file_path_sans_ext(((strsplit(out_reads, "_")) [[1]])[[3]])
  
  write.table(binding_region, file = paste0(outfilename, ".txt") , sep = "\t", col.names = c("reference", "position", "count", "nucleotide", "percentage", "log10" ), row.names = F )
}

