
#' RACEseqR dataframe function
#'
#' Function to generate and output the RACE seq dataframe of the specified start and end positions.
#' @param str The start position of the binding region.
#' @param end The end position of the binding region.
#' @param replicon_ref The reference sequence used in the alignment.
#' @param filename The filename that will be used for the output data file.
#' @keywords dataframe output
#' @export datafile_out
#' @examples
#' out_csv()


datafile_out<- function(str, end, replicon_ref, filename) {
  
  #setting function call conditions
  if(missing(str)) stop("str value must be set")
  if(missing(end)) stop("end value must be set")
  if(!str %in% seq(1, 10000000, by = 1)) stop("str value must be >= 1)")
  if(!end %in% seq(1, 10000000, by = 1)) stop("end value must be >= 1)")
  
  #reading ref sequence name from function call
  if(missing(replicon_ref)) stop("No input .fasta reference file available")
  
  #transforming reference sequence
  nt_reference <-strsplit((toString(readBStringSet(replicon_ref))), NULL , fixed = T)
  nt_reference<- data.frame(lapply(nt_reference, function(x) toupper(x)), stringsAsFactors = F)
  
  #read the output  file
  #input reads in .txt format
  out_reads<- list.files(".", pattern ="alignment", all.files = F, full.names = F)
  if ((length(out_reads))==0) {
    stop("No output alignment file available")
  } else if ((length(out_reads))>=2) {
    stop("More than one output alignment file available")
  } else reads<- read.delim(out_reads, header = F, quote = " ")
  
  #create dataframe with reference and reads
  
  dataframe<- data.frame(reads, nt_reference , stringsAsFactors = F)
  
  #calculating the % and log10 columns
  dataframe[,5] <- (dataframe[,3]/sum(dataframe[,3])*100)
  dataframe[,6] <- (log10(dataframe[,3]))
  dataframe[dataframe== -Inf] <-0
  
  #focusing on target region can be ajusted acording to experiment
  binding_region <- dataframe[str:end,]
  colnames(binding_region)<- c("reference", "position", "count", "nucleotide", "percentage", "log10" )
  
  out_name<- file_path_sans_ext(((strsplit(out_reads, "_")) [[1]])[[2]])
  
  write.table(binding_region, file = paste0("datafile_",filename ,"_",out_name, ".txt") , sep = "\t", col.names = c("reference", "position", "count", "nucleotide", "percentage", "log10" ), row.names = F , quote = FALSE)
  
  return(binding_region)
}
  