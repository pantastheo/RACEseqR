
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


datafile_out<- function(str, end, replicon_ref) {
  
  #setting function call conditions
  if(missing(str)) stop("str value must be set")
  if(missing(end)) stop("end value must be set")
  if(!str %in% seq(1, 10000000, by = 1)) stop("str value must be >= 1)")
  if(!end %in% seq(1, 10000000, by = 1)) stop("end value must be >= 1)")
  
  #reading ref sequence name from function call
  if(missing(replicon_ref)) {
    #reading ref sequence from  working dir
    replicon_ref<- as.character(list.files(".", pattern ="fasta", all.files = F, full.names = F))
    if ((length(replicon_ref))==0) {
      stop("No input .fasta reference file available")
    } else if ((length(replicon_ref))>=2) {
      stop("More than one .fasta reference file")
    }
  }
  
  #transforming reference sequence
  nt_reference <-strsplit((toString(readBStringSet(replicon_ref))), NULL , fixed = T)
  nt_reference<- data.frame(lapply(nt_reference, function(x) toupper(x)), stringsAsFactors = F)
  
  #read the output  file
  #input reads in .txt format
  out_reads<- list.files(".", pattern =".txt", all.files = F, full.names = F)
  if ((length(out_reads))==0) {
    stop("No output reads file available")
  } else if ((length(out_reads))>=2) {
    stop("More than one output reads file")
  }
  
  #create dataframe with reference and reads
  
  dataframe<- data.frame(reads, nt_reference , stringsAsFactors = F)
  
  #calculating the % and log10 columns
  dataframe[,5] <- (dataframe[,3]/sum(dataframe[,3])*100)
  dataframe[,6] <- (log10(dataframe[,3]))
  dataframe[dataframe== -Inf] <-0
  
  #focusing on target region can be ajusted acording to experiment
  binding_region <- dataframe[str:end,]
  
  out_reads<- "datafile_mm0.txt"
  
  out_name<- file_path_sans_ext(((strsplit(out_reads, "_")) [[1]])[[2]])
  
  write.table(binding_region, file = paste0("datafile_",out_name, ".txt") , sep = "\t", col.names = c("reference", "position", "count", "nucleotide", "percentage", "log10" ), row.names = F )
  }
  