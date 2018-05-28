
#' RACEseqR main function
#'
#' @description Function to generate the alignment from a RACE sequencing experiment
#' @param input_data The input .fastq files to be used for the alignment.
#' @param replicon_ref The reference sequence to be used in the alignment.
#' @param filename The filename that will be used for the output pdf file.
#' @param mismatch The number of mismatches to be accepted during alignment. DEFAULT = 0
#' @param RACE_adapter The RACE adapter to be trimmied using cutadapt. 
#' @param str The start position of the binding region.
#' @param end The end position of the binding region.
#' @param tmap The option to use the Tmap aligner instead of the default BOWTIE aligner. DEFAULT = FALSE
#' @param plot The option to output the results in plot. The default output is an R dataframe object. DEFAULT = FALSE
#' @keywords binding_region
#' @export RACEseq
#' @examples
#' RACEseq(filename = "foo", str = 9478, end = 9498, plot = T)
#' region_foo <- RACEseq(input_data = foo.fastq, replicon_ref = "foo.fasta", filename = "bar", RACE_adapter = "GGACACTGACATGGACTGAAGGAGTAGAAA", str = 9478, end = 9498)



RACEseq<- function(input_data, 
                   replicon_ref,
                   filename,
                   mismatch = 0, 
                   RACE_adapter, 
                   str=NULL, 
                   end=NULL, 
                   tmap =NULL, 
                   plot=NULL) {
  
  if(missing(filename)) stop("Please specify an output filename")
  if(missing(str)) stop("str value must be set")
  if(missing(end)) stop("end value must be set")
  if(!str %in% seq(1, 10000000, by = 1)) stop("str value must be >= 1)")
  if(!end %in% seq(1, 10000000, by = 1)) stop("end value must be >= 1)")

  if(missing(input_data)) {
    #reading input data from  working dir
    input_data<- as.character(list.files(".", pattern ="fastq", all.files = F, full.names = F))
    if ((length(input_data))==0) {
      stop("No input .fastq data file available")
    } else if ((length(input_data))>=2) {
      stop("More than one input .fastq data file available")
    }
  }
  if(missing(replicon_ref)) {
    #reading ref sequence from  working dir
    replicon_ref<- as.character(list.files(".", pattern ="fasta", all.files = F, full.names = F))
    if ((length(replicon_ref))==0) {
      stop("No input .fasta reference file available")
    } else if ((length(replicon_ref))>=2) {
      stop("More than one .fasta reference file")
    }
  }

  if(missing(tmap)) {
    bowtie_align(input_data, replicon_ref, mismatch, RACE_adapter, out_name = "alignment")
  } else if (tmap==F) {
    bowtie_align(input_data, replicon_ref, mismatch, RACE_adapter, out_name = "alignment")
    } else 
    tmap_align(input_data, replicon_ref, mismatch, RACE_adapter, out_name = "alignment")


  binding_region <- datafile_out(str, end, replicon_ref, filename)

  del_files("read_counts")
  del_files("fasta.tmap.")
  del_files("aligned.bam")
  del_files("out.sam")
  del_files("index")

  if(!missing(plot)) {
    if(plot==T) {
      plot_out(binding_region, filename)
    } 
  }
  return(binding_region)
}






