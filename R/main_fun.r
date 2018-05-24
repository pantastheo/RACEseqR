

#Main function to call all sub founctions in main source directory R.
source_dir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}

source_dir("R/alignment")
source_dir("R/scripts")

#List of required libraries to be loaded
suppressMessages(packages(Biostrings))
suppressMessages(packages(tools))

#' RACEseq main function
#'
#' Function to generate the alignment from a RACE sequencing experiment
#' @param input_data The input .fastq files to be used for the alignment.
#' @param replicon_ref The reference sequence to be used in the alignment.
#' @param str The start position of the binding region.
#' @param end The end position of the binding region.
#' @param mismatch The number of mismatches to be accepted during alignment.
#' @param RACE_adapter The RACE adapter to be trimmied using cutadapt .
#' @param tmap The option to use the Tmap aligner. Defaults NULL
#' @param filename The filename that will be used for the output pdf file.
#' @keywords binding_region
#' @export
#' @examples
#' RACEseq()



RACEseq<- function(input_data, replicon_ref, mismatch = 0, RACE_adapter=NULL, str, end, filename, tmap ) {

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
  } else {
    tmap_align(input_data, replicon_ref, mismatch, RACE_adapter, out_name = "alignment")
  }

  binding_region <- datafile_out(str, end, replicon_ref)

  del_files("read_counts")
  del_files("fasta.tmap.")
  del_files("aligned.bam")
  del_files("out.sam")
  del_files("index")

  plot_out(binding_region, filename)

  return(binding_region)
}

