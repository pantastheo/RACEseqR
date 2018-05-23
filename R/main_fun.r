

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

main_RACE_seq_R<- function(input_data, replicon_ref, mismatch, RACE_adapter, str, end, filename, aligner ) {

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

  if(missing(mismatch)) mismatch<- 0
  out_name<- "alignment"

  if(missing(aligner)) {
    bowtie_align(input_data, replicon_ref, mismatch, RACE_adapter, out_name)
  } else {
    tmap_align(input_data, replicon_ref, mismatch, RACE_adapter, out_name)
  }

  binding_region <- datafile_out(str, end, replicon_ref)

  del_files("read_count_")
  del_files("fasta.tmap.")
  del_files("aligned.bam")
  del_files("out.sam")
  del_files("index")

  plot_out(binding_region, filename)

  return(binding_region)
}

