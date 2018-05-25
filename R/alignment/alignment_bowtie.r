
#' RACEseq BOWTIE alignment function
#'
#' Function to perform the trimming of the RACE adapter using CUTADAPT and the alignment of the sequencing data from the .fastq files using the BOWTIE aligner.
#' @param input_data The input .fastq files to be used for the alignment.
#' @param replicon_ref The reference sequence to be used in the alignment.
#' @param mismatch The number of mismatches to be allowed during alignment.
#' @param RACE_adapter The RACE adapter nucleotide sequence to be trimmed from the reads using cutadapt.
#' @param out_name The output filename from the alignment.
#' @keywords alignment bowtie cutadapt
#' @export bowtie_align
#' @examples
#' bowtie_align()
#' @seealso tmap_align

bowtie_align <- function(input_data, replicon_ref, mismatch, RACE_adapter, out_name){


  if(missing(mismatch)) mismatch<- 0
  if(missing(input_data)) stop("No input .fastq data file available")
  if(missing(replicon_ref)) stop("No input .fasta reference file available")

  if(missing(RACE_adapter)) {
    #build the index
    CMD_bowindex<- paste("bowtie-build -q -f", replicon_ref, "index", sep=" ")
    system(CMD_bowindex)

    #perform alignment with bowtie and read count using bedtools with no adapter trimming
    print(paste0("Performing alignment with ", mismatch, " mismatch using bowtie"))
    CMD_bow<- paste("bowtie -p 2 -S -k 1 -v ", mismatch, " index ", input_data," | samtools view -bS - | genomeCoverageBed -d -5 -ibam stdin > read_counts_mm",mismatch,"_",out_name,".txt", sep="")
    system(CMD_bow)

  } else {
    #build the index
    CMD_bowindex<- paste("bowtie-build -q -f", replicon_ref, "index", sep=" ")
    system(CMD_bowindex)

    #adapter trimming using cutadapt and alignment with bowtie and read count using bedtools
    print(paste0("Performing adapter trimming and alignment with ", mismatch, " mismatch using bowtie"))
    CMD_bow<- paste("cutadapt -g ", RACE_adapter, " -e0 --no-indels -m10 --discard-untrimmed --quiet ", input_data," |bowtie -p 8 -S -k 1 -v ", mismatch, " index - | samtools view -bS - | genomeCoverageBed -d -5 -ibam stdin > read_counts_mm",mismatch,"_",out_name,".txt", sep="")
    system(CMD_bow)
  }

}

