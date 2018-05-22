
#' RACEseq TMAP alignment function
#'
#' Function to perform the trimming of the RACE adapter using CUTADAPT and the alignment of the sequencing data from the .fastq files using the TMAP aligner.
#' @param input_data The input .fastq files to be used for the alignment.
#' @param replicon_ref The reference sequence to be used in the alignment.
#' @param mismatch The number of mismatches to be allowed during alignment.
#' @param RACE_adapter The RACE adapter nucleotide sequence to be trimmed from the reads using cutadapt. 
#' @param out_name The output filename from the alignment.
#' @keywords alignment tmap cutadapt
#' @export
#' @examples
#' tmap_align()
#' @seealso bowtie_align

tmap_align <- function(input_data, replicon_ref,mismatch,  RACE_adapter, out_name){
  
  if(missing(mismatch)) mismatch<- 0
  if(missing(input_data)) stop("No input .fastq data file available")
  if(missing(replicon_ref)) stop("No input .fasta reference file available")
  
  if(missing(RACE_adapter)) {
    #build the index
    CMD_tmapindex<- paste("tmap index -f", replicon_ref , sep=" ")
    system(CMD_tmapindex)
    
    #perform alignment with tmap and read count using bedtools with no adapter trimming
    print(paste0("Performing alignment with ", mismatch, " mismatch using tmap"))
    CMD_tmap<- paste("tmap map1 -a 0 -g 3 --max-mismatches ",mismatch," -f ", replicon_ref," -r ", input_data, " | samtools view -bt ", replicon_ref," - | genomeCoverageBed -d -5 -ibam stdin > read_counts_mm",mismatch,"_",out_name,".txt", sep="")
    system(CMD_tmap)
    
  } else {
    #build the index
    CMD_tmapindex<- paste("tmap index -f", replicon_ref , sep=" ")
    system(CMD_tmapindex)
    
    #adapter trimming using cutadapt and alignment with tmap and read count using bedtools
    print(paste0("Performing adapter trimming and alignment with ", mismatch, " mismatch using tmap"))
    CMD_tmap<- paste("cutadapt -g ", RACE_adapter, " -e0 --no-indels -m10 --discard-untrimmed --quiet ", input_data," |tmap map1 -a 0 -g 3 --max-mismatches ",mismatch," -f ", replicon_ref," -i fastq | samtools view -bt ", replicon_ref," - | genomeCoverageBed -d -5 -ibam stdin > read_counts_mm",mismatch,"_",out_name,".txt", sep="")
    system(CMD_tmap) 
  }
}




