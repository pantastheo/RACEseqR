
#' RACEseqR TMAP alignment function
#'
#' Function to perform the alignment of the sequencing data from the .fastq files using the TMAP aligner.
#' @param input_data The input .fastq files to be used for the alignment.
#' @param replicon_ref The reference sequence to be used in the alignment.
#' @param mismatch The number of mismatches to be allowed during alignment.
#' @param out_name The output filename from the alignment.
#' @keywords alignment tmap
#' @export
#' @examples
#' tmap_align()
#' @seealso bowtie_align
#' @seealso bowtie_trim_align
#' @seealso tmap_trim_align

tmap_align<- function(input_data, replicon_ref, mismatch, out_name){
  
  if(missing(mismatch)) mismatch<- 0
  if(missing(out_name)) out_name<- "alignment"
  if(missing(input_data)) {
    #reading ref sequence from  working dir
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
  
  #build the index
  CMD_tmapindex<- paste("tmap index -f", replicon_ref , sep=" ")
  system(CMD_tmapindex)

  #perform alignment with tmap and read count using bedtools with no adapter trimming
  print(paste0("Performing alignment with ", mismatch, " mismatch using tmap"))
  CMD_tmap<- paste("tmap map1 -a 0 -g 3 --max-mismatches ",mismatch," -f ", replicon_ref," -r ", input_data, " | samtools view -bt ", replicon_ref," - | genomeCoverageBed -d -5 -ibam stdin > read_counts_mm",mismatch,"_",out_name,".txt", sep="")
  system(CMD_tmap)
}




#' RACEseqR TMAP alignment function
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
#' tmap_trim_align()
#' @seealso bowtie_align
#' @seealso bowtie_trim_align
#' @seealso tmap_align

tmap_trim_align <- function(input_data, replicon_ref,mismatch,  RACE_adapter, out_name){
  
  if(missing(RACE_adapter)) stop("No RACE adapter sequence available for trimming")
  if(missing(mismatch)) mismatch<- 0
  if(missing(out_name)) out_name<- "alignment"
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
  
  #build the index
  CMD_tmapindex<- paste("tmap index -f", replicon_ref , sep=" ")
  system(CMD_tmapindex)
  
  #adapter trimming using cutadapt and alignment with tmap and read count using bedtools
  print(paste0("Performing adapter trimming and alignment with ", mismatch, " mismatch using tmap"))
  CMD_tmap<- paste("cutadapt -g ", RACE_adapter, " -e0 --no-indels -m10 --discard-untrimmed --quiet ", input_data," |tmap map1 -a 0 -g 3 --max-mismatches ",mismatch," -f ", replicon_ref," -i fastq | samtools view -bt ", replicon_ref," - | genomeCoverageBed -d -5 -ibam stdin > read_counts_mm",mismatch,"_",out_name,".txt", sep="")
  system(CMD_tmap) 
}



#' RACEseqR BOWTIE alignment function
#'
#' Function to perform the alignment of the sequencing data from the .fastq files using the BOWTIE aligner.
#' @param input_data The input .fastq files to be used for the alignment.
#' @param replicon_ref The reference sequence to be used in the alignment.
#' @param mismatch The number of mismatches to be allowed during alignment.
#' @param out_name The output filename from the alignment.
#' @keywords alignment bowtie 
#' @export
#' @examples
#' bowtie_align()
#' @seealso bowtie_trim_align
#' @seealso tmap_align
#' @seealso tmap_trim_align

bowtie_align<- function(input_data, replicon_ref, mismatch, out_name){

  if(missing(mismatch)) mismatch<- 0
  if(missing(out_name)) out_name<- "alignment"
  if(missing(input_data)) {
    #reading ref sequence from  working dir
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
  
  #build the index 
  CMD_bowindex<- paste("bowtie-build -q -f", replicon_ref, "index", sep=" ")
  system(CMD_bowindex)

  #perform alignment with bowtie and read count using bedtools with no adapter trimming
  print(paste0("Performing alignment with ", mismatch, " mismatch using bowtie"))
  CMD_bow<- paste("bowtie -p 2 -S -k 1 -v ", mismatch, " index ", input_data," | samtools view -bS - | genomeCoverageBed -d -5 -ibam stdin > read_counts_mm",mismatch,"_",out_name,".txt", sep="")
  system(CMD_bow)
}



#' RACEseqR BOWTIE alignment function
#'
#' Function to perform the trimming of the RACE adapter using CUTADAPT and the alignment of the sequencing data from the .fastq files using the BOWTIE aligner.
#' @param input_data The input .fastq files to be used for the alignment.
#' @param replicon_ref The reference sequence to be used in the alignment.
#' @param mismatch The number of mismatches to be allowed during alignment.
#' @param RACE_adapter The RACE adapter nucleotide sequence to be trimmed from the reads using cutadapt. 
#' @param out_name The output filename from the alignment.
#' @keywords alignment bowtie cutadapt
#' @export
#' @examples
#' bowtie_trim_align()
#' @seealso bowtie_align
#' @seealso tmap_align
#' @seealso tmap_trim_align

bowtie_trim_align <- function(input_data, replicon_ref, mismatch, RACE_adapter, out_name){

  if(missing(RACE_adapter)) stop("No RACE adapter sequence available for trimming")
  if(missing(mismatch)) mismatch<- 0
  if(missing(out_name)) out_name<- "alignment"
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
  
  #build the index 
  CMD_bowindex<- paste("bowtie-build -q -f", replicon_ref, "index", sep=" ")
  system(CMD_bowindex)
  
  #adapter trimming using cutadapt and alignment with bowtie and read count using bedtools
  print(paste0("Performing adapter trimming and alignment with ", mismatch, " mismatch using bowtie"))
  CMD_bow<- paste("cutadapt -g ", RACE_adapter, " -e0 --no-indels -m10 --discard-untrimmed --quiet ", input_data," |bowtie -p 8 -S -k 1 -v ", mismatch, " index - | samtools view -bS - | genomeCoverageBed -d -5 -ibam stdin > read_counts_mm",mismatch,"_",out_name,".txt", sep="")
  system(CMD_bow)
}

