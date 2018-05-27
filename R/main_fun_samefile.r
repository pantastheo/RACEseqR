
# all the scripts together just testing the functions



#Function to download, install and load the required libraries only when needed.

#' Package download function
#'
#' Function to download, install and load the required libraries only when needed
#' @param x The name of the package you wish to download and install? Defaults NULL.
#' @keywords cats
#' @export
#' @examples
#' packages()

packages <- function(x) {
  x <- as.character(match.call()[[2]])
  if (!require(x, character.only = TRUE)) {
    install.packages(pkgs = x, repos = "http://cran.r-project.org")
    require(x, character.only = TRUE)
  }
}



#Function to remove read count and index files created from main script.

#' File delete function
#'
#' Function to delete files generated in the process.
#' @param pattern The filename extension to be deleted.
#' @keywords pattern
#' @export
#' @examples
#' del_files()

del_files<- function(pattern){
  fl_rm<-list.files(".", pattern = pattern, all.files = F, full.names = F)
  for(i in fl_rm){
    i<-paste("rm", i , sep = " ")
    system(i)
  }
}


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
  if(missing(replicon_ref)) stop("No input .fasta reference file available")

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

  out_name<- file_path_sans_ext(((strsplit(out_reads, "_")) [[1]])[[2]])

  write.table(binding_region, file = paste0("datafile_",out_name, ".txt") , sep = "\t", col.names = c("reference", "position", "count", "nucleotide", "percentage", "log10" ), row.names = F )

  return(binding_region)
}



#' RACEseqR plot function
#'
#' @description Function to generate and output the RACE seq plot of the specified start and end positions.
#' @param binding_region The region that will be plotted in the graph. Can be an R dataframe or a tab delim file in your working directory.
#' @param filename The filename that will be used for the output pdf file.
#' @keywords binding_region
#' @export
#' @examples
#' out_plot()

plot_out<- function(binding_region, filename) {

  #setting function call conditions
  if(missing(filename)) filename<- "RACE_graph"
  if(missing(binding_region)){
    #read the output binding region tab delim file
    #input reads from working dir
    datafile_input<- list.files(".", pattern ="datafile_", all.files = F, full.names = F)
    if ((length(datafile_input))==0) {
      stop("No output reads .txt file available")
    } else if ((length(datafile_input))>=2) {
      stop("More than one output reads .txt file")
    }
    binding_region<- read.csv(datafile_input, header = F )
    out_name<- file_path_sans_ext(((strsplit(datafile_input, "_")) [[1]])[[2]])
    filename<- paste0(filename,"_", out_name)
  }



  #create wildtype linear & log scale graph

  pdf(paste0(filename,".pdf"), width=15)

  mp <- barplot(binding_region[,5],
                ylab="Novel 5\' Ends (%)",
                xlab="Binding site",
                names.arg=(binding_region[,4]),
                las=1,
                cex.names = 2.2,
                col="darkgrey" ,
                main="Novel 5' Ends in linear",
                cex.main=2.3,
                cex.lab=1.3,
                ylim=c(0,100))
  text(mp,binding_region[,5]+5 ,cex = 1.3, adj = 0 ,labels=binding_region[,3] ,srt=90)

  #in log10  logarithmic scale
  mp <- barplot(binding_region[,6],
                ylab=expression("Novel 5\' Ends (log"[10]*")"),
                xlab="Binding site",
                names.arg=(binding_region[,4]),
                las=1,
                cex.names = 2.2,
                col="darkgrey",
                main="Novel 5' Ends in logarithmic",
                cex.main=2.3,
                cex.lab=1.3,
                ylim=c(0,10))
  text(mp,binding_region[,6]+0.5 ,cex = 1.3, adj = 0 ,labels=binding_region[,3] ,srt=90)
  dev.off()

}


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



RACEseq<- function(input_data, replicon_ref, mismatch = 0, RACE_adapter=NULL, str=NULL, end=NULL, filename, tmap =NULL) {


  #Main function to call all sub founctions in main RACEseqR directory R.
  source_dir <- function(path, trace = TRUE, ...) {
    for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
      if(trace) cat(nm,":")
      source(file.path(path, nm), ...)
      if(trace) cat("\n")
    }
  }

  # source_dir("R/alignment")
  # source_dir("R/scripts")

  #List of required libraries to be loaded
  suppressMessages(packages(Biostrings))
  suppressMessages(packages(tools))

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






