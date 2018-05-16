
#Function to download, install and load the required libraries only when needed

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
#List of required libraries to be loaded
suppressMessages(packages(Biostrings))


# description = "This is a custom R script for the downstream analysis of RACE-seq data.",
# epilogue = "Thank you for using RACE-SEQ lite"))



#' RACEseq main function
#'
#' Function to generate the alignment from a RACE sequencing experiment
#' @param str The start position of the binding region.
#' @param end The end position of the binding region.
#' @param mismatch The number of mismatches to be accepted during alignment.
#' @param RACE_adapter The RACE adapter to be trimmied using cutadapt .
#' @param tmap_opt The option to use the Tmap aligner. Defaults FALSE
#' @keywords binding_region
#' @export
#' @examples
#' RACEseq()

RACEseq <- function(str, end, mismatch, RACE_adapter, tmap_opt) {


if(!is.na(str) & !is.na(end)) {
  #input the reference sequence in .fasta format
  refname<- list.files(".", pattern ="fasta", all.files = F, full.names = F)
    if ((length(refname))==0) {
      stop("No input .fasta reference file available")
    } else if ((length(refname))>=2) {
      stop("More than one reference file")
    } else if ((length(refname))==1) {
      replicon_ref<- as.character(refname)
    }

  #input the data in .fastq or .fastq.gz format
  data_fastq<- list.files(".", pattern="fastq", all.files = F, full.names = F)
    if ((length(data_fastq))==0) {
      stop("No input .fastq file available")
    } else if ((length(data_fastq))>=2) {
      stop("More than one .fastq file")
    } else if ((length(data_fastq))==1) {
      input_data<- data_fastq
    }
}else {stop("Please input Start and End nucleotide positions \n Or type [option] -h for help")}

#reading and transforming reference sequence
nt_reference <-strsplit((toString(readBStringSet(replicon_ref))), NULL , fixed = T)
nt_reference<- data.frame(lapply(nt_reference, function(x) toupper(x)), stringsAsFactors = F)

#set output names
filename <- paste("mm", mismatch, sep = "")
out_name <- paste("read_count_", filename, sep="")

if (tmap_opt==TRUE){

  prefix<-"tmap"

  #build the index
  CMD_tmapindex<- paste("tmap index -f", replicon_ref , sep=" ")
  system(CMD_tmapindex)

  #perform alignment with tmap and read count using bedtools
  if (is.na(RACE_adapter)){

    #no adapter trimming
    print(paste0("Performing alignment with ", mismatch, " mismatch using tmap"))
    CMD_tmap<- paste("tmap map1 -a 0 -g 3 --max-mismatches ",mismatch," -f ", replicon_ref," -r ", input_data, " | samtools view -bt ", replicon_ref," - | genomeCoverageBed -d -5 -ibam stdin > ",out_name, sep="")
  system(CMD_tmap)
  } else {

    #adapter trimming using cutadapt
    print(paste0("Performing adapter trimming and alignment with ", mismatch, " mismatch using tmap"))
    CMD_tmap<- paste("cutadapt -g ", RACE_adapter, " -e0 --no-indels -m10 --discard-untrimmed --quiet ", input_data," |tmap map1 -a 0 -g 3 --max-mismatches ",mismatch," -f ", replicon_ref," -i fastq | samtools view -bt ", replicon_ref," - | genomeCoverageBed -d -5 -ibam stdin > ",out_name, sep="")
    system(CMD_tmap)
  }
} else {

  prefix<-"bowtie"

  #build the index
  CMD_bowindex<- paste("bowtie-build -q -f", replicon_ref, "index", sep=" ")
  system(CMD_bowindex)

  #perform alignment with bowtie and read count using bedtools
  if (is.na(RACE_adapter)){

    #no adapter trimming
    print(paste0("Performing alignment with ", mismatch, " mismatch using bowtie"))
    CMD_bow<- paste("bowtie -p 2 -S -k 1 -v", mismatch, "index", input_data," | samtools view -bS - | genomeCoverageBed -d -5 -ibam stdin >", out_name, sep=" ")
  system(CMD_bow)

  } else {

    #adapter trimming using cutadapt
    print(paste0("Performing adapter trimming and alignment with ", mismatch, " mismatch using bowtie"))
    CMD_bow<- paste("cutadapt -g", RACE_adapter, "-e0 --no-indels -m10 --discard-untrimmed --quiet ", input_data,"|bowtie -p 8 -S -k 1 -v", mismatch, "index - | samtools view -bS - | genomeCoverageBed -d -5 -ibam stdin >", out_name, sep=" ")
    system(CMD_bow)
  }
}

#read and merge ref and reads
reads<- read.delim(out_name, header = F )
dataframe<- data.frame(reads, nt_reference , stringsAsFactors = F)

#calculating the % and log10 columns
dataframe[,5] <- (dataframe[,3]/sum(dataframe[,3])*100)
dataframe[,6] <- (log10(dataframe[,3]))
dataframe[dataframe== -Inf] <-0

#focusing on target region can be ajusted acording to experiment
binding_region <- dataframe[str:end,]

#remove read_count and index file created



#' File delete function
#'
#' Function to delete files generated in the process
#' @param pattern The filename extension to be deleted
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

counts<- list.files(".", pattern="read_count", all.files = F, full.names = F)
for(i in counts ){
  i<- paste("rm", i , sep=" ")
  system(i)
}
index_files<- list.files(".", pattern=".fasta.tmap.", all.files = F, full.names = F)
for(i in index_files ){
  i<- paste("rm", i , sep=" ")
  system(i)
}
index_files<- list.files(".", pattern="aligned.bam", all.files = F, full.names = F)
for(i in index_files ){
  i<- paste("rm", i , sep=" ")
  system(i)
}
index_files<- list.files(".", pattern="out.sam", all.files = F, full.names = F)
for(i in index_files ){
  i<- paste("rm", i , sep=" ")
  system(i)
}
index_files<- list.files(".", pattern="index", all.files = F, full.names = F)
for(i in index_files ){
  i<- paste("rm", i , sep=" ")
  system(i)
}

return(binding_region)
}

#remove read_count and index file created

#' File delete function
#'
#' Function to delete files generated in the process
#' @param pattern The filename extension to be deleted
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



#' RACEseq dataframe function
#'
#' Function to generate and output the RACE seq dataframe of the specified start and end positions.
#' @param binding_region The region that will be written in the csv file. Defaults NULL.
#' @param filename The filename of the dataset
#' @keywords binding_region
#' @export
#' @examples
#' out_csv()

out_csv<- function(binding_region, filename){
  write.table(binding_region, file = paste0(filename, ".csv") , sep = "\t", col.names = c("reference", "position", "count", "nucleotide", "percentage", "log10" ), row.names = F )
}




#' RACEseq plot function
#'
#' Function to generate and output the RACE seq plot of the specified start and end positions.
#' @param binding_region The region that will be plotted in the graph. Defaults NULL.
#' @param filename The filename of the dataset
#' @keywords binding_region
#' @export
#' @examples
#' out_plot()

out_plot<- function(binding_region, filename){

  #create wildtype linear & log scale graph
  pdf(paste0(filename, ".pdf"), width=15)

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
