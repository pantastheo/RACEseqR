

#' RACEseqR plot function
#'
#' @description Function to generate and output the RACE seq plot of the specified start and end positions.
#' @param binding_region The region that will be plotted in the graph. Can be an R dataframe or a tab delim file in your working directory.
#' @param filename The filename that will be used for the output pdf file.
#' @keywords binding_region
#' @export plot_out
#' @examples
#' out_plot()

plot_out<- function(binding_region, filename) {
  
  #setting function call conditions
  if(missing(filename)) filename<- "plot"
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
