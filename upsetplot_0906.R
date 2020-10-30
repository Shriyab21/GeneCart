RunupsetplotF <- function(input)
{
library(data.table)
library(dplyr)

universe_filename <- input$upsetplot_universegeneids
geneids_filename_1 <- input$upsetplot_geneids_1
geneids_filename_2 <- input$upsetplot_geneids_2
geneids_filename_3 <- input$upsetplot_geneids_3
df1 <- read.table(file = universe_filename$datapath, header=FALSE)
df1$V1 = toupper(df1$V1)   
df2 <- read.table(file = geneids_filename_1$datapath, header=FALSE)
df2$V1 = toupper(df2$V1)
if (!is.null(geneids_filename_2$datapath)){
    df3 <- read.table(file = geneids_filename_2$datapath, header=FALSE)
    df3$V1 = toupper(df3$V1)
} 
if (!is.null(geneids_filename_3$datapath)){
  df4 <- read.table(file = geneids_filename_3$datapath, header=FALSE)
  df4$V1 = toupper(df4$V1)
}

plot_list <- list(a=as.vector(df1$V1),
                  b=as.vector(df2$V1),
                  if (!is.null(geneids_filename_2$datapath)){
                    c=as.vector(df3$V1)},
                  if (!is.null(geneids_filename_3$datapath)){
                    d=as.vector(df4$V1)})
                                    
return(plot_list[lapply(plot_list,length)>0])
}  

RunupsetplotT <- function(input)
{
  library(data.table)
  library(dplyr)
  
  plot_list <- list(a=as.vector(unlist(strsplit(input$list1, "\n"))), 
                    b=as.vector(unlist(strsplit(input$list2, "\n"))),
                    c=as.vector(unlist(strsplit(input$list3, "\n"))),
                    d=as.vector(unlist(strsplit(input$list4, "\n"))))
  
  return(plot_list[lapply(plot_list,length)>0])
}  
