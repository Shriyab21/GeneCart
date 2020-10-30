#' @title GeneCart
#'
#' @description This package reads the Genes file and runs the simulation against the main Gene file and convert descriptions into Monogram and Bigram and calculates PVal and -log10 val.
#'
#' @author Shriya Bansal
#' @param
#' the TEXT file loaded by the researchers.
#'
#' @return the dataframe
#
#' @examples GeneCart("C:/GeneCart/ExampleList.TXT")
#'
#' @export
  GeneCart <- function(input)
  #GeneCart <- function(Filename, Bigram_flg, Rmv_stopwrd)
{
  # Load the package required to read JSON files
  # library("rjson")
  
  library(tidytext)
  library(tibble)
  library(dplyr)
  library(wordcloud2)
  library(stopwords)
  library(data.table)
  library(tidyr)
  library(tm)
  
    Filename <- input$Wordcloud_file$datapath
    Bigram_flg <- input$wordfreq
    Rmv_stopwrd <- input$stopwrd
    
  df1 <- read.table(file = Filename,
                    #  df1 <- read.table(file = "C:/GeneCart/GeneCartApp/DataSet/Wang.txt",
                    header=FALSE)
  #change data into uppercase
  df1$V1 = toupper(df1$V1)             
  
  #Load Athaliana.text file in dataframe
  df2 <- read.table(file = "C:/GeneCart/GeneCartApp/DataSet/Athaliana.txt",
                    header=FALSE)
  
  #change data into uppercase
  df2$V1 = toupper(df2$V1)                  
  
  #Load GeneNames file
  df3 <- read.table(file = "C:/GeneCart/GeneCartApp/DataSet/Athaliana_gene_names.txt",
                    header=FALSE, sep="|",
                    na.strings = "",
                    stringsAsFactors = FALSE,
                    comment.char = "",
                    quote = "\"",
                    fill = FALSE)
  
  # Com_Gene <- intersect(df1,df2)
  Com_Gene <- left_join(df1,df2)
  #Merged the result of the intersection to the GeneNames file
  #Com_Gene <- merge(Com_Gene,df3, by.x = "V1", by.y = "V1")
  Com_Gene <- left_join(Com_Gene,df3)
  Gene_count = as.integer(nrow(Com_Gene))
  
  #Com_Gene <- select(Com_Gene, -c(V1))
  if(Rmv_stopwrd == "Y")
  {
    Com_Gene$V2 = removeWords(Com_Gene$V2, stopwords("english"))
    Com_Gene$V2 = removeNumbers(Com_Gene$V2)
    Com_Gene$V2 = removePunctuation(Com_Gene$V2)  
  }
  
  if(Bigram_flg == "Bigram")
  {
    Gene_bigrams <- Com_Gene %>% 
    unnest_tokens(bigram, V2, token = "ngrams", n = 2)
  }
  else{
    Gene_bigrams <- Com_Gene %>% 
      unnest_tokens(bigram, V2, token = "ngrams", n = 1)
    }
  
  Gene_bigrams_cnt <- Gene_bigrams %>%
    count(bigram, sort = TRUE)
  
  #Simulation_func <- function(df2, df3, Gene_count, Gene_Token){
  i <- 1
  Sim_Matrix <- setNames(data.frame(sample_n(df2,Gene_count)), toString(i))
  i <- i + 1
  while (i <= 100) {
    Temp_df <- setNames(data.frame(sample_n(df2,Gene_count)), toString(i))
    Sim_Matrix <- cbind(Sim_Matrix,Temp_df)
    i <- i + 1
  }
  
  #Melt the Matrix
  Sim_Matrix1 <- melt(Sim_Matrix, id.vars = c())
  #write.table(Sim_Matrix1, file="Sim_Matrix1.csv", sep=",")
  #Merged the result of the intersection to the GeneNames file
  Sim_genelist <- merge(Sim_Matrix1,df3, by.x = "value", by.y = "V1")
  
  #Sim_genelist <- select(Sim_genelist, -c(value))
  
  if(Rmv_stopwrd == "Y")
  {  
    Sim_genelist$V2 = removeWords(Sim_genelist$V2, stopwords("english"))
    Sim_genelist$V2 = removeNumbers(Sim_genelist$V2)
    Sim_genelist$V2 = removePunctuation(Sim_genelist$V2)
  }
  
  if(Bigram_flg == "Bigram")
  {
    Sim_Gene_bigrams <- Sim_genelist %>% 
    unnest_tokens(bigram, V2, token = "ngrams", n = 2)
  }
  else{
    Sim_Gene_bigrams <- Sim_genelist %>% 
    unnest_tokens(bigram, V2, token = "ngrams", n = 1)
  }
  Sim_bigrams_cnt <- Sim_Gene_bigrams %>%
    group_by(bigram,variable) %>%
    summarize(count=n()) %>%
    arrange(count)
  
  Sim_gene_Freq1 <- merge(Gene_bigrams_cnt, Sim_bigrams_cnt,
                          by.x = "bigram", by.y = "bigram")
  
  #Simulation gene count
  Sim_gene_Cnt <- Sim_gene_Freq1
  Sim_gene_Cnt$variable <- NULL
  
  groupword1 <- group_by(Sim_gene_Cnt, bigram,n)
  
  Sim_gene_Cnt1 <- summarise(groupword1,
                             Sim_count = sum(count))
  #clarify logic
  Sim_gene_Freq1$Comp <- Sim_gene_Freq1$n <=
    Sim_gene_Freq1$count
  
  Sim_gene_Freq1$comp1 <- ifelse(Sim_gene_Freq1$Comp == 
                                   TRUE, 1,0)
  
  Sim_gene_Freq1$count = NULL 
  Sim_gene_Freq1$Comp = NULL
  Sim_gene_Freq1$variable = NULL
  
  groupword <- group_by(Sim_gene_Freq1, bigram,n)
  
  Sim_gene_Freq2 <- summarise(groupword,
                              Freq = sum(comp1))
  
  # P-Val  
  Sim_gene_Freq2$Freq <- Sim_gene_Freq2$Freq /100   
  
  #Remove records with 0 count, else they change to "inf" after -log10
  Sim_gene_Freq2 <- Sim_gene_Freq2[Sim_gene_Freq2$Freq > 0,]
  # Sim_gene_Freq2$count <- NULL
  #Pval
  sim_gene_Pval <- Sim_gene_Freq2
  colnames(Sim_gene_Freq2)[2] <- "count"
  
  # -log10  
  Sim_gene_log10 <- Sim_gene_Freq2
  Sim_gene_log10$Freq <- -log10(Sim_gene_log10$Freq)
  
  Generate_Data <- merge(Gene_bigrams_cnt, sim_gene_Pval, 
                         by.x = "bigram", by.y = "bigram")
  colnames(Generate_Data)[4] <- "Pval"
  colnames(Generate_Data)[2] <- "Gene_Count"
  Generate_Data$n.y=NULL
  
  Generate_Data <- merge(Generate_Data,  
                         Sim_gene_log10,
                         by.x = "bigram", by.y = "bigram") 
  
  colnames(Generate_Data)[5] <- "Log10"
  Generate_Data$count=NULL
  
  Generate_Data <- merge(Generate_Data, Sim_gene_Cnt1,
                         by.x = "bigram", by.y = "bigram")
  Generate_Data$n = NULL
  
  #Generate_Data2 <- Generate_Data
  #Generate_Data2 <- Generate_Data[order(Generate_Data$Log10)]
  #mtcars[order(mpg, cyl),]
  Generate_Data2 <- Generate_Data %>% arrange(desc(Log10))
  Generate_Data2$Gene_Count=NULL
  Generate_Data2$Pval=NULL
  Generate_Data2$Sim_count=NULL
  
  Generate_Data1 <- merge(Generate_Data, Gene_bigrams, 
                          by.x = "bigram", by.y = "bigram")    
  
  #  write.table(Generate_Data, file="WangAnalysis.csv", sep=",")
  
  return_list <- list("df1" = Generate_Data2, "df2" = Generate_Data1)
  return(return_list)
}