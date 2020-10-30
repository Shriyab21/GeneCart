
#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

# In order to run Genecart, please make sure the following packages are
# installed:

if (!requireNamespace("shiny", quietly = TRUE))
    install.packages("shiny")
if (!requireNamespace("DT", quietly = TRUE))
    install.packages("DT")
if (!requireNamespace("shinydashboard", quietly = TRUE))
    install.packages("shinydashboard")
if (!requireNamespace("RSQLite", quietly = TRUE))
    install.packages("RSQLite")
if (!requireNamespace("rjson", quietly = TRUE))
    install.packages("rjson")
if (!requireNamespace("readr", quietly = TRUE))
    install.packages("readr")
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
if (!requireNamespace("UpSetR", quietly = TRUE))
    install.packages("UpSetR")


library(shiny)
library(DT)
library(shinydashboard)
library(RSQLite)
library(rjson)
library(readr)
library(wordcloud2)
library(webshot)
library(data.table)
library(dplyr)
#library(UpSetR)
library(upsetjs)
#library(GeneCart)

source("C:/TestApp/tools/annotation.R")
source("C:/TestApp/tools/biomaps.R")
source("C:/TestApp/tools/genecartdb.R")
source("C:/TestApp/tools/genecart_sqlite.R")
source("C:/TestApp/tools/genenetwork.R")
source("C:/TestApp/tools/genesect.R")
source("C:/TestApp/tools/getspecies.R")
source("C:/TestApp/tools/heatmap.R")
#source("c:/TestApp/tools/GeneCart.R")
source("c:/TestApp/tools/GeneCart_1.R")
#source("c:/TestApp/tools/upsetplot.R")
#source("c:/TestApp/tools/upsetplot_0815.R")
#source("c:/TestApp/tools/upsetplot_0824.R")
#source("c:/TestApp/tools/upsetplot_0826.R")
source("c:/TestApp/tools/upsetplot_0906.R")

#library(rcytoscapejs)

# Define UI for application 
ui <- dashboardPage(
    
    ############################ Dashboard Header #########################
    dashboardHeader(
        title="Gene Cart"),
    ############################ Dashboard Sidebar ########################
    
    dashboardSidebar(
        sidebarMenu(
            
            menuItem("GeneCart", 
                     tabName = "genecart", 
                     icon=icon("cart-plus")
            ), 
            menuItem("Gene Lists", 
                     tabName = "tools", 
                     icon=icon("stats",
                               lib="glyphicon")
            ),
            menuItem("Experiments", 
                     tabName = "experiments", 
                     icon=icon("braille")
            ),
            menuItem("Gene Networks", 
                     tabName = "networks", 
                     icon=icon("code-branch")
            )
        ) 
    ),
    ############################ Dashboard Body #########################
    
    dashboardBody(
        fluidRow(
            # A static infoBox
            infoBox("Species"),
            # Dynamic infoBoxes
            infoBoxOutput("Species")
        ),
        tabItems(
            
            #  ),
            
            ########################### Tab Item - Tools ########################
            tabItem(tabName = "tools",
                    h2("Gene Lists"),
                    # A panel to place all the tabs
                    tabsetPanel(type="tabs", 
                                
                                ############# QUERY ########################
                                tabPanel("Query",value = "query" ,
                                         # Within a tab panel we can have 
                                         # an entire layout
                                         # Sidebar layout includes 
                                         # a sidebar panel and a mainpanel
                                         
                                         sidebarLayout(
                                             # Where the inputs are accepted
                                             sidebarPanel(
                                                 selectInput("species", 
                                                             "Choose a species:",
                                                             choices = c("Arabidopsis",
                                                                         "Rice", 
                                                                         "Maize")),
                                                 selectInput("feature", 
                                                             "What do you want to query?",
                                                             choices = c("gene", 
                                                                         "description")),
                                                 textInput("queryterm",
                                                           label="Query:"),
                                                 actionButton("Query", 
                                                              label="Submit", 
                                                              icon("submit"))                                               
                                             ),
                                             #Where the output will be displayed
                                             #mainPanel = tableOutput("queryresult")
                                             mainPanel = dataTableOutput("queryresult")
                                         )
                                ),
                                
                                ########## BIOMAPS #########################
                                
                                tabPanel("Biomaps", value="biomaps",
                                         
                                         # Within a tab panel we can have 
                                         # an entire layout
                                         # Sidebar layout includes a sidebar panel and
                                         # a mainpanel
                                         
                                         sidebarLayout(
                                             # Where the inputs are accepted
                                             sidebarPanel(
                                                 selectInput("bpspecies", 
                                                             "Choose a species:",
                                                             choices = c("Arabidopsis", 
                                                                         "Rice", 
                                                                         "Maize")),
                                                 fileInput("geneids", 
                                                           label="Gene Ids"),
                                                 fileInput("universegeneids",
                                                           label="Universe Gene Ids"),
                                                 selectInput("ont", "Choose GO namespace",
                                                             choices = c("BP", 
                                                                         "MF", 
                                                                         "CC")),
                                                 selectInput("golevel", 
                                                             "Restrict output to GO-term level",
                                                             choices=c("All",
                                                                       "1",
                                                                       "2",
                                                                       "3")
                                                 ),
                                                 selectInput("pvalcutoff", 
                                                             "Choose adj P-value cutoff",
                                                             choices = c(1,
                                                                         0.05, 
                                                                         0.01, 
                                                                         0.001)),
                                                 actionButton("Biomaps", 
                                                              label="Submit", 
                                                              icon("submit"))
                                             ),
                                             
                                             #Where the output will be displayed
                                             #mainPanel = tableOutput("biomapsresult")
                                             mainPanel = dataTableOutput("biomapsresult")
                                         )),
                                
                                
                                ##### GENESECT #########################
                                
                                tabPanel("GeneSect", value="genesect",
                                         # Within a tab panel we can have 
                                         # an entire layout
                                         # Sidebar layout includes a sidebar panel 
                                         # and a mainpanel
                                         sidebarLayout(
                                             
                                             # Where the inputs are accepted
                                             sidebarPanel(
                                                 selectInput("gsspecies", 
                                                             "Choose a species:",
                                                             choices = c("Arabidopsis",
                                                                         "Rice", 
                                                                         "Maize")),
                                                 fileInput("genesect_geneids_1", 
                                                           label="Gene List 1"),
                                                 fileInput("genesect_geneids_2", 
                                                           label="Gene List 2"),
                                                 fileInput("genesect_universegeneids",
                                                           label="Universe Gene Ids"),
                                                 actionButton("Genesect", label="Submit",
                                                              icon("submit"))
                                             ),
                                             
                                             #Where the output will be displayed
                                             mainPanel = tableOutput("genesectresult")
                                         )),
                                
                                ##### WORDCLOUD #########################
                                
                                tabPanel("Wordcloud", value="wordcloud",
                                         # Within a tab panel we can have 
                                         # an entire layout
                                         # Sidebar layout includes a sidebar panel 
                                         # and a mainpanel
                                         sidebarLayout(
                                             
                                             # Where the inputs are accepted
                                             sidebarPanel(
                                                 fileInput("Wordcloud_file", 
                                                           label="Select a file"),
                                                 selectInput("wordfreq", 
                                                             "Choose Word frequency:",
                                                             choices = c("Bigram",
                                                                         "Monogram")),
                                                 selectInput("stopwrd", 
                                                             "Remove Stopwords:",
                                                             choices = c("Y",
                                                                         "N")),
                                                 
                                                 actionButton("Wordcloud", label="Submit",
                                                              icon("submit"))
                                             ),
                                             
                                             #Where the output will be displayed
                                             mainPanel (
                                                 p(" "),
                                                 wordcloud2Output("cloud"),
                                                 downloadButton(outputId = "savecloud"),
                                                 p(" "),
                                                 dataTableOutput("table"),
                                                 downloadButton(outputId = "savetable")
                                             ) #mainPanel
                                         ) #sidebar
                                ), #tabPanel
                                ##### UPSET PLOT #########################
                                tabPanel("Upset Plot", value="upsetplot",
                                         # Sidebar layout includes a sidebar panel 
                                         # and a mainpanel
                                         sidebarLayout(
                                             # Where the inputs are accepted
                                             sidebarPanel(
                                                 selectInput("Upset_Input_Type", 
                                                             "Choose Input type:",
                                                             choices = c("File",
                                                                         "Text")),
                                                 conditionalPanel(
                                                     condition = "input.Upset_Input_Type == 'File'",
                                                     
                                                     fileInput("upsetplot_universegeneids",
                                                               label="Universe Gene Ids"),
                                                     fileInput("upsetplot_geneids_1", 
                                                               label="Gene List 1"),
                                                     conditionalPanel(
                                                         condition = "input.addinput1 >= 1",
                                                         fileInput("upsetplot_geneids_2", 
                                                                   label="Gene List 2")),
                                                     conditionalPanel(
                                                         condition = "input.addinput1 >= 2",
                                                         fileInput("upsetplot_geneids_3", 
                                                                   label="Gene List 3")),
                                                     actionButton("addinput1","Add another gene set"),
                                                     actionButton("upsetplotF", label="Run Upset",
                                                                  icon("submit"))),
                                                 
                                                 conditionalPanel(
                                                     condition = "input.Upset_Input_Type == 'Text'",
                                                     textAreaInput("list1", "Please input the first list"),
                                                     textAreaInput("list2", "Please input the second list"),
                                                     conditionalPanel(
                                                         condition = "input.addinput >= 1",
                                                         textAreaInput("list3", "Please input the third list")
                                                     ),
                                                     conditionalPanel(
                                                         condition = "input.addinput >= 2",
                                                         textAreaInput("list4", "Please input the forth list")
                                                     ),
                                                     actionButton("addinput","Add another gene set"),
                                                     
                                                     actionButton("upsetplotT", label="Run Upset",
                                                                  icon("submit")))
                                             ),
                                             mainPanel (
                                                 p(" "), 
                                                 conditionalPanel(
                                                     condition = "input.Upset_Input_Type == 'File'",
                                                     upsetjsOutput("upsetplotresult")),
                                                 conditionalPanel(
                                                     condition = "input.Upset_Input_Type == 'Text'",
                                                     upsetjsOutput("upsetplotresult1")),
                                                 #plotOutput("upsetplotresult")
                                                 #          hover = "plot_hover",
                                                 #  hover = "plot_hover" , click = "plot_click"),
                                                 p(" "),
                                                 #           verbatimTextOutput("plot_info")
                                                 
                                             )  #mainPanel
                                         )#sidebar
                                ) #tabPanel
                    ) #tabsetPanel
            ), #tab item
            ########################### Tab Item - Genecart ########################
            
            tabItem(tabName="genecart", 
                    h2("Gene Cart")),
            
            tabItem(tabName="experiments", 
                    h2("Experiments"),
                    tabsetPanel(type="tabs",
                                tabPanel("Experiments", value="HeatMap",
                                         # Within a tab panel we can have an entire layout
                                         # Sidebar layout includes a sidebar panel and a mainpanel
                                         sidebarLayout(
                                             
                                             # Where the inputs are accepted
                                             sidebarPanel(
                                                 fileInput("expgeneids", label="Gene Ids"),
                                                 fileInput("expvalues", label="Experiment Values"),
                                                 actionButton("Heatmap", label="Submit", icon("submit"))
                                             ),
                                             
                                             #Where the output will be displayed
                                             mainPanel = plotOutput("heatmapresult")
                                         )
                                ) #tabPanel
                    ) #tabsetpanel
                    
                    
            ),
            
            tabItem(tabName="networks", 
                    h2("Gene Networks"),
                    ############ GENENETWORK #########################
                    tabsetPanel(type="tabs",
                                tabPanel("GeneNetwork", value="genenetwork",
                                         # Within a tab panel we can have an entire layout
                                         # Sidebar layout includes a sidebar panel and a mainpanel
                                         sidebarLayout(
                                             
                                             # Where the inputs are accepted
                                             sidebarPanel(
                                                 selectInput("gnspecies", "Choose a species:",
                                                             choices = c("Arabidopsis", "Rice", "Maize")),
                                                 fileInput("gngeneids", label="Gene Ids"),
                                                 actionButton("Genenetwork", label="Submit", icon("submit"))
                                             ),
                                             
                                             #Where the output will be displayed
                                             mainPanel = tableOutput("genenetworkresult")
                                         )
                                ) #tabPanel
                    ) #tabsetpanel
            ) #tab item
        ) #tab items
    ) #dashboard body
) #ui

#########################SERVER ########################################
server <- function(input, output, session) {
    output$progressBox <- renderInfoBox({
        infoBox(
            input$species
            
        )
    })
    print(str(input))
    ############################ QUERY #########################
    queryevent = eventReactive(input$Query, {
        as.data.frame(getAnnotation(input))
    })
    
    output$queryresult = renderDataTable( queryevent()  )
    
    ############################ BIOMAPS #########################
    biomapsevent = eventReactive(input$Biomaps, {
        as.data.frame(runBiomaps(input))
    })
    
    output$biomapsresult = renderDataTable( biomapsevent()  )
    ############################ GENESECT #########################
    genesectevent = eventReactive(input$Genesect, {
        as.data.frame(runGenesect(input))
    })
    ############################ WORDCLOUD-func #########################
    wordcloudevent = eventReactive(input$Wordcloud, {
        data_source <- ({
            if (is.null(input$Wordcloud_file$datapath)) {
                return(NULL)
            }
            ((GeneCart(input)))
        })
        return(data_source)
    })
    ############################ WORDCLOUD-Cloud #########################
    wordcloudevent1 = eventReactive(input$Wordcloud, {
        req(wordcloudevent())
        data_source <- wordcloudevent()
        
        wordcloud2(data.frame(data_source$df1), size=.2,
                   backgroundColor = "White")
    })
    output$cloud = renderWordcloud2(wordcloudevent1())
    
    output$savecloud <- downloadHandler(
        filename = paste("wordcloud", '.png', sep=''),
        content = function(file) {
            owd <- setwd(tempdir())
            on.exit(setwd(owd))
            saveWidget(wordcloudevent1(), "temp.html", selfcontained = FALSE)
            webshot("temp.html", delay =15, file = file, cliprect = "viewport")
        })
    ############################ WORDCLOUD-Table #########################
    wordcloudevent2 = eventReactive(input$Wordcloud, {
        req(wordcloudevent())
        data_source <- wordcloudevent()
        as.data.frame(data_source$df2)
    })
    output$table = renderDataTable(wordcloudevent2())
    
    output$savetable <- downloadHandler(
        filename = paste("wordcloudtable", '.csv', sep=''),
        content = function(file) {
            owd <- setwd(tempdir())
            on.exit(setwd(owd))
            saveWidget(wordcloudevent2(), "temp.csv", selfcontained = FALSE)
            webshot("temp.csv", delay =15, file = file, cliprect = "viewport")
        }) 
    ############################ UPSET Function - File #########################
    upsetplotevent = eventReactive(input$upsetplotF, {
        upset_results <- as.list(RunupsetplotF(input))
        return(upset_results)             
    })
    ############################ UPSET Function - Plot (File)#########################
    upsetplotevent1 = eventReactive(input$upsetplotF, {
        req(upsetplotevent())
        upset_results <- upsetplotevent()
        upsetjs() %>% upsetjs::fromList(upset_results) %>% interactiveChart()
    })
    output$upsetplotresult = renderUpsetjs(upsetplotevent1())
    
    ############################ UPSET Function - Text #########################
    upsetplotevent2 = eventReactive(input$upsetplotT, {
        upset_results <- as.list(RunupsetplotT(input))
        return(upset_results)             
    })
    ############################ UPSET Function - Plot (File)#########################
    upsetplotevent3 = eventReactive(input$upsetplotT, {
        req(upsetplotevent2())
        upset_results <- upsetplotevent2()
        upsetjs() %>% upsetjs::fromList(upset_results) %>% interactiveChart()
    })
    output$upsetplotresult1 = renderUpsetjs(upsetplotevent3())
    
    ############################ GENENETWORK #########################
    genenetworkevent = eventReactive(input$Genenetwork, {
        as.data.frame(runGenenetwork(input))
    })
    output$genenetworkresult1 = renderTable( genenetworkevent()  )
    ############################ HEATMAP #########################
    heatmapevent = eventReactive(input$Heatmap, {
        as.data.frame(runHeatmap(input))
    })
    
    output$heatmapresult = renderPlot( heatmapevent()  )
    
}

# Run the application 
shinyApp(ui = ui, server = server)

