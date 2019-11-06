library(shiny)
library(ggplot2)
library(ggpubr)

options(shiny.maxRequestSize=200*1024^2, shiny.launch.browser = TRUE)
# Define server logic
server <- function(input, output) {
  # Render drop down menu to select gene
  observe({
    output$gene1<-renderUI({
      expTable <- input$file1
      if(is.null(expTable)) return(NULL)
      
      expTable_data=read.table(expTable$datapath, header=TRUE,sep="\t")
      selectInput("gene1", "Select Gene", choices = rownames(expTable_data))
      
    })
  })
  
  # Render drop down menu to select from several columns
  observe({
    output$waves<-renderUI({
      gTable <- input$file2
      if(is.null(gTable)) return(NULL)
      
      gTable_data=read.table(gTable$datapath, header=TRUE,sep="\t")
      selectInput("waves", "Select the class", choices = rownames(gTable_data))
      
    })
  })
  
  # Violin plot
  observe({
    output$vioPlot<-renderPlot({
    expTable <- input$file1
    gTable <- input$file2
    gene <- input$gene1
    wave <- input$waves
    
    if(is.null(expTable)) return(NULL)
    if(gene=="") return(NULL)
    
    if(is.null(gTable)){
      expTable_data=read.table(expTable$datapath, header=TRUE,
                               sep="\t", row.names = 1)
      plot_data=data.frame("Category" = rep("All Subjects", ncol(expTable_data)),
                           "Expression" = as.numeric(expTable_data[gene,]))
    } else{
      expTable_data=read.table(expTable$datapath, header=TRUE,sep="\t", row.names = 1)
      gTable_data=read.table(gTable$datapath, header=TRUE,sep="\t", row.names = 1)
      
      plot_data=data.frame("Oscillation" = as.numeric(gTable_data[wave,]),
                           "Expression" = as.numeric(expTable_data[gene,]))
      
    }

    vPlot=ggstatsplot::ggscatterstats(
  data = plot_data,                                            
  x = Expression,                                                   
  y = Oscillation,                                                  
  xlab = "Adjusted Expression",                 
  ylab = "Oscillation",                                     
  point.alpha = 0.7,
  point.size = 4,
  point.color = "grey50",
  marginal = TRUE,                                             # show marginal distribution 
  marginal.type = "density",                                   # type of plot for marginal distribution
  #centrality.para = "median",                                    # centrality parameter to be plotted
  margins = "both",                                            # marginal distribution on both axes
  xfill = "#CC79A7",                                           # fill for marginals on the x-axis
  yfill = "#009E73",                                           # fill for marginals on the y-axis
  xalpha = 0.5,                                                # transparency for the x-axis marginals
  yalpha = 0.5,                                               # transparency for the y-axis marginals
  xsize = 1,                                                   # size for the x-axis marginals
  ysize = 1,                                                   # size for the y-axis marginals
  type = "spearman",                                            # type of linear association
  title = "Relationship between SME and Gene Expression",
  messages = FALSE)

    return(vPlot)
    
  })
  })
  
}

# UI
library(shiny)

# Define UI for application
ui <- fluidPage(
  
  # Application title
  titlePanel(h1("Within Subject Correlation")),
  
  # Sidebar with a slider input for the number of bins
  sidebarLayout(
    sidebarPanel(
        # Expression Matrix
        fileInput('file1', 'Choose Expression Table',
                  accept=c('text/csv', 
                           'text/comma-separated-values,text/plain', 
                           '.csv')),
        
        # Subgroups (optional)
        fileInput('file2', 'Choose SME Table',
                  accept=c('text/csv', 
                           'text/comma-separated-values,text/plain', 
                           '.csv')),
        
        # Choosing which subgroup
        textInput("waves", 'Choose a WAVE (e.g. Delta)', value = "", 
                  width = NULL, placeholder = "Input SME here"),

        # Gene name
        textInput("gene1", 'Choose a Gene (e.g. IL1RAPL2)', value = "", 
                  width = NULL, placeholder = "Input gene here"),
        width = 3),

    # Show plots and tables
    mainPanel(
      
        tabsetPanel(  
        tabPanel("Scatterplot", plotOutput("vioPlot",width = "500px", height = "500px"))
      )
    )
  )
)

shinyApp(ui = ui, server = server)