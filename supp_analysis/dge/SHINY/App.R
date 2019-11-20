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
      
      expTable_data=read.table(expTable$datapath, header=TRUE,sep="\t", row.names = 1)
      selectInput("gene1", "Select Gene", choices = rownames(expTable_data))
      
    })
  })
  
  # Render drop down menu to select from several columns
  observe({
    output$cell_type<-renderUI({
      gTable <- input$file2
      if(is.null(gTable)) return(NULL)
      
      gTable_data=read.table(gTable$datapath, header=TRUE,sep="\t", row.names = 1)
      selectInput("cell_type", "Select the class", choices = colnames(gTable_data))
      
    })
  })
  
  # Violin plot
  observe({
  output$vioPlot<-renderPlot({
    expTable <- input$file1
    gTable <- input$file2
    gene <- input$gene1
    cell_cl <- input$cell_type
    
    if(is.null(expTable)) return(NULL)
    if(gene=="") return(NULL)
    
    if(is.null(gTable)){
      expTable_data=read.table(expTable$datapath, header=TRUE,
                               sep="\t", row.names = 1)
      plot_data=data.frame("Category" = rep("All Subjects", ncol(expTable_data)),
                           "Expression" = as.numeric(expTable_data[gene,]))
    } else{
      expTable_data=read.table(expTable$datapath, header=TRUE,
                               sep="\t", row.names = 1)
      gTable_data=read.table(gTable$datapath, header=TRUE,
                             sep="\t", row.names = 1)
      
      plot_data=data.frame("Category" = gTable_data[colnames(expTable_data),cell_cl],
                           "Expression" = as.numeric(expTable_data[gene,]))
      
    }

    tocompare <- list( c("Type1", "Type2"))
    vPlot=ggviolin(plot_data, x = "Category", y = "Expression", fill = "Category",
      palette = c("#E69F00", "#56B4E9"),
      add = "boxplot", add.params = list(fill = "white"))+
      theme_classic() +
      theme(legend.position="none",plot.title = element_text(lineheight=.8, face="bold"))+ 
      scale_x_discrete(labels=c("Type1" = "Ep", "Type2" = "Tamm"))+
      stat_compare_means(comparisons = tocompare, label = "p.signif")+
      stat_compare_means(method = "t.test",label.y = max(plot_data$Expression)*1.1,size = 6)+
      xlab("")+
      ylab("Adjusted Expression")+
      theme(axis.text.x=element_text(size=14,face="bold",angle = 45, hjust = 1))+
      theme(axis.text.y=element_text(size=14,face="bold"),axis.title.y=element_text(size=16,face="bold"))

    return(vPlot)
    
  })
  })
  
  
  # Box plot
  observe({
    output$boxPlot<-renderPlot({
      expTable <- input$file1
      gTable <- input$file2
      gene <- input$gene1
      cell_cl <- input$cell_type
      
      if(is.null(expTable)) return(NULL)
      if(gene=="") return(NULL)
      
      if(is.null(gTable)){
        expTable_data=read.table(expTable$datapath, header=TRUE,
                                 sep="\t", row.names = 1)
        plot_data=data.frame("Category" = rep("All Subjects", ncol(expTable_data)),
                             "Expression" = as.numeric(expTable_data[gene,]))
      } else{
        expTable_data=read.table(expTable$datapath, header=TRUE,
                                 sep="\t", row.names = 1)
        gTable_data=read.table(gTable$datapath, header=TRUE,
                               sep="\t", row.names = 1)
        
        plot_data=data.frame("Category" = gTable_data[colnames(expTable_data),cell_cl],
                             "Expression" = as.numeric(expTable_data[gene,]))
        
      }
  
      tocompare <- list( c("Type1", "Type2"))
      bPlot=ggboxplot(plot_data, "Category", "Expression",
      color = "Category", palette =c("#E69F00", "#56B4E9"),
      add = "jitter", shape = "Category",notch = FALSE)+
      theme_classic() +
      theme(legend.position="none",plot.title = element_text(lineheight=.8, face="bold"))+ 
      scale_x_discrete(labels=c("Type1" = "Ep", "Type2" = "Tamm"))+
      stat_compare_means(comparisons = tocompare, label = "p.signif")+
      stat_compare_means(method = "t.test",label.y = max(plot_data$Expression)*1.1,size = 6)+
      xlab("")+
      ylab("Adjusted Expression")+
      theme(axis.text.x=element_text(size=14,face="bold",angle = 45, hjust = 1))+
      theme(axis.text.y=element_text(size=14,face="bold"),axis.title.y=element_text(size=16,face="bold"))

      return(bPlot)
      
    })
  })
  
}

# UI
library(shiny)

# Define UI for application
ui <- fluidPage(
  
  # Application title
  titlePanel(h1("WS DGE")),
  
  # Sidebar with a slider input for the number of bins
  sidebarLayout(
    sidebarPanel(
        # Expression Matrix
        fileInput('file1', 'Choose Expression Table',
                  accept=c('text/csv', 
                           'text/comma-separated-values,text/plain', 
                           '.csv')),
        
        # Subgroups (optional)
        fileInput('file2', 'Choose Subgroups Table',
                  accept=c('text/csv', 
                           'text/comma-separated-values,text/plain', 
                           '.csv')),
        
        # Choosing which subgroup
        uiOutput("cell_type"),
        
        # Gene name
        textInput("gene1", 'Choose a Gene', value = "", 
                  width = NULL, placeholder = "Input gene here"),
        width = 3),

    # Show plots and tables
    mainPanel(
      
        tabsetPanel(
        tabPanel("Violin Plot", plotOutput("vioPlot",width = "500px", height = "500px")), 
        tabPanel("Boxplot", plotOutput("boxPlot",width = "500px", height = "500px"))
      )
    )
  )
)

shinyApp(ui = ui, server = server)