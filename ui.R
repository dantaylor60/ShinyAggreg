ui <- shinyUI(fluidPage(
  headerPanel('NEON data aggregation by species'),
  sidebarPanel(
    selectInput('group','data type',c('smallMammal','carabid','plantPresenceCover')),
    selectInput('space.Agg','spatial scale',c('byDomain','bySite','byPlot')),
    selectInput('time.Agg','temporal scale',c("aggAll",'byYear','byMonthYear','byBout')),
    selectInput('typeFiles','file type for plant presence',c("none",'1m2Data','10m2Data','100m2Data','400m2Data')),
    numericInput('numcols', label = "number of columns to display",value=4),
    h3("Download data"),
    downloadButton('downloadData', 'Download')
    ),
  
  mainPanel(
    tabsetPanel(
      tabPanel('preview data',
               dataTableOutput("aggData")),
      tabPanel('summary',
               dataTableOutput("summary"))
      )
    )
  
))

#source("server.R")
#shinyApp(ui = ui, server = serverAgg)