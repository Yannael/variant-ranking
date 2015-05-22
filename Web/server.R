library(shiny)
requireNamespace('htmlwidgets')
library(DT)
library(d3heatmap)

shinyServer(function(input, output,session) {
  sessionvalues <- reactiveValues()
  sessionvalues$phenotypes<-phenotypesAll
  
  sessionvalues$currentResults<-res
  
  getSampleIDFromGroup<-function(groupName) {
    index.group<-which(mySampleGroups[[1]]==paste0(groupName))
    mySampleGroups[[2]][[index.group]]
  }
  
  observe({
    selectedSampleGroup<-input$selectedSampleGroup
    #browser()
    sessionvalues$phenotypes<-
      switch(selectedSampleGroup,
             all=phenotypesAll,
             erasme=phenotypesErasme,
             genomes1000=phenotypes1000Gen,
             phenotypesAll[which(phenotypesAll[,'Sample ID'] %in% getSampleIDFromGroup(selectedSampleGroup)),]
      )
  })
  
  output$phenotypesTable<-DT::renderDataTable({
    data<-as.data.frame(sessionvalues$phenotypes[,input$showVarPhenotype])
    #data},
    datatable(
      data,escape= -1, rownames=F
      ,filter='top', 
      options = list(
        dom='fltip', 
        lengthMenu = list(c(10, 25, -1), c('10', '25','All')),pageLength = 10,
        autoWidth = T,columnDefs = list(list(sClass="alignRight",aTargets="_all"))
      )
    )
  }
  )
  
  output$resultsTable2<-DT::renderDataTable({
    
    data<-as.data.frame(sessionvalues$currentResults$matRes)
    action = dataTableAjax(session, data)
    datatable(data, 
              #rownames = checkboxRows(data),
              extensions = c('ColVis','Scroller'),
              escape=F,
              filter = 'top',
              options = list(
                #sscrollX = TRUE,
                #colVis=list(exclude=c(0)),
                dom= 'CrtiS',
                ajax = list(url = action),
                deferRender = TRUE,
                scrollY = 500,
                scrollCollapse = T,
                lengthMenu = list(c(10, 25, -1), c('10', '25','All')),pageLength = 10,
                autoWidth = T,
                columnDefs = list(list(className="dt-right",targets="_all"),
                                  list(visible=F,targets=c(4,6:7,9:14)) 
                )
              )
    )
  })
  
  output$resultsTable<-DT::renderDataTable({
    data<-sessionvalues$currentResults$matRes
    if (!is.null(dim(data))) {
      action = dataTableAjax(session, data)
      
      widget = datatable(data, 
                         server = TRUE, 
                         extensions = c('ColVis','Scroller'),
                         escape=F,
                         #rownames=F,
                         #filter = 'top',
                         options = list(
                           dom= 'CfrtSi',
                           scrollX = TRUE,
                           ajax = list(url = action),
                           deferRender = TRUE,
                           scrollY = 335,
                           scrollCollapse = T,
                           colVis=list(exclude=c(0)),
                           autoWidth = F,
                           columnDefs = list(list(className="dt-right",targets="_all"), 
                                             list(width='100px',targets="_all"),
                                             list(visible=F,targets=c(0,4,6:7,10:15)) 
                                             
                           )
                         )
      )
      widget
    }
  })
  
  output$heatMapGenotypes<-renderD3heatmap({
    #data<-input$resultsTable_rows_all
    subMat<-sessionvalues$currentResults$snpsMat[1:10,]
    rownames(subMat)<-sessionvalues$currentResults$matRes$Locus[1:10]
    
    patient_class<-colnames(subMat)
    patient_class<-sapply(patient_class,strsplit,':')
    patient_class<-unlist(patient_class,use.name=F)
    patient_class<-matrix(patient_class,ncol(subMat),2,byrow=T)
    
    controls<-which(patient_class[,2]=="Control")
    cases<-which(patient_class[,2]=="Patho")
    
    i.order<-sort(patient_class[,2],index=T)$ix
    subMat<-subMat[,i.order]
    
    d3heatmap(subMat,cluster=F)
  })
  
  output$resultsMetadata<-renderUI({
    fluidRow(
      column(3,
             strong("Analysis name: "),br(),
             strong("Control group: "),br(),
             strong("Pathological group: "),br(),
             strong("Start time: "),br(),
             strong("End time: ")
      ),
      column(4,
             sessionvalues$currentResults$name[[1]][1],br(),
             sessionvalues$currentResults$metadata$controlgroup,br(),
             sessionvalues$currentResults$metadata$pathogroup,br(),
             sessionvalues$currentResults$metadata$timestart,br(),
             sessionvalues$currentResults$metadata$timeend
      )
      
    )
  })
  
  #observe({
    #input$refreshHeatMaps
    #cdata<-input$resultsTable_rows_current
  #  browser()
  #  print("hey"
  #        )
  #})
  
  output$resultsPanel<-renderUI({
    fluidRow(
      column(10,offset=1,
             uiOutput("resultsMetadata"),
             hr(),
             h3("Ranking results"),
             DT::dataTableOutput('resultsTable'),
             hr(),
             h3("Genotypes SNPs/samples"),
             actionButton("refreshHeatMaps","Refresh"),
             d3heatmapOutput("heatMapGenotypes")
             #DT::dataTableOutput('tt')
      )
    )
  })
  
  
})

