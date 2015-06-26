library(shiny)
requireNamespace('htmlwidgets')
library(DT)
library(d3heatmap)
library(plyr)
library(ggplot2)

shinyServer(function(input, output,session) {
  sessionvalues <- reactiveValues()
  sessionvalues$phenotypes<-phenotypesAll
  
  #sessionvalues$currentResults<-results['Trios_Compound_Heterozygous'][[1]][[1]]
  
  getSampleIDFromGroup<-function(groupName) {
    index.group<-which(mySampleGroups[[1]]==paste0(groupName))
    mySampleGroups[[2]][[index.group]]
  }
  
  observe({
    sessionvalues$currentResults<-results[input$selectedResultGroup][[1]][[1]]
  })
  
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
  
  output$resultsTable<-DT::renderDataTable({
    data<-sessionvalues$currentResults$scoreSummary
    if (!is.null(dim(data))) {
      widget = datatable(data, 
                         server = FALSE, 
                         escape=F,
                         selection = 'single',
                         filter = 'top',
                         rownames=T,
                         options = list(
                           dom= 'liptlpi',
                           lengthMenu = list(c(10, 25, 100), c('10', '25','100')),pageLength = 10,
                           autoWidth = F,
                           columnDefs = list(
                             list(className="dt-right",targets="_all")
                           )
                         )
      )
      widget
    }
  })
  
  output$variantsMetadataTable<-DT::renderDataTable({
    if (length(input$resultsTable_rows_selected)>0) {
      res<-NULL
      if (sessionvalues$currentResults$type=="singleVariant") {
        infovariants<-sessionvalues$currentResults$infovariants[input$resultsTable_rows_selected,]
        data<-cbind(names(infovariants),t(infovariants))
        #browser()
        colnames(data)<-c("Variable","Value")
        #action = dataTableAjax(session, data)
        if (!is.null(dim(data))) {
          widget = datatable(data, 
                             server = FALSE, 
                             extensions = c('Scroller'),
                             escape=F,
                             selection = 'single',
                             filter = 'bottom',
                             rownames=F,
                             options = list(
                               dom= 'tS',
                               deferRender = TRUE,
                               scrollY = 335,
                               #ajax = list(url = action),
                               #colVis=list(exclude=c(0)),
                               scrollCollapse = T,
                               #lengthMenu = list(c(10, 25, 100), c('10', '25','100')),pageLength = 10,
                               autoWidth = F,
                               columnDefs = list(
                                 list(className="dt-right",targets="_all")
                               )
                             )
          )
          res<-widget
        }
      }
      if (sessionvalues$currentResults$type=="pairVariantsMonogenic") {
        infovariants1<-sessionvalues$currentResults$infovariants1[input$resultsTable_rows_selected,]
        infovariants2<-sessionvalues$currentResults$infovariants2[input$resultsTable_rows_selected,]
        data<-cbind(names(infovariants1),t(infovariants1),t(infovariants2))
        #browser()
        data[is.na(data)]<-"NA"
        colnames(data)<-c("Variable","Value1","Value2")
        action = dataTableAjax(session, data,rownames=F)
        if (!is.null(dim(data))) {
          widget = datatable(data, 
                             server = TRUE, 
                             #extensions = c('Scroller'),
                             selection = 'single',
                             rownames=F,
                             filter = 'bottom',
                             escape=T,
                             options = list(
                               dom= 't',
                               #lengthMenu = list(c(100), c('100')),
                               pageLength = 100,
                               #deferRender = FALSE,
                               scrollY = 335,
                               ajax = list(url = action),
                               columnDefs = list(
                                 list(
                                   targets = c(1,2),
                                   render = JS(
                                     "function(data, type, row, meta) {",
                                     "return type === 'display' && data.length > 15 ?",
                                     "'<span title=\"' + data + '\">' + data.substr(0, 11) + '...</span>' : data;",
                                     "}")
                                 ),
                                 #list(width=c(59),targets=c(1)),
                                 #list(visible=F,targets=c(0)),
                                 list(className="dt-right",targets="_all")
                               )
                             )
          )
          res<-widget
        }
      }
      res
    }
  })
  
  
  output$heatMapGenotypes<-renderPlot({
    if (length(input$resultsTable_rows_selected)>0) {
      res<-NULL
      if (sessionvalues$currentResults$type=="singleVariant") {
        genotypes<-sessionvalues$currentResults$genotypes[input$resultsTable_rows_selected,]
        genotypes<-factor(genotypes,levels=c("0","1","2","NA"))
        genotypes<-mapvalues(genotypes,from=c(0,1,2,"NA"),to=c("Ref/Ref","Ref/Alt","Alt/Alt","NA"))
        
        DF<-data.frame(rep(paste(rep("Trio ",11),1:11,sep=""),each=3),rep(c("Child","Father","Mother"),11),genotypes)
        DF[,1]<-factor(DF[,1],levels=rev(unique(DF[,1])),ordered=T)
        colnames(DF)<-c("Trios","FMC","Zygosity")
        
        base_size=20
        p<-ggplot(DF, aes(y=Trios,x=FMC))
        p <- p + geom_tile(aes(fill=Zygosity),colour = "white")
        p <- p +  theme_grey(base_size = base_size) + labs(x = "", y = "") + 
          scale_x_discrete(expand = c(0, 0)) +
          scale_y_discrete(expand = c(0, 0)) +
          theme(aspect.ratio=3)+
          theme(axis.title.x = element_blank(),
                axis.text.x = element_text(size = base_size *  0.8, angle = 330,hjust = 0, colour = "grey50")) 
        p<-p+ scale_fill_manual(breaks=c("Ref/Ref","Ref/Alt","Alt/Alt","NA"),
                                values=c("#fff7bc","#fec44f","#d95f0e","000000"),
                                #values=c("#ece7f2","#a6bddb","#2b8cbe","000000"),
                                name="Zigosity")
        
        p
        res<-p
      }
      if (sessionvalues$currentResults$type=="pairVariantsMonogenic") {
        genotypes1<-sessionvalues$currentResults$genotypes1[input$resultsTable_rows_selected,]
        genotypes1<-factor(genotypes1,levels=c("0","1","2","NA"))
        genotypes1<<-mapvalues(genotypes1,from=c(0,1,2,"NA"),to=c("Ref/Ref","Ref/Alt","Alt/Alt","NA"))
        genotypes2<-sessionvalues$currentResults$genotypes2[input$resultsTable_rows_selected,]
        genotypes2<-factor(genotypes2,levels=c("0","1","2","NA"))
        genotypes2<<-mapvalues(genotypes2,from=c(0,1,2,"NA"),to=c("Ref/Ref","Ref/Alt","Alt/Alt","NA"))
        
        DF1<-data.frame(rep(paste(rep("Trio ",11),1:11,sep=""),each=3),rep(c("Child variant 1","Father variant 1","Mother variant 1"),11),genotypes1)
        colnames(DF1)<-c("Trios","FMC","Zygosity")
        DF2<-data.frame(rep(paste(rep("Trio ",11),1:11,sep=""),each=3),rep(c("Child variant 2","Father variant 2","Mother variant 2"),11),genotypes2)
        colnames(DF2)<-c("Trios","FMC","Zygosity")
        DF<-rbind(DF1,DF2)
        
        DF[,1]<-factor(DF[,1],levels=rev(unique(DF[,1])),ordered=T)
        DF[,2]<-factor(DF[,2],levels=sort(levels(DF[,2])),ordered=T)
        
        base_size=20
        p<-ggplot(DF, aes(y=Trios,x=FMC))
        p <- p + geom_tile(aes(fill=Zygosity),colour = "white")
        p <- p +  theme_grey(base_size = base_size) + labs(x = "", y = "") + 
          scale_x_discrete(expand = c(0, 0)) +
          scale_y_discrete(expand = c(0, 0)) +
          theme(aspect.ratio=1.5)+
          theme(axis.title.x = element_blank(),
                axis.text.x = element_text(size = base_size *  0.8, angle = 310,hjust = 0, colour = "grey50")) 
        p<-p+ scale_fill_manual(breaks=c("Ref/Ref","Ref/Alt","Alt/Alt","NA"),
                                values=c("#fff7bc","#fec44f","#d95f0e","000000"),
                                #values=c("#ece7f2","#a6bddb","#2b8cbe","000000"),
                                name="Zigosity")
        
        res<-p
      }
      res
    }
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
             sessionvalues$currentResults$name,br(),
             sessionvalues$currentResults$metadata$controlgroup,br(),
             sessionvalues$currentResults$metadata$pathogroup,br(),
             sessionvalues$currentResults$metadata$timestart,br(),
             sessionvalues$currentResults$metadata$timeend
      )
      
    )
  })
  
  output$resultsPanel<-renderUI({
    fluidRow(
      column(12,
             uiOutput("resultsMetadata"),
             hr(),
             h3("Ranking results"),
             DT::dataTableOutput('resultsTable'),
             hr(),
             h3("Result details"),
             fluidRow(
               column(5,offset=1,
                      fluidRow(
                        h4("Variant metadata"),
                        DT::dataTableOutput('variantsMetadataTable')
                      )
               ),
               column(5,offset=1,
                      fluidRow(
                        h4("Genotypes Variants/samples"),
                        plotOutput("heatMapGenotypes")
                      )
               )
             )
      )
      
      #DT::dataTableOutput('tt')
    )
    
  })
  
  
})

