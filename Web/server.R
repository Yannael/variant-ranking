library(shiny)
library(DT)
library(plyr)
library(ggplot2)
library(queryBuildR)
library(shinyBS)

shinyServer(function(input, output,session) {
  sessionvalues <- reactiveValues()
  sessionvalues$phenotypes<-loadData("phenotypes","")
  sessionvalues$variants<-loadData("variants","")
  sessionvalues$samplesSets<-loadData("samplesSets","")
  sessionvalues$variantsSets<-loadData("variantsSets","")
  
  ####################################################
  #Sample group manager
  ####################################################
  
  output$selectSampleGroupUI<-renderUI({
    selectInput('selectedSampleGroup', 'Select sample set', 
                choices = list("Groups"=sessionvalues$samplesSets$group), 
                selected=sessionvalues$samplesSets$group[1],
                selectize = FALSE)
  })
  
  observe({
    if (length(input$selectedSampleGroup)>0) {
      selectedSampleGroupIndex<-which(sessionvalues$samplesSets$group==input$selectedSampleGroup)
      if (length(selectedSampleGroupIndex)==0) selectedSampleGroupIndex<-1
      sql<-sessionvalues$samplesSets[selectedSampleGroupIndex,2]
      session$sendCustomMessage(type='callbackHandlerSelectSampleGroup', sql)
      sessionvalues$phenotypes<-loadData("phenotypes",sql)
    }
  })
  
  observe({
    if (length(input$deleteConfirmYesButtonSampleGroup)>0 & input$deleteConfirmYesButtonSampleGroup) {
      isolate({
        sessionvalues$samplesSets<-sessionvalues$samplesSets[-which(sessionvalues$samplesSets$group==input$selectedSampleGroup),]
        groupsdb<-dbConnect(RSQLite::SQLite(), "../samplesSets.db")
        dbWriteTable(groupsdb,"samplesSets",sessionvalues$samplesSets,overwrite=T,row.names=F)
        dbDisconnect(groupsdb)
        toggleModal(session, "deleteConfirmSampleGroup", toggle = "close")
      })
      sessionvalues$samplesSets<-loadData("samplesSets","")
    }
  })
  
  observe({
    if (length(input$sqlQuerySamplesValue))
      sessionvalues$phenotypes<-loadData("phenotypes",input$sqlQuerySamplesValue)
  })
  
  observe({
    input$samplesQuerySave2
    isolate({
      if ((length(input$sampleGoupNameSave)>0) & length(input$sqlQuerySamplesValue)>0) {
        groupsdb<-dbConnect(RSQLite::SQLite(), "../samplesSets.db")
        data<-data.frame(input$sampleGoupNameSave,input$sqlQuerySamplesValue)
        colnames(data)<-c("group","sql")
        dbWriteTable(groupsdb,"samplesSets",data,append=T)
        dbDisconnect(groupsdb)
        toggleModal(session, "modalSamplesQuerySave", toggle = "close")
        sessionvalues$samplesSets<-loadData("samplesSets","")
      }
    })
  })
  
  output$queryBuilderSamples<-renderQueryBuildR({
    data<-sessionvalues$phenotypes
    data[,1]<-as.factor(data[,1])
    data[,3]<-as.factor(data[,3])
    data[,4]<-as.factor(data[,4])
    data[,5]<-as.factor(data[,5])
    data[,6]<-as.factor(data[,6])
    filters<-getFiltersFromTable(data)
    rules=NULL
    queryBuildR(rules,filters)
  })
  
  output$showVarPhenotypeUI<-renderUI({
    isolate({
      niceNames<-as.vector(sapply(colnames(sessionvalues$phenotypes),idToName))
      selectInput('showVarPhenotype', 'Display variables (max. 6)', niceNames, 
                  selected=niceNames[1:5],multiple=TRUE, selectize=TRUE)
    })
  })
  
  output$phenotypesTable<-DT::renderDataTable({
    if (length(input$showVarPhenotype)>0) {
      data<-sessionvalues$phenotypes[,sapply(input$showVarPhenotype,nameToId)]
      colnames(data)<-input$showVarPhenotype
      getWidgetTable(data,session)
    }
  })
  
  ####################################################
  #Variant filtering
  ####################################################
  
  output$selectVariantGroupUI<-renderUI({
    #browser()
    selectInput('selectedVariantGroup', 'Select variant set', 
                choices = list("Groups"=c(sessionvalues$variantsSets$group)), 
                selected=sessionvalues$variantsSets$group[1],
                selectize = FALSE)
  })
  
  observe({
    if (length(input$selectedVariantGroup)>0) {
      selectedVariantGroupIndex<-which(sessionvalues$variantsSets$group==input$selectedVariantGroup)
      if (length(selectedVariantGroupIndex)==0) selectedVariantGroupIndex<-1
      sql<-sessionvalues$variantsSets[selectedVariantGroupIndex,2]
      session$sendCustomMessage(type='callbackHandlerSelectVariantGroup', sql)
      sessionvalues$variants<-loadData("variants",sql)
    }
  })
  
  observe({
    if (length(input$sqlQueryVariantsValue))
      sessionvalues$variants<-loadData("variants",input$sqlQueryVariantsValue)
  })
  
  observe({
    input$variantsQuerySave2
    isolate({
      if ((length(input$variantGoupNameSave)>0) & length(input$sqlQueryVariantsValue)>0) {
        groupsdb<-dbConnect(RSQLite::SQLite(), "../variantsSets.db")
        data<-data.frame(input$variantGoupNameSave,input$sqlQueryVariantsValue)
        colnames(data)<-c("group","sql")
        dbWriteTable(groupsdb,"variantsSets",data,append=T)
        dbDisconnect(groupsdb)
        sessionvalues$variantsSets<-loadData("variantsSets","")
        toggleModal(session, "modalVariantsQuerySave", toggle = "close")
      }
    })
  })
  
  output$queryBuilderVariants<-renderQueryBuildR({
    data<-as.data.frame(sessionvalues$variants)[,-1]
    rules<-NULL
    queryBuildR(rules,filtersVariantsTypes)
  })
  
  output$showVarVariantsUI<-renderUI({
    isolate({
      niceNames<-as.vector(sapply(colnames(sessionvalues$variants),idToName))
    selectizeInput('showVarVariants', 'Display variables (max. 6)', niceNames[2:ncol(sessionvalues$variants)], 
                selected=niceNames[2:6],multiple=TRUE, options = list(maxItems = 6))
    })
  })

  output$variantsTable<-DT::renderDataTable({
    if (length(input$showVarVariants)>0) {
    data<-as.data.frame(sessionvalues$variants)[,sapply(input$showVarVariants,nameToId)]
    colnames(data)<-input$showVarVariants
    getWidgetTable(data,session)
    }
  })
  
  ####################################################
  #Results explorer
  ####################################################
  
  
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
      widget<-NULL
      if (sessionvalues$currentResults$type=="singleVariant") {
        infovariants<-sessionvalues$currentResults$infovariants[input$resultsTable_rows_selected,]
        data<-cbind(names(infovariants),t(infovariants))
        colnames(data)<-c("Variable","Value")
      }
      if (sessionvalues$currentResults$type=="pairVariantsMonogenic") {
        infovariants1<-sessionvalues$currentResults$infovariants1[input$resultsTable_rows_selected,]
        infovariants2<-sessionvalues$currentResults$infovariants2[input$resultsTable_rows_selected,]
        data<-cbind(names(infovariants1),t(infovariants1),t(infovariants2))
        colnames(data)<-c("Variable","Value1","Value2")
      }
      data[is.na(data)]<-"NA"
      action = dataTableAjax(session, data,rownames=F)
      if (!is.null(dim(data))) {
        widget <- datatable(data, 
                            server = TRUE, 
                            selection = 'single',
                            rownames=F,
                            filter = 'bottom',
                            escape=T,
                            options = list(
                              dom= 't',
                              pageLength = 100,
                              scrollY = 335,
                              ajax = list(url = action),
                              columnDefs = list(
                                list(
                                  targets = c(1),
                                  render = JS(
                                    "function(data, type, row, meta) {",
                                    "return type === 'display' && data.length > 15 ?",
                                    "'<span title=\"' + data + '\">' + data.substr(0, 11) + '...</span>' : data;",
                                    "}")
                                ),
                                list(className="dt-right",targets="_all")
                              )
                            )
        )
      }
      widget 
    }
  })
  
  getMatGenoTrio<-function(genotypes,namecols) {
    genotypes<-factor(genotypes,levels=c("0","1","2","NA"))
    genotypes<-mapvalues(genotypes,from=c(0,1,2,"NA"),to=c("Ref/Ref","Ref/Alt","Alt/Alt","NA"))
    
    DF<-data.frame(rep(paste(rep("Trio ",11),1:11,sep=""),each=3),rep(namecols,11),genotypes)
    DF[,1]<-factor(DF[,1],levels=rev(unique(DF[,1])),ordered=T)
    colnames(DF)<-c("Trios","FMC","Zygosity")
    
    DF
  }
  
  output$heatMapGenotypes<-renderPlot({
    if (length(input$resultsTable_rows_selected)>0) {
      res<-NULL
      DF<-NULL
      genotypes<-sessionvalues$currentResults$genotypes[[input$resultsTable_rows_selected]]
      if (sessionvalues$currentResults$type=="singleVariant") {
        DF<-getMatGenoTrio(genotypes,c("Child","Father","Mother"))
        aspect.ratio<-3
      }
      if (sessionvalues$currentResults$type=="pairVariantsMonogenic") {
        DF1<-getMatGenoTrio(genotypes[1,],c("Child variant 1","Father variant 1","Mother variant 1"))
        DF2<-getMatGenoTrio(genotypes[2,],c("Child variant 2","Father variant 2","Mother variant 2"))
        DF<-rbind(DF1,DF2)
        aspect.ratio<-1.7
      }
      
      #DF[,1]<-factor(DF[,1],levels=rev(unique(DF[,1])),ordered=T)
      DF[,2]<-factor(DF[,2],levels=sort(levels(DF[,2])),ordered=T)
      
      base_size=20
      p<-ggplot(DF, aes(y=Trios,x=FMC))
      p <- p + geom_tile(aes(fill=Zygosity),colour = "white")
      p <- p +  theme_grey(base_size = base_size) + labs(x = "", y = "") + 
        scale_x_discrete(expand = c(0, 0)) +
        scale_y_discrete(expand = c(0, 0)) +
        theme(aspect.ratio=aspect.ratio)+
        theme(axis.title.x = element_blank(),
              axis.text.x = element_text(size = base_size *  0.8, angle = 330,hjust = 0, colour = "grey50")) 
      p<-p+ scale_fill_manual(breaks=c("Ref/Ref","Ref/Alt","Alt/Alt","NA"),
                              values=c("#fff7bc","#fec44f","#d95f0e","000000"),
                              #values=c("#ece7f2","#a6bddb","#2b8cbe","000000"),
                              name="Zigosity")
      p
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
                        h4("Genotypes trios"),
                        plotOutput("heatMapGenotypes")
                      )
               )
             )
      )
      
      #DT::dataTableOutput('tt')
    )
    
  })
  
  
})

