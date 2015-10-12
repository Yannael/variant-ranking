library(shiny)
library(DT)
library(plyr)
library(ggplot2)
library(queryBuildR)
library(shinyBS)
library(rpivotTable)
library(httr)

shinyServer(function(input, output,session) {
  sessionvalues <- reactiveValues()
  sessionvalues$data<-loadData("")$data
  
  sessionvalues$samplesSets<-loadSet("samplesSets","")
  sessionvalues$selectedSamplesSet<-"All"
  sessionvalues$noUpdateFollowingSave<-F
  
  sessionvalues$nbRowsExceededWarningMessage<-""
  sessionvalues$analysesNames<-analysesNames
  sessionvalues$checkboxAddMetadata<-F
  sessionvalues$variantDataGene<-NULL
  
  
  ####################################################
  #Sample group manager
  ####################################################
  
  output$nbRowsExceededWarningMessage<-renderText({
    sessionvalues$nbRowsExceededWarningMessage
  })
  
  output$selectSampleGroupUI<-renderUI({
    selectInput('selectedSampleGroup', 'Select filter', 
                choices = list("Groups"=sessionvalues$samplesSets$group), 
                selected=sessionvalues$selectedSamplesSet,
                selectize = FALSE)
  })
  
  #Select group
  observe({
    if (length(input$selectedSampleGroup)>0) {
      withProgress(min=1, max=3, expr={
        setProgress(message = 'Retrieving data, please wait...',
                    value=2)
        
      #if (!sessionvalues$noUpdateFollowingSave) {
      selectedSampleGroupIndex<-which(sessionvalues$samplesSets$group==input$selectedSampleGroup)
      if (length(selectedSampleGroupIndex)==0) selectedSampleGroupIndex<-1
      sql<-sessionvalues$samplesSets[selectedSampleGroupIndex,2]
      sqlQuery<-sql
      if (sqlQuery=="") sqlQuery<-"reset"
      session$sendCustomMessage(type='callbackHandlerSelectSampleGroup', sqlQuery)
      data<-loadData(sql)
      sessionvalues$data<-data$data
      sessionvalues$nbRowsExceededWarningMessage<-data$nbRowsExceededWarningMessage
      #}
      #else {
      #  isolate({sessionvalues$noUpdateFollowingSave<-F})
      #}
      
      })
    }
  })
  
  #Delete group
  observe({
    if (length(input$deleteConfirmYesButtonSampleGroup)>0 & input$deleteConfirmYesButtonSampleGroup) {
      isolate({
        sessionvalues$samplesSets<-sessionvalues$samplesSets[-which(sessionvalues$samplesSets$group==input$selectedSampleGroup),]
        groupsdb<-dbConnect(RSQLite::SQLite(), "../samplesSets.db")
        dbWriteTable(groupsdb,"samplesSets",sessionvalues$samplesSets,overwrite=T,row.names=F)
        dbDisconnect(groupsdb)
      })
      toggleModal(session, "deleteConfirmSampleGroup", toggle = "close")
      sessionvalues$samplesSets<-loadSet("samplesSets","")
      #session$sendCustomMessage(type='callbackHandlerSelectSampleGroup', "")
    }
  })
  
  #Apply filters
  observe({
    if (length(input$sqlQuerySamplesValue)) {
      withProgress(min=1, max=3, expr={
        setProgress(message = 'Retrieving data, please wait...',
                    value=2)
        data<-loadData(input$sqlQuerySamplesValue)
        sessionvalues$data<-data$data
        sessionvalues$nbRowsExceededWarningMessage<-data$nbRowsExceededWarningMessage
      })
    }
  })
  
  #Save
  observe({
    input$samplesQuerySave2
    isolate({
      sessionvalues$noUpdateFollowingSave<-T
      if ((length(input$sampleGoupNameSave)>0) & length(input$sqlQuerySamplesValue)>0) {
        groupsdb<-dbConnect(RSQLite::SQLite(), "../samplesSets.db")
        data<-data.frame(input$sampleGoupNameSave,input$sqlQuerySamplesValue)
        colnames(data)<-c("group","sql")
        dbWriteTable(groupsdb,"samplesSets",data,append=T)
        dbDisconnect(groupsdb)
        toggleModal(session, "modalSamplesQuerySave", toggle = "close")
        sessionvalues$samplesSets<-loadSet("samplesSets","")
        sessionvalues$selectedSamplesSet<-input$sampleGoupNameSave
      }
    })
  })
  
  output$queryBuilderSamples<-renderQueryBuildR({
    filters<-filtersTypes
    rules<-NULL
    queryBuildR(rules,filters)
  })
  
  output$showVarPhenotypeUI<-renderUI({
    isolate({
      niceNames<-as.vector(sapply(colnames(sessionvalues$data),idToName))
      selectInput('showVarPhenotype', 'Select variables to display', niceNames, 
                  selected=niceNames[c(1,36:37,39,2:7)],multiple=TRUE, selectize=TRUE,width='1050px')
    })
  })
  
  output$phenotypesTable<-DT::renderDataTable({
    if (length(input$showVarPhenotype)>0) {
      data<-sessionvalues$data[,sapply(input$showVarPhenotype,nameToId)]
      data[is.na(data)]<-''
      colnames(data)<-input$showVarPhenotype
      getWidgetTable(data,session)
    }
  },server=T)
  
  output$downloadSelection <- downloadHandler(
    filename = function() {
      paste('variantSelection.zip', sep='')
    },
    content = function(con) {
      write.csv(sessionvalues$data, file="variantSelection.csv", row.names=F)
      zip(con,c('variantSelection.zip'))
    }
  )
  output$pivotTable<-renderRpivotTable({
    if (length(input$showVarPhenotype)>0) {
      data<-sessionvalues$data[,sapply(input$showVarPhenotype,nameToId)]
      colnames(data)<-input$showVarPhenotype
      rpivotTable(data,width='1050px')
    }
  })
  
  ####################################################
  #Ranking engine
  ####################################################
  
  output$selectSampleGroup1UI<-renderUI({
    selectInput('selectSampleGroup1', 'Control group', 
                choices = list("Groups"=sessionvalues$samplesSets$group), 
                selected=sessionvalues$samplesSets$group[1],
                selectize = FALSE)
  })
  
  output$selectSampleGroup2UI<-renderUI({
    selectInput('selectSampleGroup2', 'Case group', 
                choices = list("Groups"=sessionvalues$samplesSets$group), 
                selected=sessionvalues$samplesSets$group[1],
                selectize = FALSE)
  })
  
  observe({
    input$startAnalysisButton
    isolate({
      analysisName<-input$analysisName
      scope<-input$rankingScope
      scale<-input$rankingScale
    })
    
    if (analysisName!='') {
      analysis<-list()
      
      isolate({
        sampleGroup1name<-input$selectSampleGroup1
        sampleGroup2name<-input$selectSampleGroup2
        selectSampleGroupIndex1<-which(sessionvalues$samplesSets$group==input$selectSampleGroup1)
        selectSampleGroupIndex2<-which(sessionvalues$samplesSets$group==input$selectSampleGroup2)
      })
      
      analysis$group1name<-sampleGroup1name
      analysis$group2name<-sampleGroup2name
      
      analysis$group1<-sessionvalues$samplesSets[selectSampleGroupIndex1,]
      analysis$group2<-sessionvalues$samplesSets[selectSampleGroupIndex2,]
      
      controlData<-loadData(sessionvalues$samplesSets[selectSampleGroupIndex1,'sql'],noLimit=T,excludeID=F)[[1]]
      caseData<-loadData(sessionvalues$samplesSets[selectSampleGroupIndex2,'sql'],noLimit=T,excludeID=F)[[1]]
      
      samplesID<-c(unique(caseData[,'Sample_ID']),unique(controlData[,'Sample_ID']))
      nCases<-length(unique(caseData[,'Sample_ID']))
      
      variantData<-rbind(controlData,caseData)
      variantData[is.na(variantData)]<-''
      analysis$variantData<-variantData
      analysis$samplesID<-samplesID
      analysis$nCases<-nCases
      
      save(file=paste("analyses/",analysisName,".Rdata",sep=""),analysis)
      createGenotypeMatrix(analysisName)
      
      startCommand<-paste("spark-submit --master local --conf spark.eventLog.enabled=true GVR.py",analysisName,nCases,scope,scale)
      system(startCommand)
      
    }  
  })
  
  
  ####################################################
  #Results explorer
  ####################################################
  
  observe({
    input$refreshResultsButton
    analysesFiles<-dir("analyses/","*.Rdata")
    sessionvalues$analysesNames<-as.vector(unlist(sapply(analysesFiles,strsplit,'.Rdata')))
  })
  
  output$selectAnalysisUI<-renderUI({
    selectInput('selectAnalysis', 'Select analysis', choices = list(
      "Available analyses" = sessionvalues$analysesNames
    ), selected=sessionvalues$analysesNames[1],selectize = FALSE)
  })
  
  dummy<-function() {
    r <- POST("http://api.omim.org/api/apiKey",query = list(apiKey = "DE2F57B0E7899D4D81351AFCFE44914C610ABF03"))
    res <- GET("http://api.omim.org/api/gene",query = list(search = "mental"))  
  }
  
  output$showVarResultsUI<-renderUI({
    nameAnalysis<-input$selectAnalysis
    
    if (sessionvalues$results$scale=="variant") {
      if (sessionvalues$results$scope=="monogenic") {
        niceNames<-as.vector(sapply(colnames(sessionvalues$results$scoreSummary),idToName))[-4]
        initialSelect<-niceNames[c(1:5,9)]
      }
    }
    
    if (sessionvalues$results$scale=="gene") {
      if (sessionvalues$results$scope=="monogenic") {
        niceNames<-as.vector(sapply(colnames(sessionvalues$results$scoreSummary),idToName))
        initialSelect<-niceNames[1:4]
      }
      if (sessionvalues$results$scope=="digenic") {
        niceNames<-as.vector(sapply(colnames(sessionvalues$results$scoreSummary),idToName))
        initialSelect<-niceNames[c(1:4,8)]
      }
    }
    selectInput('showVarResults', 'Select variables to display', niceNames, 
                selected=initialSelect,multiple=TRUE, selectize=TRUE,width='1050px')
  })
  
  
  output$resultsTable<-DT::renderDataTable({
    if (length(input$showVarResults)>0) {
      nameAnalysis<-input$selectAnalysis
      if (length(setdiff(sapply(input$showVarResults,nameToId),colnames(sessionvalues$results$scoreSummary)))==0) {
        data<-sessionvalues$results$scoreSummary[,sapply(input$showVarResults,nameToId)]
        targetsShort<-which(colnames(data)!="Gene_Symbol")
        colnames(data)<-input$showVarResults
        getWidgetTable(data,session,selection='single',targetsShort=targetsShort)
      }
    }
  },server=TRUE)
  
  observe({
    nameAnalysis<-input$selectAnalysis
    if (length(nameAnalysis)>0) {
      load(paste("analyses/",nameAnalysis,".Rdata",sep=""))
      sessionvalues$analysis<-analysis
      sessionvalues$results<-procRes(fromJSON(txt=paste0("analyses/",nameAnalysis,".txt")),sessionvalues$analysis)
      if (length(input$resultsTable_rows_selected)) {
        
        if (sessionvalues$results$scale=="variant") {
          if (sessionvalues$results$scope=="monogenic") {
            variantsLocus<-sessionvalues$results$scoreSummary[input$resultsTable_rows_selected,'ID']
            i.match<-which(sessionvalues$analysis$variantData[,'ID']==variantsLocus)
            sessionvalues$variantDataGene<-sessionvalues$analysis$variantData[i.match,-c(1)]
          }
        }
        
        if (sessionvalues$results$scale=="gene") {
          if (sessionvalues$results$scope=="monogenic") {
            geneID<-sessionvalues$results$scoreSummary[input$resultsTable_rows_selected,'Gene_Ensembl']
            data<-sessionvalues$analysis$variantData[,-c(1)]
            sessionvalues$variantDataGene<-data[which(data[,'Gene_Ensembl']==geneID),]
          }
          if (sessionvalues$results$scope=="digenic") {
            geneID1<-sessionvalues$results$scoreSummary[input$resultsTable_rows_selected,'Gene_Ensembl1']
            geneID2<-sessionvalues$results$scoreSummary[input$resultsTable_rows_selected,'Gene_Ensembl2']
            data<-sessionvalues$analysis$variantData[,-c(1)]
            sessionvalues$variantDataGene<-data[which(data[,'Gene_Ensembl']==geneID1 | data[,'Gene_Ensembl']==geneID2),]
          }
        }
        
      }
    }
  })
  
  observe({
    data<-NULL
    nameAnalysis<-input$selectAnalysis
    if (!is.null(input$checkboxAddMetadata)) {
      
      if (input$checkboxAddMetadata & !sessionvalues$checkboxAddMetadata){
        sessionvalues$checkboxAddMetadata<-T
        isolate({
          samplesID<-sessionvalues$analysis$samplesID
          if (sessionvalues$results$scale=="variant") {
            if (sessionvalues$results$scope=="monogenic") {
              variantsLocus<-sessionvalues$results$scoreSummary[input$resultsTable_rows_selected,'ID']
              data<-loadData(paste0("ID='",variantsLocus,"'"),noLimit=T)[[1]]
              data<-data[which(data[,'Sample_ID'] %in% samplesID),]
              data[is.na(data)]<-''
            }
          }
          
          if (sessionvalues$results$scale=="gene") {
            if (sessionvalues$results$scope=="monogenic") {
              geneID<-sessionvalues$results$scoreSummary[input$resultsTable_rows_selected,'Gene_Ensembl']
              data<-loadData(paste0("Gene_Ensembl='",geneID,"'"),noLimit=T)[[1]]
              data<-data[which(data[,'Sample_ID'] %in% samplesID),]
              data[is.na(data)]<-''
            }
            if (sessionvalues$results$scope=="digenic") {
              geneID1<-sessionvalues$results$scoreSummary[input$resultsTable_rows_selected,'Gene_Ensembl1']
              geneID2<-sessionvalues$results$scoreSummary[input$resultsTable_rows_selected,'Gene_Ensembl2']
              data1<-loadData(paste0("Gene_Ensembl='",geneID1,"'"),noLimit=T)[[1]]
              data2<-loadData(paste0("Gene_Ensembl='",geneID2,"'"),noLimit=T)[[1]]
              data<-rbind(data1,data2)
              data<-data[which(data[,'Sample_ID'] %in% samplesID),]
              data[is.na(data)]<-''
            }
          }
        })
        sessionvalues$variantDataGene<-data
        
        
      }
      if (!input$checkboxAddMetadata & sessionvalues$checkboxAddMetadata){
        sessionvalues$checkboxAddMetadata<-F
        isolate({
          if (sessionvalues$results$scale=="variant") {
            if (sessionvalues$results$scope=="monogenic") {
              variantsLocus<-sessionvalues$results$scoreSummary[input$resultsTable_rows_selected,'ID']
              i.match<-which(sessionvalues$analysis$variantData[,'ID']==variantsLocus)
              data<-sessionvalues$analysis$variantData[i.match,-c(1)]
            }
          }
          
          if (sessionvalues$results$scale=="gene") {
            if (sessionvalues$results$scope=="monogenic") {
              geneID<-sessionvalues$results$scoreSummary[input$resultsTable_rows_selected,'Gene_Symbol']
              i.match<-which(sessionvalues$analysis$variantData[,'Gene_Symbol']==geneID)
              data<-sessionvalues$analysis$variantData[i.match,-c(1)]
            }
            if (sessionvalues$results$scope=="digenic") {
              geneID1<-sessionvalues$results$scoreSummary[input$resultsTable_rows_selected,'Gene_Symbol1']
              geneID2<-sessionvalues$results$scoreSummary[input$resultsTable_rows_selected,'Gene_Symbol2']
              i.match<-which(sessionvalues$analysis$variantData[,'Gene_Symbol']==geneID1 | sessionvalues$analysis$variantData[,'Gene_Symbol']==geneID2)
              data<-sessionvalues$analysis$variantData[i.match,-c(1)]
            }
          }
        })
        sessionvalues$variantDataGene<-data
        
      }
    }
  })
  
  output$showVarMetadataUI<-renderUI({
    nameAnalysis<-input$selectAnalysis
    if (length(sessionvalues$variantDataGene)>0) {
      niceNames<-as.vector(sapply(colnames(sessionvalues$variantDataGene),idToName))
      selectInput('showVarMetadata', 'Select variables to display', niceNames, 
                  selected=niceNames[c(1:7,23:25,16,37)],multiple=TRUE, selectize=TRUE,width='1050px')
    }
  })
  
  output$variantsMetadataTable<-DT::renderDataTable({
    nameAnalysis<-input$selectAnalysis
    if (length(sessionvalues$variantDataGene)>0 & length(input$showVarMetadata)>0) {
      data<-sessionvalues$variantDataGene[,sapply(input$showVarMetadata,nameToId)]
      getWidgetTable(data,session)
    }
  },server=TRUE)
  
  output$variantsMetadataPivotTable<-renderRpivotTable({
    input$checkboxAddMetadata
    nameAnalysis<-input$selectAnalysis
    if (length(sessionvalues$variantDataGene)>0 & length(input$showVarMetadata)>0) {
      data<-sessionvalues$variantDataGene[,sapply(input$showVarMetadata,nameToId)]
      rpivotTable(data)
    }
  })
  
  
  output$resultsMetadata<-renderUI({
    fluidRow(
      column(3,
             strong("Control group: "),br(),
             strong("Pathological group: "),br(),
             strong("Start time: "),br(),
             strong("End time: "),br(),
             strong("Total run time: ")
      ),
      column(4,
             sessionvalues$analysis$group1name,br(),
             sessionvalues$analysis$group2name,br(),
             sessionvalues$results$start_time,br(),
             sessionvalues$results$end_time,br(),
             paste(sessionvalues$results$run_time, "seconds")
      )
      
    )
  })
  
  output$resultsPanel<-renderUI({
    fluidRow(
      column(12,
             uiOutput("resultsMetadata"),
             hr(),
             h3("Ranking results"),
             uiOutput("showVarResultsUI"),
             DT::dataTableOutput('resultsTable'),
             hr(),
             h3("Advanced browsing"),
             fluidRow(
               column(12,
                      checkboxInput("checkboxAddMetadata", label = "Include all variants matching sample IDs ", value = FALSE),
                      fluidRow(
                        uiOutput("showVarMetadataUI"),
                        dataTableOutput('variantsMetadataTable')
                      ),
                      fluidRow(
                        h4("Pivot Table"),
                        rpivotTableOutput('variantsMetadataPivotTable')
                      )
               )
             )
      )
      
      #DT::dataTableOutput('tt')
    )
    
  })
  
  
})

