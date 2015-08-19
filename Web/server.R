library(shiny)
library(DT)
library(plyr)
library(ggplot2)
library(queryBuildR)
library(shinyBS)
library(rpivotTable)

shinyServer(function(input, output,session) {
  sessionvalues <- reactiveValues()
  sessionvalues$data<-loadData("")$data
  sessionvalues$samplesSets<-loadSet("samplesSets","")
  sessionvalues$nbRowsExceededWarningMessage<-""
  sessionvalues$results<-resultsAll
  sessionvalues$analyses<-analyses
  sessionvalues$analysesNames<-analysesNames
  
  ####################################################
  #Sample group manager
  ####################################################
  
  output$nbRowsExceededWarningMessage<-renderText({
    sessionvalues$nbRowsExceededWarningMessage
  })
  
  output$selectSampleGroupUI<-renderUI({
    selectInput('selectedSampleGroup', 'Select filter', 
                choices = list("Groups"=sessionvalues$samplesSets$group), 
                selected=sessionvalues$samplesSets$group[1],
                selectize = FALSE)
  })
  
  #Select group
  observe({
    if (length(input$selectedSampleGroup)>0) {
      selectedSampleGroupIndex<-which(sessionvalues$samplesSets$group==input$selectedSampleGroup)
      if (length(selectedSampleGroupIndex)==0) selectedSampleGroupIndex<-1
      sql<-sessionvalues$samplesSets[selectedSampleGroupIndex,2]
      sqlQuery<-sql
      if (sqlQuery=="") sqlQuery<-"reset"
      session$sendCustomMessage(type='callbackHandlerSelectSampleGroup', sqlQuery)
      data<-loadData(sql)
      sessionvalues$data<-data$data
      sessionvalues$nbRowsExceededWarningMessage<-data$nbRowsExceededWarningMessage
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
      session$sendCustomMessage(type='callbackHandlerSelectSampleGroup', "")
    }
  })
  
  #Apply filters
  observe({
    if (length(input$sqlQuerySamplesValue)) {
      data<-loadData(input$sqlQuerySamplesValue)
      sessionvalues$data<-data$data
      sessionvalues$nbRowsExceededWarningMessage<-data$nbRowsExceededWarningMessage
    }
  })
  
  #Save
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
        isolate({sessionvalues$samplesSets<-loadSet("samplesSets","")})
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
  })
  
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
    isolate({analysisName<-input$analysisName})
    if (analysisName!='') {
      #load("analyses.Rdata")
      #browser()
      analyses<-sessionvalues$analyses
      analyses[[analysisName]]<-list()
      
      selectSampleGroupIndex1<-which(sessionvalues$samplesSets$group==input$selectSampleGroup1)
      selectSampleGroupIndex2<-which(sessionvalues$samplesSets$group==input$selectSampleGroup2)
      
      analyses[[analysisName]]$group1<-sessionvalues$samplesSets[selectSampleGroupIndex1,]
      analyses[[analysisName]]$group2<-sessionvalues$samplesSets[selectSampleGroupIndex2,]
      
      controlData<-loadData(sessionvalues$samplesSets[selectSampleGroupIndex1,'sql'],noLimit=T,excludeID=F)[[1]]
      caseData<-loadData(sessionvalues$samplesSets[selectSampleGroupIndex2,'sql'],noLimit=T,excludeID=F)[[1]]
      
      samplesID<-c(unique(caseData[,'Sample_ID']),unique(controlData[,'Sample_ID']))
      nCases<-length(unique(caseData[,'Sample_ID']))
      
      variantData<-rbind(controlData,caseData)
      variantData[is.na(variantData)]<-''
      analyses[[analysisName]]$variantData<-variantData
      analyses[[analysisName]]$samplesID<-samplesID
      analyses[[analysisName]]$nCases<-nCases
      
      save(file="../analyses.Rdata",analyses)
      sessionvalues$analyses<-analyses
      sessionvalues$analysesNames<-names(analyses)
    }  
  })
  
  
  ####################################################
  #Results explorer
  ####################################################
  
  output$selectAnalysisUI<-renderUI({
    
    selectInput('selectAnalysis', 'Select analysis', choices = list(
      "Available analyses" = sessionvalues$analysesNames
    ), selected=sessionvalues$analysesNames[1],selectize = FALSE)
  })
  
  output$showVarResultsUI<-renderUI({
    nameAnalysis<-input$selectAnalysis
    niceNames<-as.vector(sapply(colnames(sessionvalues$results[[nameAnalysis]][[1]]$scoreSummary),idToName))[-2]
    selectInput('showVarResults', 'Select variables to display', niceNames, 
                selected=niceNames[c(1:3,7)],multiple=TRUE, selectize=TRUE,width='1050px')
  })
  
  
  output$resultsTable<-renderDataTable({
    if (length(input$showVarResults)>0) {
      nameAnalysis<-input$selectAnalysis
      data<-sessionvalues$results[[nameAnalysis]][[1]]$scoreSummary[,sapply(input$showVarResults,nameToId)]
      getWidgetTable(data,session,'single')
    }
  },server=TRUE)
  
  observe({
    nameAnalysis<-input$selectAnalysis
    if (length(input$resultsTable_rows_selected)) {
      if (sessionvalues$results[[nameAnalysis]][[1]]$type=="singleVariant") {
        variantsLocus<-sessionvalues$results[[nameAnalysis]][[1]]$scoreSummary[input$resultsTable_rows_selected,'ID']
        i.match<-which(sessionvalues$analyses[[nameAnalysis]]$variantData[,'ID']==variantsLocus)
        sessionvalues$analyses[[nameAnalysis]]$variantDataGene<-sessionvalues$analyses[[nameAnalysis]]$variantData[i.match,-c(1)]
      }
      if (sessionvalues$results[[nameAnalysis]][[1]]$type=="geneUnivariate") {
        geneID<-sessionvalues$results[[nameAnalysis]][[1]]$scoreSummary[input$resultsTable_rows_selected,'Gene']
        data<-sessionvalues$analyses[[nameAnalysis]]$variantData
        sessionvalues$analyses[[nameAnalysis]]$variantDataGene<-data[which(data[,'Gene_Ensembl']==geneID),]
      }
    }
  })
  
  output$showVarMetadataUI<-renderUI({
    nameAnalysis<-input$selectAnalysis
    if (length(sessionvalues$analyses[[nameAnalysis]]$variantDataGene)>0) {
      niceNames<-as.vector(sapply(colnames(sessionvalues$analyses[[nameAnalysis]]$variantDataGene),idToName))
      selectInput('showVarMetadata', 'Select variables to display', niceNames, 
                  selected=niceNames[c(1:7,23:25)],multiple=TRUE, selectize=TRUE,width='1050px')
    }
  })
  
  output$variantsMetadataTable<-DT::renderDataTable({
    isolate({nameAnalysis<-input$selectAnalysis})
    if (length(sessionvalues$analyses[[nameAnalysis]]$variantDataGene)>0 & length(input$showVarMetadata)>0) {
      data<-sessionvalues$analyses[[nameAnalysis]]$variantDataGene[,sapply(input$showVarMetadata,nameToId)]
      getWidgetTable(data,session)
    }
  })
  
  output$variantsMetadataPivotTable<-renderRpivotTable({
    isolate({nameAnalysis<-input$selectAnalysis})
    if (length(sessionvalues$analyses[[nameAnalysis]]$variantDataGene)>0 & length(input$showVarMetadata)>0) {
      data<-sessionvalues$analyses[[nameAnalysis]]$variantDataGene[,sapply(input$showVarMetadata,nameToId)]
      rpivotTable(data)
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

