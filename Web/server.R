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
  sessionvalues$variants<-loadData("")$data
  sessionvalues$phenotypes<-loadPhenotypes("")
  
  sessionvalues$variantsSets<-loadSet("variantsSets","")
  sessionvalues$phenotypesSets<-loadSet("phenotypesSets","")
  sessionvalues$selectedVariantsSet<-"All"
  sessionvalues$selectedPhenotypesSet<-"All"
  sessionvalues$noUpdateFollowingSave<-F
  
  sessionvalues$nbRowsExceededWarningMessage<-""
  sessionvalues$analysesNames<-analysesNames
  sessionvalues$checkboxAddMetadata<-F
  sessionvalues$variantDataGene<-NULL
  
  
  ####################################################
  #Phenotypes group manager
  ####################################################
  
  output$selectPhenotypesGroupUI<-renderUI({
    selectInput('selectedPhenotypesGroup', 'Select sample group', 
                choices = list("Groups"=sessionvalues$phenotypesSets$group), 
                selected=sessionvalues$selectedPhenotypesSet,
                selectize = FALSE)
  })
  
  #Select group
  observe({
    if (length(input$selectedPhenotypesGroup)>0) {
      selectedPhenotypesGroupIndex<-which(sessionvalues$phenotypesSets$group==input$selectedPhenotypesGroup)
      if (length(selectedPhenotypesGroupIndex)==0) selectedPhenotypesGroupIndex<-1
      sql<-sessionvalues$phenotypesSets[selectedPhenotypesGroupIndex,2]
      sqlQuery<-sql
      if (sqlQuery=="") sqlQuery<-"reset"
      session$sendCustomMessage(type='callbackHandlerSelectPhenotypesGroup', sqlQuery)
      sessionvalues$phenotypes<-loadPhenotypes(sql)
    }
  })
  
  #Delete group
  observe({
    if (length(input$deleteConfirmYesButtonPhenotypesGroup)>0 & input$deleteConfirmYesButtonPhenotypesGroup) {
      isolate({
        sessionvalues$phenotypesSets<-sessionvalues$phenotypesSets[-which(sessionvalues$phenotypesSets$group==input$selectedPhenotypesGroup),]
        groupsdb<-dbConnect(RSQLite::SQLite(), "../phenotypesSets.db")
        dbWriteTable(groupsdb,"phenotypesSets",sessionvalues$phenotypesSets,overwrite=T,row.names=F)
        dbDisconnect(groupsdb)
      })
      toggleModal(session, "deleteConfirmPhenotypesGroup", toggle = "close")
      sessionvalues$phenotypes<-loadSet("phenotypesSets","")
      #session$sendCustomMessage(type='callbackHandlerSelectSampleGroup', "")
    }
  })
  
  #Apply filters
  observe({
    if (length(input$sqlQueryPhenotypesValue)) {
      sessionvalues$phenotypes<-loadPhenotypes(input$sqlQueryPhenotypesValue)
    }
  })
  
  #Save
  observe({
    input$phenotypesQuerySave2
    isolate({
      sessionvalues$noUpdateFollowingSave<-T
      if ((length(input$phenotypesGoupNameSave)>0) & length(input$sqlQueryPhenotypesValue)>0) {
        groupsdb<-dbConnect(RSQLite::SQLite(), "../phenotypesSets.db")
        data<-data.frame(input$phenotypesGoupNameSave,input$sqlQueryPhenotypesValue)
        colnames(data)<-c("group","sql")
        dbWriteTable(groupsdb,"phenotypesSets",data,append=T)
        dbDisconnect(groupsdb)
        toggleModal(session, "modalPhenotypesQuerySave", toggle = "close")
        sessionvalues$phenotypesSets<-loadSet("phenotypesSets","")
        sessionvalues$selectedPhenotypesSet<-input$phenotypesGoupNameSave
      }
    })
  })
  
  #Get IDs samples
  output$listSamplesIDs<-renderText({
    listIDs<-sessionvalues$phenotypes$Sample_ID
    result<-paste(listIDs,sep="",collapse=" , ")
    result
  })
  
  output$queryBuilderPhenotypes<-renderQueryBuildR({
    filters<-filtersPhenotypesTypes
    rules<-NULL
    queryBuildR(rules,filters)
  })
  
  output$showVarPhenotypesUI<-renderUI({
    isolate({
      niceNames<-as.vector(sapply(colnames(sessionvalues$phenotypes),idToName))
      selectInput('showVarPhenotypes', 'Select variables to display', niceNames, 
                  selected=niceNames,multiple=TRUE, selectize=TRUE,width='1050px')
    })
  })
  
  output$phenotypesTable<-DT::renderDataTable({
    if (length(input$showVarPhenotypes)>0) {
      data<-sessionvalues$phenotypes[,sapply(input$showVarPhenotypes,nameToId)]
      data[is.na(data)]<-''
      colnames(data)<-input$showVarPhenotypes
      getWidgetTable(data,session)
    }
  },server=T)
  
  output$downloadPhenotypesSelection <- downloadHandler(
    filename = function() {
      paste('phenotypes.zip', sep='')
    },
    content = function(con) {
      write.csv(sessionvalues$phenotypes, file="phenotypes.csv", row.names=F,quote=T)
      zip(con,c('phenotypes.csv'))
    }
  )
  
  output$pivotTablePhenotypes<-renderRpivotTable({
    if (length(input$showVarPhenotypes)>0) {
      data<-sessionvalues$phenotypes[,sapply(input$showVarPhenotypes,nameToId)]
      colnames(data)<-input$showVarPhenotypes
      rpivotTable(data,width='1050px')
    }
  })
  
  ####################################################
  #Sample group manager
  ####################################################
  
  output$nbRowsExceededWarningMessage<-renderText({
    sessionvalues$nbRowsExceededWarningMessage
  })
  
  output$selectVariantsGroupUI<-renderUI({
    selectInput('selectedVariantsGroup', 'Select filter', 
                choices = list("Groups"=sessionvalues$variantsSets$group), 
                selected=sessionvalues$selectedVariantsSet,
                selectize = FALSE)
  })
  
  #Select group
  observe({
    if (length(input$selectedVariantsGroup)>0) {
      withProgress(min=1, max=3, expr={
        setProgress(message = 'Retrieving data, please wait...',
                    value=2)
        
        #if (!sessionvalues$noUpdateFollowingSave) {
        selectedVariantsGroupIndex<-which(sessionvalues$variantsSets$group==input$selectedVariantsGroup)
        if (length(selectedVariantsGroupIndex)==0) selectedVariantsGroupIndex<-1
        sql<-sessionvalues$variantsSets[selectedVariantsGroupIndex,2]
        sqlQuery<-sql
        if (sqlQuery=="") sqlQuery<-"reset"
        session$sendCustomMessage(type='callbackHandlerSelectVariantsGroup', sqlQuery)
        data<-loadData(sql)
        sessionvalues$variants<-data$data
        sessionvalues$nbRowsExceededWarningMessage<-data$nbRowsExceededWarningMessage
      })
    }
  })
  
  #Delete group
  observe({
    if (length(input$deleteConfirmYesButtonVariantsGroup)>0 & input$deleteConfirmYesButtonVariantsGroup) {
      isolate({
        sessionvalues$variantsSets<-sessionvalues$variantsSets[-which(sessionvalues$variantsSets$group==input$selectedVariantsGroup),]
        groupsdb<-dbConnect(RSQLite::SQLite(), "../variantsSets.db")
        dbWriteTable(groupsdb,"variantsSets",sessionvalues$variantsSets,overwrite=T,row.names=F)
        dbDisconnect(groupsdb)
      })
      toggleModal(session, "deleteConfirmVariantsGroup", toggle = "close")
      sessionvalues$variantsSets<-loadSet("variantsSets","")
    }
  })
  
  #Apply filters
  observe({
    if (length(input$sqlQueryVariantsValue)) {
      withProgress(min=1, max=3, expr={
        setProgress(message = 'Retrieving data, please wait...',
                    value=2)
        data<-loadData(input$sqlQueryVariantsValue)
        sessionvalues$variants<-data$data
        sessionvalues$nbRowsExceededWarningMessage<-data$nbRowsExceededWarningMessage
      })
    }
  })
  
  #Save
  observe({
    input$variantsQuerySave2
    isolate({
      sessionvalues$noUpdateFollowingSave<-T
      if ((length(input$variantsGoupNameSave)>0) & length(input$sqlQueryVariantsValue)>0) {
        groupsdb<-dbConnect(RSQLite::SQLite(), "../variantsSets.db")
        data<-data.frame(input$variantsGoupNameSave,input$sqlQueryVariantsValue)
        colnames(data)<-c("group","sql")
        dbWriteTable(groupsdb,"variantsSets",data,append=T)
        dbDisconnect(groupsdb)
        toggleModal(session, "modalVariantsQuerySave", toggle = "close")
        sessionvalues$variantsSets<-loadSet("variantsSets","")
        sessionvalues$selectedvariantsSet<-input$variantsGoupNameSave
      }
    })
  })
  
  output$queryBuilderVariants<-renderQueryBuildR({
    filters<-filtersTypes
    rules<-NULL
    queryBuildR(rules,filters)
  })
  
  output$showVarVariantsUI<-renderUI({
    isolate({
      niceNames<-as.vector(sapply(colnames(sessionvalues$variants),idToName))
      selectInput('showVarVariants', 'Select variables to display', niceNames, 
                  #              selected=niceNames[c(1,36:37,39,2:7)],multiple=TRUE, selectize=TRUE,width='1050px')
                  selected=niceNames[c(1,2:7)],multiple=TRUE, selectize=TRUE,width='1050px')
    })
  })
  
  output$variantsTable<-DT::renderDataTable({
    if (length(input$showVarVariants)>0) {
      data<-sessionvalues$variants[,sapply(input$showVarVariants,nameToId)]
      data[is.na(data)]<-''
      colnames(data)<-input$showVarVariants
      getWidgetTable(data,session)
    }
  },server=T)
  
  output$downloadVariantsSelection <- downloadHandler(
    filename = function() {
      paste('variantSelection.zip', sep='')
    },
    content = function(con) {
      write.csv(sessionvalues$variants, file="variantSelection.csv", row.names=F,quote=T)
      zip(con,c('variantSelection.csv'))
    }
  )
  output$pivotTableVariants<-renderRpivotTable({
    if (length(input$showVarVariants)>0) {
      data<-sessionvalues$variants[,sapply(input$showVarVariants,nameToId)]
      colnames(data)<-input$showVarVariants
      rpivotTable(data,width='1050px')
    }
  })
  
  ####################################################
  #Ranking engine
  ####################################################
  
  output$selectSampleGroup1UI<-renderUI({
    selectInput('selectSampleGroup1', 'Control group', 
                choices = list("Groups"=sessionvalues$variantsSets$group), 
                selected=sessionvalues$variantsSets$group[1],
                selectize = FALSE)
  })
  
  output$selectSampleGroup2UI<-renderUI({
    selectInput('selectSampleGroup2', 'Case group', 
                choices = list("Groups"=sessionvalues$variantsSets$group), 
                selected=sessionvalues$variantsSets$group[1],
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
        selectSampleGroupIndex1<-which(sessionvalues$variantsSets$group==input$selectSampleGroup1)
        selectSampleGroupIndex2<-which(sessionvalues$variantsSets$group==input$selectSampleGroup2)
      })
      
      analysis$group1name<-sampleGroup1name
      analysis$group2name<-sampleGroup2name
      
      analysis$group1<-sessionvalues$variantsSets[selectSampleGroupIndex1,]
      analysis$group2<-sessionvalues$variantsSets[selectSampleGroupIndex2,]
      
      controlData<-loadData(sessionvalues$variantsSets[selectSampleGroupIndex1,'sql'],noLimit=T,excludeID=F)[[1]]
      caseData<-loadData(sessionvalues$variantsSets[selectSampleGroupIndex2,'sql'],noLimit=T,excludeID=F)[[1]]
      
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

