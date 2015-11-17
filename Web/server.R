library(shiny)
library(DT)
library(plyr)
library(ggplot2)
library(queryBuildR)
library(shinyBS)
library(rpivotTable)
library(httr)

source("filterPhenotypes.R")
source("filterVariants.R")

shinyServer(function(input, output,session) {
  sessionvalues <- reactiveValues()
  sessionvalues$variants<-loadData("")$data
  sessionvalues$phenotypes<-loadPhenotypes("")
  
  sessionvalues$nbRowsExceededWarningMessage<-""
  sessionvalues$analysesNames<-analysesNames
  sessionvalues$checkboxAddMetadata<-F
  sessionvalues$variantDataGene<-NULL
  
  output<-createFilterPhenotype(input,output,session)
  output<-createFilterVariant(input,output,session)
  
  
  ####################################################
  #Phenotypes group manager
  ####################################################
  
  
  output$answerCliniPhenome<-renderText({
    a<-POST(CliniPhenomeAPI,body=list(ontid="Autistic"))
    try(doc <- htmlParse(a),silent=T)
    tableNodes <- getNodeSet(doc[3]$children$html, "//sample")
    
    data<-paste0(sapply(tableNodes,xmlValue),sep='<br>',collapse="")
    data
  })
  
  
  #Get IDs samples
  output$listSamplesIDs<-renderText({
    listIDs<-sessionvalues$phenotypes$Sample_ID
    result<-paste(listIDs,sep="",collapse=" , ")
    result
  })
  
  observe({
    if (length(input$filterPhenotypeQueryBuilderSQL)>0)
      sessionvalues$phenotypes<-loadPhenotypes(input$filterPhenotypeQueryBuilderSQL)
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
  #Variant group manager
  ####################################################
  
  output$nbRowsExceededWarningMessage<-renderText({
    sessionvalues$nbRowsExceededWarningMessage
  })
  
  #Apply filters
  observe({
    if (length(input$filterVariantQueryBuilderSQL)) {
      withProgress(min=1, max=3, expr={
        setProgress(message = 'Retrieving data, please wait...',
                    value=2)
        data<-loadData(input$filterVariantQueryBuilderSQL)
        sessionvalues$variants<-data$data
        sessionvalues$nbRowsExceededWarningMessage<-data$nbRowsExceededWarningMessage
      })
    }
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
    variantsGroup<-read.table(("filterVariant.csv"),header=T,stringsAsFactors=F,colClasses=c("character","character"))
    selectInput('selectSampleGroup1', 'Control group', 
                choices = list("Groups"=variantsGroup$Name), 
                selected=variantsGroup$Name[1],
                selectize = FALSE)
  })
  
  output$selectSampleGroup2UI<-renderUI({
    variantsGroup<-read.table(("filterVariant.csv"),header=T,stringsAsFactors=F,colClasses=c("character","character"))
    selectInput('selectSampleGroup2', 'Case group', 
                choices = list("Groups"=variantsGroup$Name), 
                selected=variantsGroup$Name[1],
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
      variantsGroup<-read.table(("filterVariant.csv"),header=T,stringsAsFactors=F,colClasses=c("character","character"))
      
      isolate({
        sampleGroup1name<-input$selectSampleGroup1
        sampleGroup2name<-input$selectSampleGroup2
      })
      
      selectSampleGroupIndex1<-which(variantsGroup$Name==sampleGroup1name)
      selectSampleGroupIndex2<-which(variantsGroup$Name==sampleGroup2name)
      
      group1sql<-variantsGroup$SQL[selectSampleGroupIndex1]
      group2sql<-variantsGroup$SQL[selectSampleGroupIndex2]
      
      group1sql<-preprocSQL(group1sql)
      group2sql<-preprocSQL(group2sql)
      
      startCommand<-paste('spark-submit --master local --conf spark.eventLog.enabled=true GVR.py','"',analysisName,"\"",scope,scale,"\"",group1sql,"\"","\"",group2sql,'"')
      browser()
      system(startCommand)
      
    }  
  })
  
  
  ####################################################
  #Results explorer
  ####################################################
  
  observe({
    input$refreshResultsButton
    analysesFiles<-dir("analyses/","*.txt")
    sessionvalues$analysesNames<-as.vector(unlist(sapply(analysesFiles,strsplit,'.txt')))
  })
  
  output$selectAnalysisUI<-renderUI({
    selectInput('selectAnalysis', 'Select analysis', choices = list(
      "Available analyses" = sessionvalues$analysesNames
    ), selected=sessionvalues$analysesNames[1],selectize = FALSE)
  })
  
  output$showVarResultsUI<-renderUI({
    nameAnalysis<-input$selectAnalysis
    
    if (length(sessionvalues$results)>0) {
      
      if (sessionvalues$results$scale=="variant") {
        if (sessionvalues$results$scope=="monogenic") {
          niceNames<-as.vector(sapply(colnames(sessionvalues$results$scoreSummary),idToName))
          initialSelect<-niceNames[c(1:5,8)]
        }
      }
      
      if (sessionvalues$results$scale=="gene") {
        if (sessionvalues$results$scope=="monogenic") {
          niceNames<-as.vector(sapply(colnames(sessionvalues$results$scoreSummary),idToName))
          initialSelect<-niceNames[1:4]
        }
        if (sessionvalues$results$scope=="digenic") {
          niceNames<-as.vector(sapply(colnames(sessionvalues$results$scoreSummary),idToName))
          initialSelect<-niceNames[c(1:5)]
        }
      }
      selectInput('showVarResults', 'Select variables to display', niceNames, 
                  selected=initialSelect,multiple=TRUE, selectize=TRUE,width='1050px')
    }
  })
  
  
  output$resultsTable<-DT::renderDataTable({
    if ((length(input$showVarResults)>0) & (length(sessionvalues$results)>0)) {
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
      sessionvalues$results<-procRes(fromJSON(txt=paste0("analyses/",nameAnalysis,".txt")))
      if (length(input$resultsTable_rows_selected)) {
        withProgress(min=1, max=4, expr={
          setProgress(message = 'Retrieving control data, please wait...',
                      value=2)
          
          if (sessionvalues$results$scale=="variant") {
            if (sessionvalues$results$scope=="monogenic") {
              variantsData<-sessionvalues$results$scoreSummary[input$resultsTable_rows_selected,]
              sqlControl<-sessionvalues$variantsSets[which(sessionvalues$variantsSets[,1]==sessionvalues$results$group1name),2]
              sqlControl<-paste0(sqlControl," and chr='",variantsData['Chr'],"' and position=",variantsData['Position']," and reference='",variantsData['Reference'],"' and alternative='",variantsData['Alternative'],"'")
              variantsControl<-loadData(sqlControl)
              setProgress(message = 'Retrieving case data, please wait...',
                          value=3)
              sqlCase<-sessionvalues$variantsSets[which(sessionvalues$variantsSets[,1]==sessionvalues$results$group2name),2]
              sqlCase<-paste0(sqlCase," and chr='",variantsData['Chr'],"' and position=",variantsData['Position']," and reference='",variantsData['Reference'],"' and alternative='",variantsData['Alternative'],"'")
              variantsCase<-loadData(sqlCase)
              variants<-rbind(variantsControl$data,variantsCase$data)
              variants<-cbind("Group"=c(rep(sessionvalues$results$group1name,nrow(variantsControl$data)),rep(sessionvalues$results$group2name,nrow(variantsCase$data))),variants)
              sessionvalues$variantDataGene<-variants
            }
          }
          
          if (sessionvalues$results$scale=="gene") {
            if (sessionvalues$results$scope=="monogenic") {
              geneID<-sessionvalues$results$scoreSummary[input$resultsTable_rows_selected,'Gene_Symbol']
              geneID<-strsplit(geneID,'<|>')[[1]][3]
              sqlControl<-sessionvalues$variantsSets[which(sessionvalues$variantsSets[,1]==sessionvalues$results$group1name),2]
              sqlControl<-paste0(sqlControl," and gene_symbol='",geneID,"'")
              variantsControl<-loadData(sqlControl)
              setProgress(message = 'Retrieving case data, please wait...',
                          value=3)
              sqlCase<-sessionvalues$variantsSets[which(sessionvalues$variantsSets[,1]==sessionvalues$results$group2name),2]
              sqlCase<-paste0(sqlCase," and gene_symbol='",geneID,"'")
              variantsCase<-loadData(sqlCase)
              variants<-rbind(variantsControl$data,variantsCase$data)
              variants<-cbind("Group"=c(rep(sessionvalues$results$group1name,nrow(variantsControl$data)),rep(sessionvalues$results$group2name,nrow(variantsCase$data))),variants)
              sessionvalues$variantDataGene<-variants
            }
            if (sessionvalues$results$scope=="digenic") {
              geneID1<-sessionvalues$results$scoreSummary[input$resultsTable_rows_selected,'Gene_Symbol1']
              geneID1<-strsplit(geneID1[[1]],'<|>')[[1]][3]
              geneID2<-sessionvalues$results$scoreSummary[input$resultsTable_rows_selected,'Gene_Symbol2']
              geneID2<-strsplit(geneID2[[1]],'<|>')[[1]][3]
              sqlControl<-sessionvalues$variantsSets[which(sessionvalues$variantsSets[,1]==sessionvalues$results$group1name),2]
              sqlControl<-paste0(sqlControl," and (gene_symbol='",geneID1,"' or gene_symbol='",geneID2,"')")
              variantsControl<-loadData(sqlControl)
              setProgress(message = 'Retrieving case data, please wait...',
                          value=3)
              sqlCase<-sessionvalues$variantsSets[which(sessionvalues$variantsSets[,1]==sessionvalues$results$group2name),2]
              sqlCase<-paste0(sqlCase," and (gene_symbol='",geneID1,"' or gene_symbol='",geneID2,"')")
              variantsCase<-loadData(sqlCase)
              variants<-rbind(variantsControl$data,variantsCase$data)
              variants<-cbind("Group"=c(rep(sessionvalues$results$group1name,nrow(variantsControl$data)),rep(sessionvalues$results$group2name,nrow(variantsCase$data))),variants)
              sessionvalues$variantDataGene<-variants
            }
          }
        }
        )
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
    }
  })
  
  output$showVarMetadataUI<-renderUI({
    nameAnalysis<-input$selectAnalysis
    if (length(sessionvalues$variantDataGene)>0) {
      niceNames<-as.vector(sapply(colnames(sessionvalues$variantDataGene),idToName))
      selectInput('showVarMetadata', 'Select variables to display', niceNames, 
                  selected=niceNames[c(1:8,24:26,17)],multiple=TRUE, selectize=TRUE,width='1050px')
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
             sessionvalues$results$group1name,br(),
             sessionvalues$results$group2name,br(),
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
    )
    
  })
  
  
})

