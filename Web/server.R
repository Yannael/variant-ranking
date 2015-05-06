library(shiny)
requireNamespace('htmlwidgets')
library(DT)

shinyServer(function(input, output) {
  sessionvalues <- reactiveValues()
  sessionvalues$phenotypes<-phenotypesAll
  sessionvalues$result_filterNCBIGeneIDs<-resultsAll[[2]][[2]]
  
  getSampleIDFromGroup<-function(groupName) {
    #dataGroup<-hb.get('shinyapp_sample_groups','yleborgn')
    #browser()
    index.group<-which(mySampleGroups[[1]]==paste0(groupName))
    #hb.get('shinyapp_sample_groups','yleborgn')[[1]][[3]][[index.group]]
    mySampleGroups[[2]][[index.group]]
  }
  
  observe({
    input$result_filterNCBI
    result_filterNCBIterms<-isolate({input$result_filterNCBIterms})
    if (result_filterNCBIterms=="") {
      sessionvalues$result_filterNCBIGeneIDs<-resultsAll[[2]][[2]]
    }
    else {
      result_filterNCBIterms<-gsub(" +","+",result_filterNCBIterms)
      #browser()
      xmlResults<-html(paste0("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&retmax=1000&term=",result_filterNCBIterms))
      geneIDs<-html_nodes(xmlResults,'id') %>% html_text
      sessionvalues$result_filterNCBIGeneIDs<-geneIDs
    }
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
    datatable(
      data,rownames=checkboxRows(data),escape= -1,
      ,options = list(
        dom='fltip', 
        lengthMenu = list(c(10, 25, -1), c('10', '25','All')),pageLength = 10,
        autoWidth = FALSE
        #columns.width = list(list(width = "200px", width = "200px",
        #                          width = "200px", width = "30px"))#, bFilter=F)
      )
    )
  }
  )
  
  output$resultsTable<-renderDataTable({
    entriesToSelect<-intersect(sessionvalues$result_filterNCBIGeneIDs,resultsAll[[2]][[2]])
    i.entriesToSelect<-match(entriesToSelect,resultsAll[[2]][[2]])
    data<-as.data.frame(resultsAll[[2]][[1]][i.entriesToSelect,input$showVarResults])
    datatable(
      data,escape=F,
      ,options = list( 
        dom='fltip', 
        lengthMenu = list(c(10, 25, -1), c('10', '25','All')),pageLength = 10,
        autoWidth = FALSE
        #columns.width = list(list(width = "200px", width = "200px",
        #                          width = "200px", width = "30px"))#, bFilter=F)
      )
    )
  }
  )
  
  getOMIM_ID<-function(geneID) {
    #geneID<-60498
    #geneID<-100288069
    geneNCBIPage<-html(paste0("http://www.ncbi.nlm.nih.gov/gene/?term=",geneID))
    possibleFields<-geneNCBIPage %>% html_nodes('#summaryDl a , #summaryDl dt') %>% html_text 
    OMIM_field_i<-sapply(possibleFields, grepl,pattern='MIM:(.*)')
    OMIM_field_i<-which(OMIM_field_i==T)
    OMIM_ID<-""
    if (length(OMIM_field_i)>0) { 
      OMIM_ID<-sub("MIM:","",possibleFields[[OMIM_field_i]])
      OMIM_ID<-sub(";","",OMIM_ID)
    }
    OMIM_ID
    
  }
  
  function() {
    
    xmlResults<-html("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term=IgA+NEPHROPATHY&retmax=1000")
    geneIDs<-html_nodes(xmlResults,'id') %>% html_text
    
    OMIM_IDs<-c()
    tt<-system.time({
      for (i in 1:5) {
        OMIM_IDs<-c(OMIM_IDs,getOMIM_ID(geneIDs[i]))
      }
    }
    )
    #ncbi<-html("http://www.ncbi.nlm.nih.gov/gene/?term=IgA+nephropathy")
    #res<-ncbi %>% html_nodes('td')
    
    data<-res[(0:19)*5+1]
    
    
    
    
  }
  
  #   observe({
  #     input$saveFile
  #     isolate({
  #       name<-input$sampleGroupNameSave
  #       phenoID<-sessionvalues$phenotypes[,"Sample ID"]
  #       hb.insert("shinyapp_sample_groups",list(list("yleborgn",c(paste0("c:",name)), list(phenoID))))
  #     })
  #     #browser()
  #   })
  #   
  
})

