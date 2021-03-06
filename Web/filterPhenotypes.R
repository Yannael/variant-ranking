createFilterPhenotype<-function(input,output,session) {
  
  if (!file.exists("filterPhenotype.csv")) {
    data<-rbind(c("All",""))
    colnames(data)<-c("Name","SQL")
    write.table(file="filterPhenotype.csv",data,row.names=F,sep="\t")
  }
  
  output$filterPhenotypeQueryBuilderWidget<-renderQueryBuildR({
    load("filterPhenotypeSpec.Rdata")
    rules<-NULL
    queryBuildR(rules,filters)
  })
  
  output$filterPhenotypeSelectLoadUI<-renderUI({
    input$filterPhenotypeConfirmSave
    input$filterPhenotypeConfirmDelete
    data<-read.table(("filterPhenotype.csv"),header=T,stringsAsFactors=F,colClasses=c("character","character"))
    selectInput('filterPhenotypeSelectLoad', 'Select filter', 
                choices = data$Name, 
                selected=data$Name[1],
                selectize = FALSE)
  })
  
  output$filterPhenotypeSelectDeleteUI<-renderUI({
    input$filterPhenotypeConfirmSave
    input$filterPhenotypeConfirmDelete
    data<-read.table(("filterPhenotype.csv"),header=T,stringsAsFactors=F,colClasses=c("character","character"))
    selectInput('filterPhenotypeSelectDelete', 'Select filter', 
                choices = data$Name[-1], 
                selected=data$Name[2],
                multiple=T,
                selectize = FALSE)
  })
  
  #Confirm load
  observe({
    input$filterPhenotypeConfirmLoad
    isolate({
      if (length(input$filterPhenotypeSelectLoad)>0) {
        data<-read.table(("filterPhenotype.csv"),header=T,stringsAsFactors=F,colClasses=c("character","character"))
        selected<-input$filterPhenotypeSelectLoad
        sqlQuery<-data$SQL[which(data$Name==selected)]
        if (sqlQuery=="") sqlQuery<-"reset"
        session$sendCustomMessage(type='filterPhenotypeCallbackHandlerLoad', sqlQuery)
        toggleModal(session, "filterPhenotypeModalLoad", toggle = "close")
      }
    })
  })
  
  #Confirm Save
  observe({
    input$filterPhenotypeConfirmSave
    isolate({
      if ((length(input$filterPhenotypeNameSave)>0) & length(input$filterPhenotypeQueryBuilderSQL)>0) {
        data<-data.frame(input$filterPhenotypeNameSave,input$filterPhenotypeQueryBuilderSQL)
        colnames(data)<-c("group","sql")
        write.table(file=("filterPhenotype.csv"),data,append=T,row.names=F,col.names=F,sep="\t")
        toggleModal(session, "filterPhenotypeModalSave", toggle = "close")
      }
    })
  })
  
  #Confirm Delete
  observe({
    input$filterPhenotypeConfirmDelete
    isolate({
      if (length(input$filterPhenotypeSelectDelete)>0) {
        data<-read.table(("filterPhenotype.csv"),header=T,stringsAsFactors=F,colClasses=c("character","character"))
        selected<-input$filterPhenotypeSelectDelete
        data<-data[-match(selected,data$Name),]
        colnames(data)<-c("Name","SQL")
        write.table(file=("filterPhenotypeSpec.csv"),data,row.names=F,col.names=T,sep="\t")
        toggleModal(session, "filterPhenotypeModalDelete", toggle = "close")
      }
    })
  })
  
  
  
  output$filterPhenotype<-renderUI({
    fluidRow(column(12,
      actionButton("filterPhenotypeLoad", label = "Load filter"),  
      actionButton("filterPhenotypeSave", label = "Save filter"),
      actionButton("filterPhenotypeDelete", label = "Delete filter(s)"),
      tags$div(class="extraspace2"),
      queryBuildROutput("filterPhenotypeQueryBuilderWidget",width="970px",height="100%"),
      tags$div(class="extraspace2"),
      actionButton("filterPhenotypeApply", label = "Apply filter"),
      tags$script('
                  function filterPhenotypeGetSQLStatement() {
                  var sql = $("#filterPhenotypeQueryBuilderWidget").queryBuilder("getSQL", false);
                  Shiny.onInputChange("filterPhenotypeQueryBuilderSQL", sql);
                  };
                  document.getElementById("filterPhenotypeApply").onclick = function() {filterPhenotypeGetSQLStatement()}
                  '),
      tags$script('            
                  Shiny.addCustomMessageHandler("filterPhenotypeCallbackHandlerLoad",  function(sqlQuery) {
                  if (sqlQuery=="reset") $("#filterPhenotypeQueryBuilderWidget").queryBuilder("reset")
                  else $("#filterPhenotypeQueryBuilderWidget").queryBuilder("setRulesFromSQL",sqlQuery);
                  });
                  '),
      tags$script('
                  document.getElementById("filterPhenotypeSave").onclick = function() {filterPhenotypeGetSQLStatement()}
                  '),
      bsModal("filterPhenotypeModalLoad", "Load filter", "filterPhenotypeLoad", 
              size = "small",
              uiOutput("filterPhenotypeSelectLoadUI"),
              actionButton("filterPhenotypeConfirmLoad", label = "Load")
      ),
      bsModal("filterPhenotypeModalSave", "Save filter", "filterPhenotypeSave", 
              size = "small",
              textInput("filterPhenotypeNameSave", "Save filter as:", value = ""),
              actionButton("filterPhenotypeConfirmSave", label = "Save")
      ),
      bsModal("filterPhenotypeModalDelete", "Delete filter", "filterPhenotypeDelete", 
              size = "small",
              uiOutput("filterPhenotypeSelectDeleteUI"),
              actionButton("filterPhenotypeConfirmDelete", label = "Delete")
      )
    )
    )
  })
  
  output
}
