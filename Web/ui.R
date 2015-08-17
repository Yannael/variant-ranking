library(shiny)
library(DT)
library(queryBuildR)
library(shinyBS)

shinyUI(
  
  fluidPage(
    includeCSS('www/style.css'),
    div(
      fluidRow(
        img(src="mgbck.jpg", height = 150, width = 1000)
      ),
      tags$div(class="extraspace2"),
      fluidRow(
        column(12,
               tabsetPanel(
                 #tabPanel("Home", 
                 #        fluidRow(
                 #          h3("Welcome to BRiDGEIris Variant ranking interface"),
                 #          "Create groups using the Create group panel"
                 #        )
                 #),
                 tabPanel("Gene & variant group manager", 
                          tags$div(class="extraspace2"),
                          fluidRow(
                            column(12,
                                   uiOutput("selectSampleGroupUI"),
                                   queryBuildROutput("queryBuilderSamples",width="970px",height="100%"),
                                   actionButton("samplesQueryApply", label = "Apply filters"),
                                   tags$script('
                                      function sqlQuerySamplesFunction() {
                                      var sqlQuerySamples = $("#queryBuilderSamples").queryBuilder("getSQL", false);
                                      Shiny.onInputChange("sqlQuerySamplesValue", sqlQuerySamples);
                                      };
                                      document.getElementById("samplesQueryApply").onclick = function() {sqlQuerySamplesFunction()}
                                      document.getElementById("samplesQuerySave").onclick = function() {sqlQuerySamplesFunction()}
                                 '),
                                   tags$script('            
                                      Shiny.addCustomMessageHandler("callbackHandlerSelectSampleGroup",  function(sqlQuery) {
                                           if (sqlQuery=="") $("#queryBuilderSamples").queryBuilder("reset")
                                           else $("#queryBuilderSamples").queryBuilder("setRulesFromSQL",sqlQuery);
                                      });
                                '),
                                   actionButton("samplesQuerySave", label = "Save group"),
                                   bsModal("modalSamplesQuerySave", "Save group", "samplesQuerySave", 
                                           size = "small",
                                           textInput("sampleGoupNameSave", "Save group as :", value = ""),
                                           actionButton("samplesQuerySave2", label = "Save")
                                   ),
                                   actionButton("deleteButtonSampleGroup", label = "Delete group"),
                                   bsModal("deleteConfirmSampleGroup", "Are you sure?", "deleteButtonSampleGroup", 
                                           size = "small",
                                           actionButton("deleteConfirmYesButtonSampleGroup", label="Yes")
                                   )
                                   
                            )
                          ),
                          hr(),
                          fluidRow(
                            column(12,
                                   h5(textOutput("nbRowsExceededWarningMessage"))
                            )
                          ),
                          fluidRow(
                            column(12,
                                   uiOutput("showVarPhenotypeUI"),
                                   dataTableOutput('phenotypesTable'),
                                   hr(),
                                   h5(strong("Pivot table")),
                                   rpivotTableOutput("pivotTable"),
                                   tags$div(class="extraspace1")
                            )
                          )
                 ),
                 tabPanel("Ranking engine", 
                          tags$div(class="extraspace2"),
                          fluidRow(
                            column(3,
                                   h3("1) Sample groups"),
                                   uiOutput("selectSampleGroup1UI"),
                                   uiOutput("selectSampleGroup2UI")
                            ),
                            column(4,offset=1,
                                   h3("2) Ranking parameters"),
                                   radioButtons("rankingScale", "Ranking scale",
                                                c("Gene" = "gene",
                                                  "Variant" = "variant"
                                                )),
                                   selectInput('rankEngine_selectedGeneListVariantSet', 'Gene/Variant set ',selected="All",
                                               choices = list(
                                                 "Genes" = c('All'='allg','GeneList1' = 'l1', "GeneList2" = 'l2', "GeneList3"="l3"),
                                                 "Variants" = c('All'='allv','VariantList1' = 'v1', "VariantList2" = 'l2')
                                               ), 
                                               selectize = FALSE),
                                   checkboxGroupInput("rankingCriterion", "Ranking criterion",
                                                      c("Univariate entropy" = "entropy",
                                                        "Univariate student p-value" = "pvalue",
                                                        "Minimum Redundancy Maximum Relevance" = "mrmr"
                                                      ),
                                                      selected=c("entropy","pvalue"))
                            ),
                            column(3,
                                   h3("3) Results collection"),
                                   textInput("analysisName","Analysis name",""),
                                   radioButtons("email", "Mail notification",
                                                c("Yes" = "yes",
                                                  "No" = "no"
                                                ))
                            )
                          ),
                          fluidRow(
                            column(12,
                                   hr(),
                                   div(actionButton("startAnalysisButton","Start analysis"),align="center"),
                                   tags$div(class="extraspace1")
                            )
                          )
                 ),
                 tabPanel("Results explorer", 
                          tags$div(class="extraspace2"),
                          fluidRow(
                            column(3,
                                   uiOutput("selectAnalysisUI")
                            )
                          ),
                          hr(),
                          fluidRow(
                            column(12,
                                   uiOutput("resultsPanel"),
                                   tags$div(class="extraspace1")
                            )
                            
                          )
                 )
               )
        )
      )
    )
  )
  
)



