library(shiny)
library(DT)
library(queryBuildR)
library(shinyBS)

shinyUI(
  
  fluidPage(
    includeCSS('www/style.css'),
    div(
      fluidRow(
        img(src="mgbck2.jpg", height = 150, width = 1000)
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
                 tabPanel("Gene & variant filtering manager", 
                          tags$div(class="extraspace2"),
                          fluidRow(
                            column(12,
                                   uiOutput("selectSampleGroupUI"),
                                   queryBuildROutput("queryBuilderSamples",width="970px",height="100%"),
                                   actionButton("samplesQueryApply", label = "Apply filter"),
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
                                           if (sqlQuery=="reset") $("#queryBuilderSamples").queryBuilder("reset")
                                           else $("#queryBuilderSamples").queryBuilder("setRulesFromSQL",sqlQuery);
                                      });
                                '),
                                   actionButton("samplesQuerySave", label = "Save filter"),
                                   bsModal("modalSamplesQuerySave", "Save filter", "samplesQuerySave", 
                                           size = "small",
                                           textInput("sampleGoupNameSave", "Save filter as :", value = ""),
                                           actionButton("samplesQuerySave2", label = "Save")
                                   ),
                                   actionButton("deleteButtonSampleGroup", label = "Delete filter"),
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
                                   h3("1) Variants groups"),
                                   uiOutput("selectSampleGroup1UI"),
                                   uiOutput("selectSampleGroup2UI")
                            ),
                            column(4,offset=1,
                                   h3("2) Ranking parameters"),
                                   radioButtons("rankingScale", "Ranking scale",
                                                c("Gene" = "gene",
                                                  "Variant" = "variant"
                                                )),
                                   radioButtons("multivariateRankingRadio", "Interaction",
                                                c("Univariate" = "univariate",
                                                  "Bivariate" = "bivariate"
                                                )),
                                   checkboxGroupInput("rankingCriterion", "Scoring function",
                                                      c("Count" = "count",
                                                        "Odds ratio" = "oddsratio",
                                                        "Student p-value" = "pvalue",
                                                        "Minimum Redundancy Maximum Relevance" = "mrmr"
                                                      ),
                                                      selected=c("count","pvalue"))
                            ),
                            column(3,
                                   h3("3) Results collection"),
                                   textInput("analysisName","Analysis name",""),
                                   radioButtons("email", "Mail notification",
                                                c("Yes" = "yes",
                                                  "No" = "no"),
                                                selected=c("no"))
                            )
                          ),
                          hr(),
                          fluidRow(
                            column(2,offset=3,
                                   div(actionButton("estimateAnalysisTimeButton","Estimate analysis time"),align="center"),
                                   tags$div(class="extraspace1")
                            ),
                            column(2,
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



