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
                 tabPanel("Phenotype manager", 
                          tags$div(class="extraspace2"),
                          fluidRow(
                            column(12,
                                   uiOutput("selectPhenotypesGroupUI"),
                                   queryBuildROutput("queryBuilderPhenotypes",width="970px",height="100%"),
                                   actionButton("phenotypesQueryApply", label = "Apply filter"),
                                   tags$script('
                                               function sqlQueryPhenotypesFunction() {
                                               var sqlQuerySamples = $("#queryBuilderPhenotypes").queryBuilder("getSQL", false);
                                               Shiny.onInputChange("sqlQueryPhenotypesValue", sqlQuerySamples);
                                               };
                                               document.getElementById("phenotypesQueryApply").onclick = function() {sqlQueryPhenotypesFunction()};
                                               //   document.getElementById("phenotypesQuerySave").onclick = function() {sqlQueryPhenotypesFunction()}
                                               '),
                                   tags$script('            
                                               Shiny.addCustomMessageHandler("callbackHandlerSelectPhenotypesGroup",  function(sqlQuery) {
                                               if (sqlQuery=="reset") $("#queryBuilderPhenotypes").queryBuilder("reset")
                                               else $("#queryBuilderPhenotypes").queryBuilder("setRulesFromSQL",sqlQuery);
                                               });
                                               '),
                                   actionButton("phenotypesQuerySave", label = "Save filter"),
                                   tags$script('
                                               document.getElementById("phenotypesQuerySave").onclick = function() {sqlQueryPhenotypesFunction()}
                                               '),
                                   bsModal("modalPhenotypesQuerySave", "Save sample group", "phenotypesQuerySave", 
                                           size = "small",
                                           textInput("phenotypesGoupNameSave", "Save sample group as :", value = ""),
                                           actionButton("phenotypesQuerySave2", label = "Save")
                                   ),
                                   actionButton("deleteButtonPhenotypesGroup", label = "Delete sample group"),
                                   bsModal("deleteConfirmPhenotypesGroup", "Are you sure?", "deleteButtonPhenotypesGroup", 
                                           size = "small",
                                           actionButton("deleteConfirmYesButtonPhenotypesGroup", label="Yes")
                                   ),
                                   actionButton("getIDButtonPhenotypesGroup", label = "Get sample IDs"),
                                   bsModal("getIDPhenotypesGroup", "List of sample IDs", "getIDButtonPhenotypesGroup", 
                                           size = "large",textOutput('listSamplesIDs')
                                   )
                                   
                                   )
                                   ),
                          hr(),
                          fluidRow(
                            column(12,
                                   div(downloadButton('downloadPhenotypesSelection', label = "Download selection (CSV)",class = NULL),
                                       align="right"),
                                   uiOutput("showVarPhenotypesUI"),
                                   dataTableOutput('phenotypesTable'),
                                   hr(),
                                   h5(strong("Pivot table")),
                                   rpivotTableOutput("pivotTablePhenotypes"),
                                   tags$div(class="extraspace1")
                            )
                          )
               ),
               tabPanel("Gene & variant filtering manager", 
                        tags$div(class="extraspace2"),
                        fluidRow(
                          column(12,
                                 uiOutput("selectVariantsGroupUI"),
                                 queryBuildROutput("queryBuilderVariants",width="970px",height="100%"),
                                 actionButton("variantsQueryApply", label = "Apply filter"),
                                 tags$script('
                                             function sqlQueryVariantsFunction() {
                                             var sqlQueryVariants = $("#queryBuilderVariants").queryBuilder("getSQL", false);
                                             Shiny.onInputChange("sqlQueryVariantsValue", sqlQueryVariants);
                                             };
                                             document.getElementById("variantsQueryApply").onclick = function() {sqlQueryVariantsFunction()};
                                             //   document.getElementById("variantsQuerySave").onclick = function() {sqlQueryVariantsFunction()}
                                             '),
                                 tags$script('            
                                             Shiny.addCustomMessageHandler("callbackHandlerSelectVariantsGroup",  function(sqlQuery) {
                                             if (sqlQuery=="reset") $("#queryBuilderVariants").queryBuilder("reset")
                                             else $("#queryBuilderVariants").queryBuilder("setRulesFromSQL",sqlQuery);
                                             });
                                             '),
                                 actionButton("variantsQuerySave", label = "Save filter"),
                                 tags$script('
                                             document.getElementById("variantsQuerySave").onclick = function() {sqlQueryVariantsFunction()}
                                             '),
                                 bsModal("modalVariantsQuerySave", "Save variant group", "variantsQuerySave", 
                                         size = "small",
                                         textInput("variantsGoupNameSave", "Save variant group as :", value = ""),
                                         actionButton("variantsQuerySave2", label = "Save")
                                 ),
                                 actionButton("deleteButtonVariantsGroup", label = "Delete variant group"),
                                 bsModal("deleteConfirmVariantsGroup", "Are you sure?", "deleteButtonVariantsGroup", 
                                         size = "small",
                                         actionButton("deleteConfirmYesButtonVariantsGroup", label="Yes")
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
                                 div(downloadButton('downloadVariantsSelection', label = "Download selection (CSV)",class = NULL),
                                     align="right"),
                                 uiOutput("showVarVariantsUI"),
                                 dataTableOutput('variantsTable'),
                                 hr(),
                                 h5(strong("Pivot table")),
                                 rpivotTableOutput("pivotTableVariants"),
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
                                   radioButtons("rankingScope", "Scope",
                                                c("Monogenic" = "monogenic",
                                                  "Digenic" = "digenic"
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
                                   uiOutput("selectAnalysisUI"),
                                   actionButton("refreshResultsButton","Refresh")
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



