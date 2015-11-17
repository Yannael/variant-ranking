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
                                   div(actionButton('connectCliniPhenome', label = "Connect CliniPhenome",class = NULL),
                                       align="right"),
                                   bsModal("cliniPhenomeBS", "Answer from CliniPhenome", "connectCliniPhenome", 
                                           size = "large",htmlOutput('answerCliniPhenome')
                                   ),
                                   uiOutput("filterPhenotype"),
                                   div(
                                     actionButton("getIDButtonPhenotypesGroup", label = "Get sample IDs"),
                                    downloadButton('downloadPhenotypesSelection', label = "Download selection (CSV)",class = NULL),
                                       align="right"),
                                   bsModal("getIDPhenotypesGroup", "List of sample IDs", "getIDButtonPhenotypesGroup", 
                                           size = "large",textOutput('listSamplesIDs')
                                   ),
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
                                 uiOutput("filterVariant"),
                                 h5(textOutput("nbRowsExceededWarningMessage")),
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



