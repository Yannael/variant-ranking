library(shiny)
library(DT)
library(queryBuildR)

shinyUI(fluidPage(
  includeCSS('www/style.css'),
  includeCSS('www/query-builder.default.min.css'),
  #includeCSS('www/bootstrap.min.css'),
  #tags$head(tags$script(src="bootstrap.min.js")),
  #tags$head(tags$script(src="jquery.min.js")),
  #tags$head(tags$script(src="query-builder.js")),
  tags$head(tags$script(src="query-builder.standalone.js")),
  fluidRow(
           img(src="mgbck.jpg", height = 150, width = 1000)
           #headerPanel("Gene & Variant Ranking Toolbox")
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
             tabPanel("Sample group manager", 
                      tags$div(class="extraspace2"),
                      fluidRow(
                        column(3,
                               selectInput('selectedSampleGroup', 'Select sample group', choices = list(
                                 "General sample groups" = c('ALL' = 'all', "Erasme" = 'erasme', "1000 genomes"="genomes1000"),
                                 "Saved sample groups" = mySampleGroups[[1]]
                               ), selectize = FALSE)
                        )
                      ),
                      fluidRow(
                        column(7,
                               queryBuildROutput("queryBuilder",width="600px",height="100%"),
                               textOutput("sqlquery"),
                               actionButton("samplesQuery", label = "Apply filters"),
                               tags$script('
                                      document.getElementById("samplesQuery").onclick = function() {
                                      var sqlQuery = $("#queryBuilder").queryBuilder("getSQL", "question_mark");

                                      Shiny.onInputChange("sqlQuery", sqlQuery);
                                      };
                                 ')
                        )
                      ),
                      hr(),
                      fluidRow(
                        column(3,
                               checkboxGroupInput('showVarPhenotype', 'Display fields:',
                                                  names(phenotypesAll), selected = c("Data source","Sample ID","Pathology", "Control","Gender","Super population")),
                               hr(),
                               textInput("sampleGroupNameSave", label = "Save current sample group as:", 
                                         value = ""),
                               actionButton("saveFile", label = "Save current sample group")
                        ),
                        column(9,
                               div(
                                 DT::dataTableOutput('phenotypesTable'),
                                 tags$style(type="text/css", '.shiny-datatable-output tfoot {display:table-header-group;}')
                                 , style = 'width:690px;'),
                               tags$div(class="extraspace1")
                               
                        )
                        
                      )
             ),
             tabPanel("Filters", 
                      tags$div(class="extraspace2"),
                      fluidRow(
                        column(8,
                               queryBuildROutput("queryBuilderVariants",width="600px",height="100%"),
                               textOutput("sqlQueryVariants"),
                               actionButton("variantsQuery", label = "Apply filters"),
                               tags$script('
                                      document.getElementById("variantsQuery").onclick = function() {
                                      var sqlQueryVariants = $("#queryBuilderVariants").queryBuilder("getSQL", "question_mark");
                                      Shiny.onInputChange("sqlQueryVariantsValue", sqlQueryVariants);
                                      };
                                 ')
                        )
                      ),
                      tags$div(class="extraspace1")
             ),
             tabPanel("Ranking engine", 
                      tags$div(class="extraspace2"),
                      fluidRow(
                        column(3,
                               h3("1) Sample groups"),
                               selectInput('selectedSampleGroup1', 'Control sample group', selected="1000genomes_EUR",
                                           choices = list(
                                             "General sample groups" = c('ALL' = 'all', "Erasme" = 'erasme', "1000 genomes"="1000genomes"),
                                             "Saved sample groups" = mySampleGroups[[1]]
                                           ), 
                                           selectize = FALSE),
                               selectInput('selectedSampleGroup2', 'Pathological sample group',selected="Erasme_Hydrocephalus",
                                           choices = list(
                                             "General sample groups" = c('ALL' = 'all', "Erasme" = 'erasme', "1000 genomes"="1000genomes"),
                                             "Saved sample groups" = mySampleGroups[[1]]
                                           ), 
                                           selectize = FALSE)
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
                               textInput("resultName","Analysis name",""),
                               radioButtons("email", "Mail notification",
                                            c("Yes" = "yes",
                                              "No" = "no"
                                            ))
                        )
                      ),
                      fluidRow(
                        column(12,
                               hr(),
                               div(actionButton("NULL","Start analysis"),align="center"),
                               tags$div(class="extraspace1")
                               )
                      )
             ),
             tabPanel("Results explorer", 
                      tags$div(class="extraspace2"),
                      fluidRow(
                        column(3,
                               selectInput('selectedResultGroup', 'Select result ID', choices = list(
                                 "Available results" = names(results)
                               ), selected="Trios_Compound_Heterozygous",selectize = FALSE)
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
  
))


