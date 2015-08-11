library(shiny)
library(DT)
library(queryBuildR)
library(shinyBS)

shinyUI(fluidPage(
  includeCSS('www/style.css'),
  #bootstrapPage(),
  #includeCSS('www/bootstrap.min.css'),
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
             tabPanel("Sample set manager", 
                      tags$div(class="extraspace2"),
                      fluidRow(
                        column(10,offset=1,
                               uiOutput("selectSampleGroupUI"),
                               queryBuildROutput("queryBuilderSamples",width="800px",height="100%"),
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
                        column(3,
                               uiOutput("showVarPhenotypeUI")
                        ),
                        column(9,
                               DT::dataTableOutput('phenotypesTable'),
                               tags$div(class="extraspace1")
                        )
                      )
             ),
             tabPanel("Variant set manager", 
                      tags$div(class="extraspace2"),
                      fluidRow(
                        column(10,offset=1,
                               uiOutput("selectVariantGroupUI"),
                               queryBuildROutput("queryBuilderVariants",width="800px",height="100%"),
                               actionButton("variantsQueryApply", label = "Apply filters"),
                               tags$script('
                                      function sqlQueryVariantsFunction() {
                                      var sqlQueryVariants = $("#queryBuilderVariants").queryBuilder("getSQL", false);
                                      Shiny.onInputChange("sqlQueryVariantsValue", sqlQueryVariants);
                                      };
                                      document.getElementById("variantsQueryApply").onclick = function() {sqlQueryVariantsFunction()}
                                      document.getElementById("variantsQuerySave").onclick = function() {sqlQueryVariantsFunction()}
                                 '),
                               tags$script('            
                                      Shiny.addCustomMessageHandler("callbackHandlerSelectVariantGroup",  function(sqlQuery) {
                                           if (sqlQuery=="") $("#queryBuilderVariants").queryBuilder("reset")
                                           else $("#queryBuilderVariants").queryBuilder("setRulesFromSQL",sqlQuery);
                                      });
                                '),
                               actionButton("variantsQuerySave", label = "Save group"),
                               bsModal("modalVariantsQuerySave", "Save group", "variantsQuerySave", 
                                       size = "small",
                                       textInput("variantGoupNameSave", "Save group as :", value = ""),
                                       actionButton("variantsQuerySave2", label = "Save")
                               ),
                               actionButton("deleteButtonVariantGroup", label = "Delete group"),
                               bsModal("deleteConfirmVariantGroup", "Are you sure?", "deleteButtonVariantGroup", 
                                       size = "small",
                                       actionButton("deleteConfirmYesButtonVariantGroup", label="Yes")
                               )
                        )
                      ),
                      hr(),
                      fluidRow(
                        column(3,
                               uiOutput("showVarVariantsUI")
                        ),
                        column(9,
                               DT::dataTableOutput('variantsTable'),
                               tags$div(class="extraspace1")
                        )
                      )#,
                      #tags$div(class="extraspace1")
             ),
             tabPanel("Ranking engine", 
                      tags$div(class="extraspace2"),
                      fluidRow(
                        column(3,
                               h3("1) Sample groups"),
                               selectInput('selectedSampleGroup1', 'Control sample group', selected="1000genomes_EUR",
                                           choices = list(
                                             "General sample groups" = c('ALL' = 'all', "Erasme" = 'erasme', "1000 genomes"="1000genomes"),
                                             "Saved sample groups" = ""
                                           ), 
                                           selectize = FALSE),
                               selectInput('selectedSampleGroup2', 'Pathological sample group',selected="Erasme_Hydrocephalus",
                                           choices = list(
                                             "General sample groups" = c('ALL' = 'all', "Erasme" = 'erasme', "1000 genomes"="1000genomes"),
                                             "Saved sample groups" = ""
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


