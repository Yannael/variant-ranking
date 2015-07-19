require("rvest")
require("magrittr")
require("RMySQL")
require("jsonlite")

username<-'yleborgn'

load("myResults.Rdata")

get4<-function(l) {l[[4]]}
get1<-function(l) {l[[1]]}
get2<-function(l) {l[[2]]}

procRes<-function(results) {
  res<-list()
  res$name<-results[[1]]
  res$type<-results[[2]]
  res$genotypes<-lapply(results[[3]],get4)
  res$scores<-sapply(results[[3]],get2)
  
  res$locus<-sapply(results[[3]],get1)
  
  uniqueid1<-apply(res$locus[2:4,],2,paste,collapse=":")
  to.keep<-match(uniqueid1,dbvariants$uniqueid)
  res$infovariants1<-dbvariants[to.keep,]
  
  scoreSummary1<-cbind(res$locus[2,],res$infovariants1[,'reference'],res$infovariants1[,'alternative'])
  colnames(scoreSummary1)<-c("Locus","Reference", "Alternative")
  
  if (results[[2]]=="singleVariant") {
    scoreSummary1<-cbind(scoreSummary1,Score=res$infovariants1[,'gene_symbol'])
  }

  scoreSummary2<-NULL
  if (results[[2]]=="pairVariantsMonogenic") {
    uniqueid2<-apply(res$locus[5:7,],2,paste,collapse=":")
    to.keep<-match(uniqueid2,dbvariants$uniqueid)
    res$infovariants2<-dbvariants[to.keep,]
    
    scoreSummary2<-cbind(res$locus[5,],res$infovariants2[,'reference'],res$infovariants2[,'alternative'],res$infovariants1[,'gene_symbol'])
    colnames(scoreSummary2)<-c("Locus2","Reference2", "Alternative2","Gene symbol")
  }
  
  if (results[[2]]=="pairVariantsDigenic") {
    uniqueid2<-apply(res$locus[6:8,],2,paste,collapse=":")
    to.keep<-match(uniqueid2,dbvariants$uniqueid)
    res$infovariants2<-dbvariants[to.keep,]
    scoreSummary2<-cbind(res$infovariants1[,'gene_symbol'],res$locus[6,],res$infovariants2[,'reference'],res$infovariants2[,'alternative'],res$infovariants2[,'gene_symbol'])
    colnames(scoreSummary2)<-c("Gene symbol1","Locus2","Reference2", "Alternative2","Gene symbol2")
    
  }
  
  res$scoreSummary<-cbind(scoreSummary1,scoreSummary2,Score=res$scores)
  
  res
}

simpleCap <- function(x) {
  s <- strsplit(x, "_")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=" ")
}


getFiltersFromTable<-function(data) {
  filters<-list()
  
  namesCol<-colnames(data)
  niceNames<-as.vector(sapply(namesCol,simpleCap))
  
  for (i in 1:ncol(data)) {
    
    filterCol<-
      switch(class(data[,i]),
             character={
               if (length(unique(data[,i]))>50) {
                 list(
                   id= tolower(gsub(" ","",namesCol[i])),
                   label= niceNames[i],
                   type= 'string',
                   default_value=data[1,i],
                   operators=list('equal','not_equal','contains', 'is_empty', 'is_not_empty'))
               }
               else {
                 values<-setdiff(unique(data[,i]),"")
                 list(
                   id= tolower(gsub(" ","",namesCol[i])),
                   label= niceNames[i],
                   type= 'string',
                   input='select',
                   values=values,
                   default_value=values[1],
                   operators=list('equal','not_equal','contains', 'is_empty', 'is_not_empty'))
               }
             },
             integer=list(
               id= tolower(gsub(" ","",namesCol[i])),
               label= niceNames[i],
               type= 'integer',
               default_value=0,
               operators=list('equal','not_equal',  'less', 'less_or_equal', 'greater','greater_or_equal','between')),
             numeric=list(
               id= tolower(gsub(" ","",namesCol[i])),
               label= niceNames[i],
               type= 'double',
               default_value=0,
               validation=list(
                 min= 0,
                 step= 0.01
               ),
               operators=list('equal','not_equal',  'less', 'less_or_equal', 'greater','greater_or_equal','between'))
      )
    filters<-c(filters,list(filterCol))
  }
  filters
}

getWidgetTable<-function(data,session) {
  action <- dataTableAjax(session, data,rownames=F)
  widget<-datatable(data, 
            extensions = 'Scroller',
            server = TRUE, 
            selection = 'single',
            rownames=F,
            #filter = 'top',
            escape=T,
            options = list(
              dom= 'tiS',
              deferRender = TRUE,
              scrollY = 335,
              ajax = list(url = action),
              columnDefs = list(
                list(
                  targets = c("_all"),
                  render = JS(
                    "function(data, type, row, meta) {",
                    "return type === 'display' && data.length > 15 ?",
                    "'<span title=\"' + data + '\">' + data.substr(0, 11) + '...</span>' : data;",
                    "}")
                ),
                list(className="dt-right",targets="_all")
              )
            )
  )
  widget
}

con <- dbConnect(RSQLite::SQLite(), "../groupsToComparePartial.db")
rs<-dbSendQuery(con,"select * from dbvariants" )
dbvariants<-fetch(rs, n=-1)

results<-fromJSON(txt="singVarDeNovo.txt")
results1<-list(procRes(results))
results<-fromJSON(txt="pairVarCompound.txt")
results2<-list(procRes(results))
results<-fromJSON(txt="pairVarDigenic.txt")
results3<-list(procRes(results))
results<-list(results1,results2,results3)
names(results)<-c("Trios_De_Novo","Trios_Compound_Heterozygous","Trios_Digenic")

dbDisconnect(con)

variantsdb <- dbConnect(RSQLite::SQLite(), "../variants.db")
variants<-dbReadTable(variantsdb,"variants")
dbDisconnect(variantsdb)

loadGroups<-function() {
groupsdb<-dbConnect(RSQLite::SQLite(), "../groups.db")
groups<-dbReadTable(groupsdb,"groups")
dbDisconnect(groupsdb)
groups
}

loadPhenotypes<-function(sql) {
  if (sql!="") sql<-paste0("where ",sql)
  phenotypesdb<-dbConnect(RSQLite::SQLite(), "../phenotypes.db")
  phenotypes<-dbGetQuery(phenotypesdb,paste0("select * from phenotypes ",sql))
  dbDisconnect(phenotypesdb)
  phenotypes
}


