#require("rvest")
#require("magrittr")
require("RMySQL")
require("jsonlite")
require('rpivotTable')
require('plyr')
library(RJDBC)

source("../processing/coverage.R")

drv <- JDBC(driverClass = "org.apache.hive.jdbc.HiveDriver",
            classPath = list.files("/Users/yalb/Downloads/impala-jdbc-0.5-2",pattern="jar$",full.names=T),
            identifier.quote="`")


Sys.setenv(SPARK_HOME="/Users/yalb/spark")
Sys.setenv(PATH=paste0("/Users/yalb/spark/bin:/Users/yalb/spark/sbin:",Sys.getenv("PATH")))

username<-'yleborgn'

#load("myResults.Rdata")

getWidgetTable<-function(data,session,selection='none',targetsShort="_all") {
  action <- dataTableAjax(session, data,rownames=F)
  widget<-datatable(data, 
                    rownames=F,
                    escape=F,
                    selection = selection,
                    options = list(
                      ajax = list(url = action),
                      dom= 'lipt',
                      lengthMenu = list(c(10, 100, 1000), c('10', '100','1000')),pageLength = 10,
                      columnDefs = list(
                        list(
                          targets = c(1),
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

loadPhenotypes<-function(sql) {
  
  if (sql!="") {
    sql<-paste0(" where ",sql)
    sql<-gsub(',',"','",sql)
  }
  condb<-dbConnect(RSQLite::SQLite(), paste0("../phenotypes.db"))
  data<-dbGetQuery(condb,paste0("select * from phenotypes ",sql))
  dbDisconnect(condb)
  data
}


loadData<-function(sql,noLimit=F,excludeID=T) {
  
  if (sql!="") {
    sql<-paste0(" where ",sql)
    sql<-gsub(' *, *',"','",sql)
  }
  #condb<-dbConnect(RSQLite::SQLite(), paste0("../data.db"))
  #condb <- dbConnect(drv, "jdbc:hive2://192.168.99.101:21050/;auth=noSasl")
  condb <- dbConnect(drv, "jdbc:hive2://127.0.0.1:21050/;auth=noSasl")
  
  #nbrows<-dbGetQuery(condb,paste0("select count(*) from variants ",sql))
  nbrows<-1000
  if (noLimit) limit<-""
  else limit<-" limit 1000"
  nbRowsExceededWarningMessage<-""
  if (nbrows>1000) {
    nbRowsExceededWarningMessage<-paste0("Query returns ",nbrows," records. First 100000 retrieved.")
  }
  
  #data<-dbGetQuery(condb,paste0("select * from variants ",sql,limit))
  data<-dbGetQuery(condb,paste0("select * from gvr4.variants ",sql,limit))
  
  dbDisconnect(condb)
  results<-list(data=data,nbRowsExceededWarningMessage=nbRowsExceededWarningMessage)
  results
}


get4<-function(l) {l[[4]]}
get1<-function(l) {l[[1]]}
get2<-function(l) {l[[2]]}

procRes<-function(results) {
  res<-list()
  res$name<-results[[1]]
  res$scale<-results[[2]]
  res$scope<-results[[3]]
  res$start_time<-as.POSIXct(results[[4]],origin = "1970-01-01",tz="Europe/Brussels")
  res$end_time<-as.POSIXct(results[[5]],origin = "1970-01-01",tz="Europe/Brussels")
  res$run_time<-round(results[[6]],digits=2)
  res$scores<-sapply(results[[7]],get2)
  res$locus<-sapply(results[[7]],get1)
  res$caseSampleID=results[[8]]
  res$controlSampleID=results[[9]]
  res$group1name=results[[10]]
  res$group2name=results[[11]]
  if (res$scale=="variant") {
    if (res$scope=="monogenic") {
      variantsID<-res$locus[1,]
      scoreSummary<-cbind(t(data.frame(sapply(variantsID,strsplit,':'))),res$locus[2,])
      colnames(scoreSummary)<-c("Chr","Position","Reference","Alternative",'Gene_Symbol')
      rownames(scoreSummary)<-NULL
      scoreSummary[,'Gene_Symbol']<-paste0("<a href='http://www.ncbi.nlm.nih.gov/omim/?term=",scoreSummary[,'Gene_Symbol'],"' target='_blank'>",scoreSummary[,'Gene_Symbol'],"</a>")
      #if (results[[2]]=="pairVariantsMonogenic") {
      #  uniqueid2<-apply(res$locus[5:7,],2,paste,collapse=":")
      #  to.keep<-match(uniqueid2,dbvariants$uniqueid)
      #  res$infovariants2<-dbvariants[to.keep,]
      #  
      #  scoreSummary2<-cbind(res$locus[5,],res$infovariants2[,'reference'],res$infovariants2[,'alternative'],res$infovariants1[,'gene_symbol'])
      #  colnames(scoreSummary2)<-c("Locus2","Reference2", "Alternative2","Gene symbol")
      #}
    }
    
    if (results[[2]]=="pairVariantsDigenic") {
      uniqueid2<-apply(res$locus[6:8,],2,paste,collapse=":")
      to.keep<-match(uniqueid2,dbvariants$uniqueid)
      res$infovariants2<-dbvariants[to.keep,]
      scoreSummary2<-cbind(res$infovariants1[,'gene_symbol'],res$locus[6,],res$infovariants2[,'reference'],res$infovariants2[,'alternative'],res$infovariants2[,'gene_symbol'])
      colnames(scoreSummary2)<-c("Gene symbol1","Locus2","Reference2", "Alternative2","Gene symbol2")
      
    }

    res$scoreSummary<-cbind(Score=t(res$scores),scoreSummary)
    colnames(res$scoreSummary)[1:3]<-c("Score","Score_Case","Score_Control")
    
  }
  
  if (res$scale=="gene") {
    if (res$scope=="monogenic") {
      geneID<-res$locus
      geneID<-paste0("<a href='http://www.ncbi.nlm.nih.gov/omim/?term=",geneID,"' target='_blank'>",geneID,"</a>")
      res$scoreSummary<-cbind(t(res$scores),geneID)
      colnames(res$scoreSummary)<-c("Score","Score_Case","Score_Control","Gene_Symbol")
    }
    
    if (res$scope=="digenic") {
      genes<-ldply(res$scores[1,])
      res$scores<-res$scores[2:4,]
      geneID1<-genes[,1]
      geneID1<-paste0("<a href='http://www.ncbi.nlm.nih.gov/omim/?term=",geneID1,"' target='_blank'>",geneID1,"</a>")
      geneID2<-genes[,2]
      geneID2<-paste0("<a href='http://www.ncbi.nlm.nih.gov/omim/?term=",geneID2,"' target='_blank'>",geneID2,"</a>")
      res$scoreSummary<-cbind(t(res$scores),geneID1,geneID2)
      colnames(res$scoreSummary)<-c("Score","Score_Case","Score_Control","Gene_Symbol1","Gene_Symbol2")
    }
  }
  res
}

analysesFiles<-dir("analyses/","*.*")
analysesNames<-as.vector(unlist(sapply(analysesFiles,strsplit,'.txt')))

#resultsAll<-list()
#for (analysisName in analysesNames) {
#  results<-fromJSON(txt=paste0("analyses/",analysisName,".txt"))
#  resultsAll[[analysisName]]<-list(procRes(results,analysisName))
#}

