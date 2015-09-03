#require("rvest")
#require("magrittr")
require("RMySQL")
require("jsonlite")
require('rpivotTable')
require('plyr')

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

loadData<-function(sql,noLimit=F,excludeID=T) {
  
  if (sql!="") {
    sql<-paste0(" where ",sql)
    sql<-gsub(',',"','",sql)
  }
  condb<-dbConnect(RSQLite::SQLite(), paste0("../data.db"))
  #nbrows<-dbGetQuery(condb,paste0("select count(*) from variants ",sql))
  nbrows<-100000
  if (noLimit) limit<-""
  else limit<-" limit 100000"
  nbRowsExceededWarningMessage<-""
  if (nbrows>100000) {
    nbRowsExceededWarningMessage<-paste0("Query returns ",nbrows," records. First 100000 retrieved.")
  }
  
  data<-dbGetQuery(condb,paste0("select * from variants ",sql,limit))
  dbDisconnect(condb)
  if (excludeID)
    results<-list(data=data[,-c(1)],nbRowsExceededWarningMessage=nbRowsExceededWarningMessage)
  else {
    results<-list(data=data,nbRowsExceededWarningMessage=nbRowsExceededWarningMessage)
  }
  results
}

loadSet<-function(db,sql) {
  
  if (sql!="") {
    sql<-paste0("where ",sql)
    sql<-gsub(',',"','",sql)
  }
  condb<-dbConnect(RSQLite::SQLite(), paste0("../",db,".db"))
  data<-dbGetQuery(condb,paste0("select * from ",db," ",sql," limit 100000"))
  dbDisconnect(condb)
  data
}

load("../filtersTypes.Rdata")

load("../analyses.Rdata")

analysesNames<-names(analyses)

get4<-function(l) {l[[4]]}
get1<-function(l) {l[[1]]}
get2<-function(l) {l[[2]]}

procRes<-function(results,analysisName) {
  res<-list()
  res$name<-results[[1]]
  res$type<-results[[2]]
  res$scores<-sapply(results[[3]],get2)
  
  res$locus<-sapply(results[[3]],get1)
  
  if (res$type=="singleVariant") {
    variantsID<-res$locus[2,]
    i.match<-match(variantsID,analyses[[analysisName]]$variantData$ID)
    scoreSummary1<-analyses[[analysisName]]$variantData[i.match,c(1,3:6,16:35)]
    scoreSummary1[,'Gene_Symbol']<-paste0("<a href='http://www.ncbi.nlm.nih.gov/omim/?term=",scoreSummary1[,'Gene_Symbol'],"' target='_blank'>",scoreSummary1[,'Gene_Symbol'],"</a>")
    #browser()
    if (results[[2]]=="singleVariant") {
      #scoreSummary1<-cbind(scoreSummary1,Score=res$infovariants1[,'gene_symbol'])
    }
    
    #scoreSummary2<-NULL
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
    
    res$scoreSummary<-cbind(Score=res$scores,scoreSummary1)
  }
  
  if (res$type=="geneUnivariate") {
    geneID<-res$locus
    i.match<-match(geneID,analyses[[analysisName]]$variantData$Gene_Ensembl)
    scoreSummary1<-analyses[[analysisName]]$variantData[i.match,c(17:20)]
    scoreSummary1[,'Gene_Symbol']<-paste0("<a href='http://www.ncbi.nlm.nih.gov/omim/?term=",scoreSummary1[,'Gene_Symbol'],"' target='_blank'>",scoreSummary1[,'Gene_Symbol'],"</a>")
    res$scoreSummary<-cbind(t(res$scores),scoreSummary1)
    colnames(res$scoreSummary)[1:3]<-c("Score","Score_Case","Score_Control")
  }
  
  if (res$type=="geneBivariate") {
    genes<-ldply(res$scores[1,])
    res$scores<-res$scores[2:4,]
    geneID<-genes[,1]
    i.match<-match(geneID,analyses[[analysisName]]$variantData$Gene_Ensembl)
    scoreSummary1<-analyses[[analysisName]]$variantData[i.match,c(17:20)]
    colnames(scoreSummary1)<-paste0(colnames(scoreSummary1),"1")
    geneID<-genes[,2]
    i.match<-match(geneID,analyses[[analysisName]]$variantData$Gene_Ensembl)
    scoreSummary2<-analyses[[analysisName]]$variantData[i.match,c(17:20)]
    colnames(scoreSummary2)<-paste0(colnames(scoreSummary2),"2")
    scoreSummary1[,'Gene_Symbol1']<-paste0("<a href='http://www.ncbi.nlm.nih.gov/omim/?term=",scoreSummary1[,'Gene_Symbol1'],"' target='_blank'>",scoreSummary1[,'Gene_Symbol1'],"</a>")
    scoreSummary2[,'Gene_Symbol2']<-paste0("<a href='http://www.ncbi.nlm.nih.gov/omim/?term=",scoreSummary2[,'Gene_Symbol2'],"' target='_blank'>",scoreSummary2[,'Gene_Symbol2'],"</a>")
    res$scoreSummary<-cbind(t(res$scores),scoreSummary1,scoreSummary2)
    colnames(res$scoreSummary)[1:3]<-c("Score","Score_Case","Score_Control")
  }
  
  res
}

resultsAll<-list()
for (analysisName in analysesNames) {
  results<-fromJSON(txt=paste0("../",analysisName,".txt"))
  resultsAll[[analysisName]]<-list(procRes(results,analysisName))
}

