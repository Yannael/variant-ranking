require("rvest")
require("magrittr")
require("RMySQL")
require("jsonlite")

load(file="phenotypes.Rdata")
phenotypesAll<-rbind(phenotypesULB,phenotypes1000Gen)

username<-'yleborgn'

load("mySampleGroups.Rdata")  

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


con <- dbConnect(RSQLite::SQLite(), "../groupsToComparePartial.db")
rs<-dbSendQuery(con,"select * from dbvariants" )
dbvariants<-fetch(rs, n=-1)

results<-fromJSON(txt="singVarDeNovo2.txt")
results1<-list(procRes(results))
results<-fromJSON(txt="pairVarCompound2.txt")
results2<-list(procRes(results))
results<-fromJSON(txt="pairVarDigenic.txt")
results3<-list(procRes(results))
results<-list(results1,results2,results3)
names(results)<-c("Trios_De_Novo","Trios_Compound_Heterozygous","Trios_Digenic")

dbDisconnect(con)

