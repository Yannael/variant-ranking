require("rvest")
require("magrittr")
require("RMySQL")
require("rjson")

load(file="phenotypes.Rdata")
phenotypesAll<-rbind(phenotypesULB,phenotypes1000Gen)

username<-'yleborgn'

load("mySampleGroups.Rdata")  

load("myResults.Rdata")

get41<-function(l) {l[[4]][[1]]}
get42<-function(l) {l[[4]][[2]]}
get4<-function(l) {l[[4]]}
get1<-function(l) {l[[1]]}
get2<-function(l) {l[[2]]}

procRes<-function(results) {
  con <- dbConnect(RSQLite::SQLite(), "../groupsToComparePartial.db")
  if (results[[2]]=="singleVariant") {
    res<-list()
    res$name<-results[[1]]
    res$type<-"singleVariant"
    res$genotypes<-t(as.data.frame(sapply(results[[3]],get4)))
    res$locus<-sapply(results[[3]],get1)
    uniqueid<-apply(res$locus[2:4,],2,paste,collapse=":")
    dbvariants<-dbReadTable(con,"dbvariants")
    to.keep<-match(uniqueid,dbvariants$uniqueid)
    res$infovariants<-dbvariants[to.keep,]
    
    res$scores<-sapply(results[[3]],get2)
    
    scoreSummary<-cbind(res$locus[2,],res$infovariants[,'reference'],res$infovariants[,'alternative'],res$infovariants[,'gene_symbol'],res$score)
    colnames(scoreSummary)<-c("Locus","Reference", "Alternative","Gene symbol","Score")
    res$scoreSummary<-scoreSummary
    
  }
  if (results[[2]]=="pairVariantsMonogenic") {
    res<-list()
    res$name<-results[[1]]
    res$type<-"pairVariantsMonogenic"
    res$genotypes1<-t(as.data.frame(sapply(results[[3]],get41)))
    res$genotypes2<-t(as.data.frame(sapply(results[[3]],get42)))
    res$locus<-sapply(results[[3]],get1)
    uniqueid1<-apply(res$locus[2:4,],2,paste,collapse=":")
    dbvariants<-dbReadTable(con,"dbvariants")
    to.keep<-match(uniqueid1,dbvariants$uniqueid)
    res$infovariants1<-dbvariants[to.keep,]
    uniqueid2<-apply(res$locus[5:7,],2,paste,collapse=":")
    dbvariants<-dbReadTable(con,"dbvariants")
    to.keep<-match(uniqueid2,dbvariants$uniqueid)
    res$infovariants2<-dbvariants[to.keep,]
    
    res$scores<-sapply(results[[3]],get2)
    
    scoreSummary<-cbind(res$locus[2,],res$infovariants1[,'reference'],res$infovariants1[,'alternative'],res$locus[5,],res$infovariants2[,'reference'],res$infovariants2[,'alternative'],res$infovariants1[,'gene_symbol'],res$score)
    colnames(scoreSummary)<-c("Locus1","Reference1", "Alternative1","Locus2","Reference2", "Alternative2","Gene symbol","Score")
    res$scoreSummary<-scoreSummary
    
  }
  dbDisconnect(con)
  res
}


con <- dbConnect(RSQLite::SQLite(), "../groupsToComparePartial.db")
rs<-dbSendQuery(con,"select * from dbvariants" )
dbvariants<-fetch(rs, n=-1)
  
results<-fromJSON(file="singVarDeNovo2.txt")
results1<-list(procRes(results))
results<-fromJSON(file="pairVarCompound2.txt")
results2<-list(procRes(results))
results<-list(results1,results2)
names(results)<-c("Trios_De_Novo","Trios_Compound_Heterozygous")


