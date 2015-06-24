require("rvest")
require("magrittr")


load(file="phenotypes.Rdata")
phenotypesAll<-rbind(phenotypesErasme,phenotypes1000Gen)

username<-'yleborgn'

load("mySampleGroups.Rdata")  

load("myResults.Rdata")

get4<-function(l) {l[[4]]}
get1<-function(l) {l[[1]]}
get2<-function(l) {l[[2]]}

procRes<-function(results) {
  con <- dbConnect(RSQLite::SQLite(), "../groupsToComparePartial.db")
  resAll<-list()
  if (results[[2]]=="singleVariant") {
    res<-list()
    res$name<-results[[1]]
    res$type<-"singleVariant"
    res$genotypes<-t(as.data.frame(sapply(results[[3]],get4)))
    res$locus<-sapply(results[[3]],get1)
    i<-1
    infSNPs<-NULL
    for (i in 1:10) {
      rs<-dbSendQuery(con,paste0("select * from dbSNPs where Locus='",res$locus[2,i],"' and reference='", 
                               res$locus[3,i],"' and alternative='",res$locus[4,i],"'"))
      infSNPs<-rbind(infSNPs,unique(fetch(rs, n=-1))[1,])
    }
    res$infoSNPs<-infSNPs
    res$scores<-sapply(results[[3]],get2)
    resAll<-c(resAll,list(res))
  }
  dbDisconnect(con)
  resAll
}


con <- dbConnect(RSQLite::SQLite(), "../groupsToComparePartial.db")
rs<-dbSendQuery(con,"select * from dbSNPs" )
dbSNPs<-fetch(rs, n=-1)
  
results<-fromJSON("singVarDeNovo.txt")

results<-procRes(results)
