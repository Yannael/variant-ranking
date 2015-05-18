library(d3heatmap)

load("../res.Rdata")

mm<-as.matrix(res.rs[[2]])

  mydist<-function(m) {
    dd<-dist(m,method="manhattan")
    dd[which(is.na(dd))]<-20
    dd
  }
  
#d3heatmap(mm,heatmap_options = list(dist=mydist),scale="none")
hm<-heatmap(mm,dist=mydist,scale="none",keep.dendro=T)
