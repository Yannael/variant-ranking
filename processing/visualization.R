library(d3heatmap)

load("res.Rdata")

patient_class<-colnames(snpsMatsub)
patient_class<-sapply(patient_class,strsplit,':')
patient_class<-unlist(patient_class,use.name=F)
patient_class<-matrix(patient_class,ncol(snpsMatsub),2,byrow=T)

controls<-which(patient_class[,2]=="Control")
cases<-which(patient_class[,2]=="Patho")

i.order<-sort(patient_class[,2],index=T)$ix
snpsMatsub<-snpsMatsub[,i.order]

sn<-snpsMatsub
sn[1,1]<-NA

d3heatmap(sn)

mm<-as.matrix(res.rs[[2]])

  mydist<-function(m) {
    dd<-dist(m,method="manhattan")
    dd[which(is.na(dd))]<-20
    dd
  }
  
d3heatmap(mm,heatmap_options = list(dist=mydist))

d3heatmap(mm)

hm<-heatmap(mm,dist=mydist,scale="none",keep.dendro=T)

hm<-heatmap(scale(mtcars),keep.dendro=T)


d3heatmap(scale(mtcars), colors = "Greens", theme = "dark")


d3heatmap(mm,heatmap_options = list(),scale="none")