require(RMySQL)
require(reshape2)
require(plyr)

createGenotypeMatrix<-function() {
  
  groupsdb <- dbConnect(RSQLite::SQLite(), "samplesSets.db")
  
  groups<-dbReadTable(groupsdb,"samplesSets")
  
  controlQuery<-groups[which(groups$group=="ULBControlRareVariants"),]$sql
  caseQuery<-groups[which(groups$group=="ULBNeuroDevRareVariants"),]$sql
  
  datadb <- dbConnect(RSQLite::SQLite(), "data.db")
  system.time(control<-dbGetQuery(datadb,paste0("select Gene_Ensembl,ID,Sample_ID,Zygosity from variants where ",controlQuery)))
  system.time(case<-dbGetQuery(datadb,paste0("select Gene_Ensembl,ID,Sample_ID,Zygosity from variants where ",caseQuery)))
  dbDisconnect(datadb)
  
  patientID<-read.table('mappingZH_ISDBM_withFMC.txt',stringsAsFactors=F)[,1]
  
  DF<-rbind(control,case)
  genoMat<-dcast(DF,Gene_Ensembl+ID~Sample_ID,fill=0)
  genoMat<-genoMat[,c(1,2,match(patientID,colnames(genoMat)))] 
  
  write.table(file="genoMat.txt",genoMat,col.names=F,row.names=F,quote=F)
  
}


