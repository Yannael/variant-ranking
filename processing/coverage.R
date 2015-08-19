require(RMySQL)
require(reshape2)
require(plyr)

createGenotypeMatrix<-function() {
  
  analysisName="ULB_Control_vs_NeuroDev_Damaging"
  
  load("analyses.Rdata")
  
  DF<-analyses[[analysisName]]$variantData
  DF<-DF[,c('ID','Sample_ID','Gene_Ensembl','Zygosity')]
  samplesID<-analyses[[analysisName]]$samplesID
  
  #Trios
  samplesID<-read.table('mappingZH_ISDBM_withFMC.txt',stringsAsFactors=F)[,1]
  
  genoMat<-dcast(DF,Gene_Ensembl+ID~Sample_ID,fill=0)
  genoMat<-genoMat[,c(1,2,match(samplesID,colnames(genoMat)))] 
  
  write.table(file="genoMatTriosDamaging.txt",genoMat,col.names=F,row.names=F,quote=F)
  
}


