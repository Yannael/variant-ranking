require(RMySQL)
require(reshape2)
require(plyr)

createGenotypeMatrix<-function(analysisName) {
  
  #analysisName="ULB_Control_vs_NeuroDev_Damaging_Gene"
  #analysisName="tol_ctrl_rare"
  #analysisName="control_tol_rare"
  
  load(paste("analyses/",analysisName,".Rdata",sep=""))
  
  DF<-analysis$variantData
  DF<-DF[,c('ID','Sample_ID','Gene_Ensembl','Zygosity')]
  samplesID<-analysis$samplesID
  
  #Trios
  #samplesID<-read.table('mappingZH_ISDBM_withFMC.txt',stringsAsFactors=F)[,1]
  
  genoMat<-dcast(DF,Gene_Ensembl+ID~Sample_ID,fill=0)
  genoMat<-genoMat[,c(1,2,match(samplesID,colnames(genoMat)))] 
  
  write.table(file=paste("analyses/",analysisName,"_genoMat.txt",sep=""),genoMat,col.names=F,row.names=F,quote=F)
  
}

andrea<-function() {
  
  genes<-read.table("listGenesAndrea.txt",stringsAsFactors=F)[,1]
  genesChar<-paste("'",genes,"'",sep="",collapse=",")
  
}
