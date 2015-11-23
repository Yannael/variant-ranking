require(RMySQL)


#filters='PASS' and
#read_depth>10 and
#mapping_quality_zero_reads=0 and
#downsampled='false' and
#allele_num=2 and
#genotype_quality>70 and
#snpeff_effect<>'SYNONYMOUS_CODING' and 

fieldsInfo<-c("patient","chr","pos","reference","alternative","zygosity","read_depth","allelic_depth_proportion_ref","allelic_depth_proportion_alt","genotype_quality",
              "filters","downsampled","mapping_quality_zero_reads","allele_num",
              "change_type", "gene_symbol", "gene_ensembl", "transcript_ensembl","num_genes","clinvar_rs", 
              "consensus_MAF","consensus_MAC", "dbsnp_id",
              "snpeff_effect", "snpeff_impact",
              "cadd_phred","cadd_raw","vest_score","pph2_hdiv_score", "pph2_hdiv_pred", 
              "pph2_hvar_score", "pph2_hvar_pred", "sift_score", "sift_pred", "short_tandem_repeat")

fieldsInfoNewNames<-c("Sample_ID","Chr","Position","Reference","Alternative",
                      "Zygosity","Read_Depth","Allelic_Depth_Ref_Prop","Allelic_Depth_Alt_Prop","Genotype_Quality"
                      ,"Filters","Downsampled","Mapping_Quality_Zero_Reads","Allele_Num",
                      "Change_Type", "Gene_Symbol", "Gene_Ensembl","Transcript_Ensembl", "Num_Genes","Clinvar_RS", 
                      "Consensus_MAF","Consensus_MAC", "dbsnp_id",
                      "Snpeff_Effect", "Snpeff_Impact",
                      "CADD_Phred","CADD_Raw","VEST_Score","PPH2_hdiv_score", "PPH2_hdiv_pred", 
                      "PPH2_hvar_score", "PPH2_hvar_pred", "SIFT_Score", "SIFT_Pred", "Short_Tandem_Repeat")

dummy<-function() {
  #To which DB to connect
  
  #Regular Highlander
  connectFile<-"../../connectHighlander.R"
  tablename<-"exomes_hc"
  
  #Highlander 1000 genomes
  connectFile<-"../../connectHighlander2.R"
  tablename<-"1000g"
  
}


createCopyHighlanderDB<-function() {
  
  source(connectFile)
  
  samples<-dbGetQuery(highlanderdb,paste0("select distinct patient from ",tablename))[,1]
  fields_select<-paste(unique(c(fieldsInfo)),collapse=",")
  
  #For 1000g
  subsamples<-paste0("('",paste(samples[1:100],sep="",collapse="','"),"')") #75s
  system.time(data <- dbGetQuery(highlanderdb,paste0("select ",fields_select," from ",tablename," where patient in ",subsamples)))
  
  #For ULB
  #Takes about 3 minutes
  system.time(data <- dbGetQuery(highlanderdb,paste0("select ",fields_select," from ",tablename)))
  uniqueID<-paste(data$chr, data$pos, data$reference, data$alternative, sep=":")
  
  #id.redundant<-tapply(data$gene_symbol,uniqueID,length)
  
  zygosity<-data$zygosity
  zygosity[zygosity=="Homozygous"]<-2
  zygosity[zygosity=="Heterozygous"]<-1
  data$zygosity<-zygosity
  data<-cbind(ID=uniqueID,data)
  
  #data[is.na(data)]<-''
  colnames(data)<-c(fieldsInfoNewNames)
  #data<-data[order(data[,3],data[,4]),]
  
  variantsdb <- dbConnect(RSQLite::SQLite(), "variants1000.db")
  dbWriteTable(variantsdb,"variants",data,overwrite=T,row.names=F)
  
  phenotypesdb<-dbConnect(RSQLite::SQLite(), paste0("phenotypes1000.db"))
  phenotypes<-dbReadTable(phenotypesdb,"phenotypes")
  dbWriteTable(variantsdb,"phenotypes",phenotypes,overwrite=T,row.names=F)
  
  #2 minutes with 15M entries
  system.time(data<-dbGetQuery(variantsdb,paste0("select * from variants,phenotypes where variants.Sample_ID=phenotypes.Sample_ID")))
  
  #Remove duplicate Sample_ID column
  data<-data[,-which(colnames(data)=="Sample_ID")[2]]
  
  
  for (i in 1:ncol(data)) {
    if (class(data[,i])=="character") data[is.na(data[,i]),i]<-''
  }
  data[is.na(data)]<-'\\N'
  
  data[,'Chr']<-factor(data[,'Chr'],levels=c(1:22,'X','Y'))
  data[,'Change_Type']<-as.factor(data[,'Change_Type'])
  data[,'Snpeff_Effect']<-as.factor(data[,'Snpeff_Effect'])
  data[,'Snpeff_Impact']<-as.factor(data[,'Snpeff_Impact'])
  data[,'SIFT_Pred']<-as.factor(data[,'SIFT_Pred'])
  data[,'PPH2_hvar_pred']<-as.factor(data[,'PPH2_hvar_pred'])
  data[,'PPH2_hdiv_pred']<-as.factor(data[,'PPH2_hdiv_pred'])
  data[,'Consensus_MAF']<-as.numeric(data[,'Consensus_MAF'])
  data[,'Consensus_MAC']<-as.numeric(data[,'Consensus_MAC'])
  data[,'CADD_Phred']<-as.numeric(data[,'CADD_Phred'])
  data[,'CADD_Raw']<-as.numeric(data[,'CADD_Raw'])
  data[,'VEST_Score']<-as.numeric(data[,'VEST_Score'])
  data[,'PPH2_hdiv_score']<-as.numeric(data[,'PPH2_hdiv_score'])
  data[,'PPH2_hvar_score']<-as.numeric(data[,'PPH2_hvar_score'])
  data[,'SIFT_Score']<-as.numeric(data[,'SIFT_Score'])
  data[,'Allelic_Depth_Ref_Prop']<-as.numeric(data[,'Allelic_Depth_Ref_Prop'])
  data[,'Allelic_Depth_Alt_Prop']<-as.numeric(data[,'Allelic_Depth_Alt_Prop'])
  data[,'Genotype_Quality']<-as.numeric(data[,'Genotype_Quality'])
  data[,'Mapping_Quality_Zero_Reads']<-as.numeric(data[,'Mapping_Quality_Zero_Reads'])
  data[,'Filters']<-as.factor(data[,'Filters'])
  
  data[,'Data_Source']<-as.factor(data[,'Data_Source'])
  data[,'Pathology']<-as.factor(data[,'Pathology'])
  data[,'Gender']<-as.factor(data[,'Gender'])
  data[,'Super_Population']<-as.factor(data[,'Super_Population'])
  data[,'Population']<-as.factor(data[,'Population'])
  
  filtersTypes<-getFiltersFromTable(data[,-1])
  save(file="filtersTypes.Rdata",filtersTypes)
  
  datadb<-dbConnect(RSQLite::SQLite(), paste0("data1000.db"))
  dbWriteTable(datadb,"variants",data,row.names=F,overwrite=T)
  dbDisconnect(datadb)
  write.table(file="variantsULB.csv",data,col.names=F,row.names=F,sep=",",quote=F)
  
  missingGenes<-datawhich(data[,'Gene_Ensembl']=='')
  ID<-data[missingGenes,'ID']
  subData<-data[which(data[,'ID'] %in% ID),c('ID','Gene_Ensembl')]
  tt<-tapply(subData[,1],subData[,1],length)
  
}

createCopyHighlanderDB2<-function() {
  
  source(connectFile)
  
  samples<-dbGetQuery(highlanderdb,paste0("select distinct patient from ",tablename))[,1]
  fields_select<-paste(unique(c(fieldsInfo)),collapse=",")
  
  #For 1000g
  subsamples<-paste0("('",paste(samples[1:100],sep="",collapse="','"),"')") #75s
  system.time(data <- dbGetQuery(highlanderdb,paste0("select ",fields_select," from ",tablename," where patient in ",subsamples)))
  
  #For ULB
  #Takes about 3 minutes
  system.time(data <- dbGetQuery(highlanderdb,paste0("select ",fields_select," from ",tablename)))
  
  for (i in 1:ncol(data)) {
    if (class(data[,i])=="character") data[is.na(data[,i]),i]<-'\\N'
  }
  data[is.na(data)]<-''
  
  data[,'Chr']<-factor(data[,'Chr'],levels=c(1:22,'X','Y'))
  data[,'Change_Type']<-as.factor(data[,'Change_Type'])
  data[,'Snpeff_Effect']<-as.factor(data[,'Snpeff_Effect'])
  data[,'Snpeff_Impact']<-as.factor(data[,'Snpeff_Impact'])
  data[,'SIFT_Pred']<-as.factor(data[,'SIFT_Pred'])
  data[,'PPH2_hvar_pred']<-as.factor(data[,'PPH2_hvar_pred'])
  data[,'PPH2_hdiv_pred']<-as.factor(data[,'PPH2_hdiv_pred'])
  data[,'Consensus_MAF']<-as.numeric(data[,'Consensus_MAF'])
  data[,'Consensus_MAC']<-as.numeric(data[,'Consensus_MAC'])
  data[,'CADD_Phred']<-as.numeric(data[,'CADD_Phred'])
  data[,'CADD_Raw']<-as.numeric(data[,'CADD_Raw'])
  data[,'VEST_Score']<-as.numeric(data[,'VEST_Score'])
  data[,'PPH2_hdiv_score']<-as.numeric(data[,'PPH2_hdiv_score'])
  data[,'PPH2_hvar_score']<-as.numeric(data[,'PPH2_hvar_score'])
  data[,'SIFT_Score']<-as.numeric(data[,'SIFT_Score'])
  data[,'Allelic_Depth_Ref_Prop']<-as.numeric(data[,'Allelic_Depth_Ref_Prop'])
  data[,'Allelic_Depth_Alt_Prop']<-as.numeric(data[,'Allelic_Depth_Alt_Prop'])
  data[,'Genotype_Quality']<-as.numeric(data[,'Genotype_Quality'])
  data[,'Mapping_Quality_Zero_Reads']<-as.numeric(data[,'Mapping_Quality_Zero_Reads'])
  data[,'Filters']<-as.factor(data[,'Filters'])
  
  
  filters<-getFiltersFromTable(data[,-1])
  save(file="filterVariantsSpec.Rdata",filters)
  
}


writeForImpala<-function() {
  data3<-data[1:100000,]
  write.table(file="data3.csv",data3,col.names=F,row.names=F,quote=F,sep=",")
  write.table(file="data.csv",data,col.names=F,row.names=F,quote=F,sep=",")
}

retrieveFromImpala<-function() {
  #Using impala.script4
  system.time(data<-read.csv("subset_variants.csv",sep=",",stringsAsFactors=F)) #1M, 30s
  
  condb<-dbConnect(RSQLite::SQLite(), paste0("phenotypes.db"))
  dbWriteTable(condb,"phenotypes",data,row.names=F,overwrite=T) #228MB
  dbDisconnect(condb)
  
  
}

