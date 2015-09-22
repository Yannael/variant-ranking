require(RMySQL)

priorlist<-c("TUBB2A","TCL1A","BRDG1","HTPAP","PPPAPDC1B","IGKV4","IGLL1",
             "IGKV1D-13","CD798","HS3ST1","SH2D1B","MS4A1","TLR5","FCRL1",
             "PNOC","SLC8A1","FCRL2",
             "FAM26F","POU4F2", "SIX4", "ADAMTS7", "GRIN2C","CHPF","HRH3",
             "ANK3","MUC17","CCDC125","TP73","PHLDA2","RAB6C",
             "FCGBP","OBSCN","UNC79","WASL","FAM26F","DDX20","NUDT16",
             "IQSEC2","LOC1002407","CD40L","POP7",
             "IQSEC2","CD40L","POP7","TDG","PRR23C","MRPS24","NARR","GCGR",
             "POU4F2","MUC17","ADAMTS7","INTS1","KRT6C",
             "PKHD1L1","SDHA","MYH3","ETFA","MYO3A","PTPRN2","CHPPF","YPEL5","HRH3","PCIF1")

geneIDs<-c("TUBB2A","TCL1A")

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
              "consensus_MAF","consensus_MAC", "dbsnp_id_137",
              "snpeff_effect", "snpeff_impact",
              "cadd_phred","cadd_raw","vest_score","pph2_hdiv_score", "pph2_hdiv_pred", 
              "pph2_hvar_score", "pph2_hvar_pred", "sift_score", "sift_pred", "short_tandem_repeat")

fieldsInfoNewNames<-c("Sample_ID","Chr","Position","Reference","Alternative",
                      "Zygosity","Read_Depth","Allelic_Depth_Ref_Prop","Allelic_Depth_Alt_Prop","Genotype_Quality"
                      ,"Filters","Downsampled","Mapping_Quality_Zero_Reads","Allele_Num",
                      "Change_Type", "Gene_Symbol", "Gene_Ensembl","Transcript_Ensembl", "Num_Genes","Clinvar_RS", 
                      "Consensus_MAF","Consensus_MAC", "dbsnp_id_137",
                      "Snpeff_Effect", "Snpeff_Impact",
                      "CADD_Phred","CADD_Raw","VEST_Score","PPH2_hdiv_score", "PPH2_hdiv_pred", 
                      "PPH2_hvar_score", "PPH2_hvar_pred", "SIFT_Score", "SIFT_Pred", "Short_Tandem_Repeat")



dummy<-function() {
  #To which DB to connect
  
  #Regular Highlander
  connectFile<-"../connectHighlander.R"
  tablename<-"exomes_hc"
  variantDBfile<-"variants.db"
  sequencingDBfile<-"sequencingULB.db"
  
  #Highlander 1000 genomes
  connectFile<-"../connectHighlander2.R"
  tablename<-"1000g"
  variantDBfile<-"variants1000g.db"
  sequencingDBfile<-"sequencing1000g.db"
  
}

createIndices<-function(dbname,tablename,cnames) {
  for (cname in cnames) {
    print(cname)
    dbGetQuery(dbname,paste("create index ",cname," on ",tablename,"(",cname,")"))
  }
}

createCopyHighlanderDB<-function() {
  
  source(connectFile)
  
  samples<-dbGetQuery(highlanderdb,paste0("select distinct patient from ",tablename))[,1]
  fields_select<-paste(unique(c(fieldsInfo)),collapse=",")
  
  #For 1000g
  subsamples<-paste0("('",paste(samples[1:3],sep="",collapse="','"),"')")
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
  
  data[is.na(data)]<-''
  colnames(data)<-c('ID',fieldsInfoNewNames)
  #data<-data[order(data[,3],data[,4]),]
  
  variantsdb <- dbConnect(RSQLite::SQLite(), "variants.db")
  dbWriteTable(variantsdb,"variants",data,overwrite=T,row.names=F)
  
  phenotypesdb<-dbConnect(RSQLite::SQLite(), paste0("phenotypes.db"))
  phenotypes<-dbReadTable(phenotypesdb,"phenotypes")
  dbWriteTable(variantsdb,"phenotypes",phenotypes,overwrite=T,row.names=F)
  
  #2 minutes with 15M entries
  system.time(data<-dbGetQuery(variantsdb,paste0("select * from variants,phenotypes where variants.Sample_ID=phenotypes.Sample_ID")))
  
  #Remove duplicate Sample_ID column
  data<-data[,-which(colnames(data)=="Sample_ID")[2]]
  
  data[data==""]<-NA
  
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
  
  datadb<-dbConnect(RSQLite::SQLite(), paste0("data.db"))
  dbWriteTable(datadb,"variants",data,row.names=F,overwrite=T)
  dbDisconnect(datadb)
  
  missingGenes<-datawhich(data[,'Gene_Ensembl']=='')
  ID<-data[missingGenes,'ID']
  subData<-data[which(data[,'ID'] %in% ID),c('ID','Gene_Ensembl')]
  tt<-tapply(subData[,1],subData[,1],length)
  
}

countLevels<-function(vec) {
  length(unique(vec))
}

dummy2<-function() {
  
  system.time(data<-dbGetQuery(data2db,paste0("select * from variants where Data_Source='1000 Genomes' limit 100000")))
  system.time(data<-dbGetQuery(datadb,paste0("select count(*) from variants,phenotypes where variants.Sample_ID=phenotypes.Sample_ID")))
  
}
