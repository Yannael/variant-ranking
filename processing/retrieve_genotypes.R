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


fieldsSequencingInfo<-c("patient","chr","pos","reference","alternative","zygosity","read_depth","allelic_depth_ref","allelic_depth_alt","genotype_quality",
                        "filters","downsampled","mapping_quality_zero_reads")

fieldsSequencingInfoNewNames<-c("Zygosity","Read_Depth","Allelic_Depth_Ref","Allelic_Depth_Alt","Genotype_Quality",
                                ,"Filters","Downsampled","Mapping_Quality_Zero_Reads")

fieldsVariantsInfo<-c("chr","pos","reference","alternative",
                      "change_type", "gene_symbol", "gene_ensembl", "num_genes","clinvar_rs", 
                      "consensus_MAF","consensus_MAC", "dbsnp_id_137",
                      "snpeff_effect", "snpeff_impact",
                      "cadd_phred","cadd_raw","vest_score","pph2_hdiv_score", "pph2_hdiv_pred", 
                      "pph2_hvar_score", "pph2_hvar_pred", "sift_score", "sift_pred", "short_tandem_repeat")

fieldsVariantsInfoNewNames<-c("Chr","Position","Reference","Alternative",
                              "Change_Type", "Gene_Symbol", "Gene_Ensembl", "Num_Genes","Clinvar_RS", 
                              "Consensus_MAF","Consensus_MAC", "dbsnp_id_137",
                              "Snpeff_Effect", "Snpeff_Impact",
                              "CADD_Phred","CADD_Raw","VEST_Score","PPH2_hdiv_score", "PPH2_hdiv_pred", 
                              "PPH2_hvar_score", "PPH2_hvar_pred", "SIFT_Score", "SIFT_Pred", "Short_Tandem_Repeat")

fieldsInfo<-c("patient","chr","pos","reference","alternative","zygosity","read_depth","allelic_depth_ref","allelic_depth_alt","genotype_quality",
              "filters","downsampled","mapping_quality_zero_reads",
                      "change_type", "gene_symbol", "gene_ensembl", "num_genes","clinvar_rs", 
                      "consensus_MAF","consensus_MAC", "dbsnp_id_137",
                      "snpeff_effect", "snpeff_impact",
                      "cadd_phred","cadd_raw","vest_score","pph2_hdiv_score", "pph2_hdiv_pred", 
                      "pph2_hvar_score", "pph2_hvar_pred", "sift_score", "sift_pred", "short_tandem_repeat")

fieldsInfoNewNames<-c("Patient","Chr","Position","Reference","Alternative",
                      "Zygosity","Read_Depth","Allelic_Depth_Ref","Allelic_Depth_Alt","Genotype_Quality"
                      ,"Filters","Downsampled","Mapping_Quality_Zero_Reads",
                              "Change_Type", "Gene_Symbol", "Gene_Ensembl", "Num_Genes","Clinvar_RS", 
                              "Consensus_MAF","Consensus_MAC", "dbsnp_id_137",
                              "Snpeff_Effect", "Snpeff_Impact",
                              "CADD_Phred","CADD_Raw","VEST_Score","PPH2_hdiv_score", "PPH2_hdiv_pred", 
                              "PPH2_hvar_score", "PPH2_hvar_pred", "SIFT_Score", "SIFT_Pred", "Short_Tandem_Repeat")


dummy<-function() {
  #To which DB to connect
  
  #Regular Highlander
  connectFile<-"../connectHighlander.R"
  tablename<-"exomes_ug"
  variantDBfile<-"variantsULB.db"
  sequencingDBfile<-"sequencingULB.db"
  
  #Highlander 1000 genomes
  connectFile<-"../connectHighlander2.R"
  tablename<-"1000g"
  variantDBfile<-"variants1000g.db"
  sequencingDBfile<-"sequencing1000g.db"
  
}

createCopyHighlanderDB<-function() {
  
  source(connectFile)
  
  variantsdb <- dbConnect(RSQLite::SQLite(), variantDBfile)
  sequencingdb <- dbConnect(RSQLite::SQLite(), sequencingDBfile)
  
  samples<-dbGetQuery(highlanderdb,paste0("select distinct patient from ",tablename))[,1]
  
  fields_select<-paste(c(fieldsSequencingInfo,fieldsVariantsInfo),collapse=",")
  
  for (sample in samples[1:5]) {
    
    print(paste0("Processing sample ",sample))
    
    sql<-paste0("select ",fields_select," from ",tablename," where 
                patient='", sample,"'")
    
    data <- dbGetQuery(highlanderdb,sql)
    
    #Unique ID of a variant is chr:pos:ref:alt
    uniqueID<-paste(data$chr, data$pos, data$reference, data$alternative, sep=":")
    
    zygosity<-data$zygosity
    zygosity[zygosity=="Homozygous"]<-2
    zygosity[zygosity=="Heterozygous"]<-1
    
    sequencingInfo<-cbind(sample=sample,ID=uniqueID,zygosity=zygosity,data[,fieldsSequencingInfo[-1]])
    
    #Remove dashes from sample names (cannot be used for a DB table name)
    #sample<-gsub('-','',sample)
    dbWriteTable(sequencingdb,"sequencing",sequencingInfo,append=T,row.names=F)
    
    print(paste0(nrow(sequencingInfo)," variant entries"))
    
    variantsInfo<-cbind(ID=uniqueID,data[,fieldsVariantsInfo])
    
    #if (length(unique(uniqueID))<length(uniqueID)) stop(paste0("Sample ",sample," has duplicate variant unique ID"))
    variantsInfo<-unique(variantsInfo)
    
    if (length(dbListTables(variantsdb))==0) {
      dbWriteTable(variantsdb,"variants",variantsInfo,row.names=F)
    }
    else {
      idInDB<-dbGetQuery(variantsdb,"select distinct ID from variants")[,1]
      toKeep<-setdiff(variantsInfo$ID,idInDB)
      if (length(toKeep)>0) {
        toKeepIndex<-which(variantsInfo$ID %in% toKeep)
        print(paste0(length(toKeepIndex)," variant entries added to variants DB"))
        dbWriteTable(variantsdb,"variants",variantsInfo[toKeepIndex,],append=T,row.names=F)
      }
    }
    
  }
  
  dbDisconnect(highlanderdb)
  dbDisconnect(variantsdb)
  dbDisconnect(sequencingdb)
  
}

postprocessing<-function() {
  variantDBfile<-"variants.db"
  variantsdb <- dbConnect(RSQLite::SQLite(), variantDBfile)
  data<-dbReadTable(variantsdb,"variants")
  data[is.na(data)]<-''
  colnames(data)<-c('ID',fieldsVariantsInfoNewNames)
  data<-data[order(data[,2],data[,3]),]
  dbWriteTable(variantsdb,"variants",data,overwrite=T,row.names=F)
  
  data[,'Chr']<-as.factor(data[,'Chr'])
  data[,'Change_Type']<-as.factor(data[,'Change_Type'])
  data[,'Snpeff_Effect']<-as.factor(data[,'Snpeff_Effect'])
  data[,'Snpeff_Impact']<-as.factor(data[,'Snpeff_Impact'])
  data[,'SIFT_Pred']<-as.factor(data[,'SIFT_Pred'])
  data[,'PPH2_hvar_pred']<-as.factor(data[,'PPH2_hvar_pred'])
  data[,'PPH2_hdiv_pred']<-as.factor(data[,'PPH2_hdiv_pred'])
  
  filtersVariantsTypes<-getFiltersFromTable(data[,-1])
  save(file="filtersVariantsTypes.Rdata",filtersVariantsTypes)
  
  variantsdb <- dbConnect(RSQLite::SQLite(), variantDBfile)
  sequencingdb <- dbConnect(RSQLite::SQLite(), sequencingDBfile)
  
  samples<-dbListTables(sequencingdb)
  write.table(file="samples.txt",samples,col=F,row=T,quote=F)
  
}

createCopyHighlanderDB2<-function() {
  
  source(connectFile)
  
  variantsdb <- dbConnect(RSQLite::SQLite(), variantDBfile)
  sequencingdb <- dbConnect(RSQLite::SQLite(), sequencingDBfile)
  
  #samples<-dbGetQuery(highlanderdb,paste0("select distinct patient from ",tablename))[,1]
  fields_select<-paste(unique(c(fieldsSequencingInfo,fieldsVariantsInfo)),collapse=",")
  
  #Takes about 3 minutes
  system.time(data <- dbGetQuery(highlanderdb,paste0("select ",fields_select," from ",tablename," limit 5000000")))
  uniqueID<-paste(data$chr, data$pos, data$reference, data$alternative, sep=":")
  
  id.redundant<-tapply(data$gene_symbol,uniqueID,length)
  
  zygosity<-data$zygosity
  zygosity[zygosity=="Homozygous"]<-2
  zygosity[zygosity=="Heterozygous"]<-1
  data$zygosity<-zygosity
  data<-cbind(ID=uniqueID,data)
  
  data[is.na(data)]<-''
  colnames(data)<-c('ID',fieldsInfoNewNames)
  data<-data[order(data[,3],data[,4]),]
  
  data[,'Chr']<-as.factor(data[,'Chr'])
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
  
  dbWriteTable(variantsdb,"variants",data,overwrite=T,row.names=F)
  
  filtersVariantsTypes<-getFiltersFromTable(data[,-1])
  save(file="filtersVariantsTypes.Rdata",filtersVariantsTypes)
  
}