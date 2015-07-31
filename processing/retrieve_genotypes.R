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


fieldsSequencingInfo<-c("zygosity","read_depth","allelic_depth_ref","allelic_depth_alt")
fieldsVariantsInfo<-c("chr","pos","reference","alternative",
                      "change_type", "gene_symbol", "gene_ensembl", "num_genes","clinvar_rs", 
                      "consensus_MAF","consensus_MAC", "dbsnp_id_137",
                      "snpeff_effect", "snpeff_impact",
                      "cadd_phred","cadd_raw","vest_score","pph2_hdiv_score", "pph2_hdiv_pred", 
                      "pph2_hvar_score", "pph2_hvar_pred", "sift_score", "sift_pred", "short_tandem_repeat")


dummy<-function() {
  #To which DB to connect
  
  #Regular Highlander
  connectFile<-"../connectHighlander.R"
  tablename<-"exomes_ug"
  variantDBfile<-"variants.db"
  sequencingDBfile<-"sequencing.db"
  
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
  
  for (sample in samples) {
    
    print(paste0("Processing sample ",sample))
    
    sql<-paste0("select ",fields_select," from ",tablename," where 
                patient='", sample,"'")
    
    data <- dbGetQuery(highlanderdb,sql)
    
    #Unique ID of a variant is chr:pos:ref:alt
    uniqueID<-paste(data$chr, data$pos, data$reference, data$alternative, sep=":")
    
    zygosity<-data$zygosity
    zygosity[zygosity=="Homozygous"]<-2
    zygosity[zygosity=="Heterozygous"]<-1
    
    sequencingInfo<-cbind(ID=uniqueID,zygosity=zygosity,read_depth=data$read_depth,
                          allelic_depth_ref=data$allelic_depth_ref,
                          allelic_depth_alt=data$allelic_depth_alt)
    
    #Remove dashes from sample names (cannot be used for a DB table name)
    sample<-gsub('-','',sample)
    dbWriteTable(sequencingdb,sample,as.data.frame(sequencingInfo),overwrite=T,row.names=F)
    
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

dummy<-function() {
  
  variantsdb <- dbConnect(RSQLite::SQLite(), variantDBfile)
  sequencingdb <- dbConnect(RSQLite::SQLite(), sequencingDBfile)
  
  samples<-dbListTables(sequencingdb)
  write.table(file="samples.txt",samples,col=F,row=T,quote=F)
  
}


