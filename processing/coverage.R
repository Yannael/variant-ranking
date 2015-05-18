library(RMySQL)

createCoverageDB<-function() {
  
  
  con <- dbConnect(RSQLite::SQLite(), "../../coverage.db")
  
  #Add X
  for (num in 1:22) {
    if (num<10) add0<-"0" else add0<-""
    
    filename<-paste0("/home/yleborgn/bridge/data/erasme/coverage/ISDBMbatch_coverage_chr",add0,num,"_trios")
    
    #Takes about 2 minutes
    system.time(CoverageTab<-read.table(filename,header=T))
    
    #30s
    system.time({
      dbWriteTable(con,paste0("chr",num),CoverageTab,overwrite=T)
    })
  }
  
  dbDisconnect(con)
  
}

getMappingZH_ISDBM<-function() {
  
  #Change control or patho
  sql<-"select * from patho, chr17
  where patho.Locus=chr17.Locus"
  
  system.time({
    rs = dbSendQuery(concoverage,sql )
    controlSNPs = fetch(rs, n=-1)
  })
  
  ids<-unique(controlSNPs[,1])
  
  ISid<-NULL
  for (i in 1:length(ids)) {
    ss<-controlSNPs[which(controlSNPs[,1]==ids[i]),c(4,22:54)]
    name<-colnames(ss[1:15,c(1,which.max(cor(ss[,1],ss[,2:ncol(ss)]))+1)])[2]
    ISid<-c(ISid,strsplit(name,"_for_")[[1]][2])
  }

}

retrieveRefEntriesChr<-function(chr,patho,control) {
  
  #2s
  system.time({
    rs<-dbSendQuery(concoverage,paste("select Average_Depth_sample from ",chr))
    avg_depth = fetch(rs, n=-1)
  })
  thresh<-median(avg_depth[,1])
  
  #30s
  system.time({
    rs<-dbSendQuery(concoverage,paste0("select * from locus_snps,",chr," where locus_snps.Locus=",chr,".Locus"))
    coverage<-fetch(rs, n=-1)
  })
  
  snpsMat<-coverage[,5:37]
  ISDBMid<-as.vector(unlist(sapply(colnames(snpsMat),substr,11,21)))
  mapping<-read.table('mappingZH_ISDBM.txt',stringsAsFactors=F)
  i<-match(ISDBMid,mapping[,2])
  colnames(snpsMat)<-mapping[i,1]
  rownames(snpsMat)<-coverage$Locus
  
  snpsMat[snpsMat<thresh]<-NA
  snpsMat[snpsMat>=thresh]<-0
  
  select.entry.homozygous<-which((patho$Locus %in% coverage$Locus) & (patho$zygosity=="Homozygous"))
  for (i in select.entry.homozygous) 
    snpsMat[patho$Locus[i],patho$patient[i]]<-2
  
  select.entry.heterozygous<-which((patho$Locus %in% coverage$Locus) & (patho$zygosity=="Heterozygous"))
  for (i in select.entry.heterozygous) 
    snpsMat[patho$Locus[i],patho$patient[i]]<-1
  
  snpsMat
}

retrieveRefEntries<-function() {
  
  #Move patho and control groups to coverage DB
  congroups <- dbConnect(RSQLite::SQLite(), "../../groupsToComparePartial.db")
  concoverage<- dbConnect(RSQLite::SQLite(), "../../coverage.db")
  
  rs<-dbSendQuery(congroups,"select * from patho" )
  patho<-fetch(rs, n=-1)
  rs<-dbSendQuery(congroups,"select * from control" )
  control<-fetch(rs, n=-1)
  
  #dbWriteTable(concoverage,"patho",patho,overwrite=T)
  #dbWriteTable(concoverage,"control",control,overwrite=T)
  
  #getMappingZH_ISDBM
  
  locusSNPs<-data.frame(unique(patho$Locus),stringsAsFactors=F)
  colnames(locusSNPs)<-"Locus"
  dbWriteTable(concoverage,"locus_snps",locusSNPs,overwrite=T)
  
  snpsMat<-NULL
  
  for (chr in 1:21) {
    print(chr)
    st<-system.time(snpsMat<-rbind(snpsMat,retrieveRefEntriesChr(paste0("chr",chr),patho,control)))
    print(st)
                   
  }
  
  rs<-dbSendQuery(congroups,"select * from dbSNPs" )
  dbSNPs<-fetch(rs, n=-1)
  
  mapping_gene<-match(rownames(snpsMat),dbSNPs$Locus)
  rownames_extended<-paste(dbSNPs$gene_ensembl[mapping_gene],rownames(snpsMat),sep="_")
  rownames(snpsMat)<-rownames_extended
  snpsMat<-snpsMat[sort(rownames(snpsMat),index.return=T)$ix,]
  
  patient_patho<-unique(patho$patient)
  match.patient<-match(patient_patho,colnames(snpsMat))
  patient_patho<-paste(patient_patho,"Patho",sep=":")
  colnames(snpsMat)[match.patient]<-patient_patho
  
  patient_control<-unique(control$patient)
  match.patient<-match(patient_control,colnames(snpsMat))
  patient_control<-paste(patient_control,"Control",sep=":")
  colnames(snpsMat)[match.patient]<-patient_control
  
  save(file="../../snpsMat.Rdata",snpsMat)
  write.table(file="../../snpsMat.txt",snpsMat)
  
}
