library(RMySQL)

createCoverageDB<-function() {
  
  
  con <- dbConnect(RSQLite::SQLite(), "../../coverage2.db")
  
  #Add X
  for (num in 1:1) {
    if (num<10) add0<-"0" else add0<-""
    
    filename<-paste0("/home/yleborgn/bridge/data/erasme/coverage2/ISDBMbatch_coverage_chr",add0,num,"_trios")
    
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
  sql<-"select * from control, chr17
  where control.Locus=chr17.Locus"
  
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

countNA<-function(x) {length(which(is.na(x)))}

retrieveRefEntriesChr<-function(chr,patho,control) {
  
  #2s
  #system.time({
  #  rs<-dbSendQuery(concoverage,paste("select Average_Depth_sample from ",chr))
  #  avg_depth = fetch(rs, n=-1)
  #})
  #thresh<-median(avg_depth[,1])
  thresh<-15
  
  #10min
  system.time({
    rs<-dbSendQuery(concoverage,paste0("select * from locus_snps,",chr," where locus_snps.Locus=",chr,".Locus"))
    coverage<-fetch(rs, n=-1)
  })
  
  snpsMat<-coverage[,6:38]
  ISDBMid<-as.vector(unlist(sapply(colnames(snpsMat),substr,11,21)))
  mapping<-read.table('mappingZH_ISDBM.txt',stringsAsFactors=F)
  i<-match(ISDBMid,mapping[,2])
  colnames(snpsMat)<-mapping[i,1]
  rownames(snpsMat)<-coverage$Locus
  
  snpsMat[snpsMat<thresh]<-NA
  snpsMat[snpsMat>=thresh]<-0
  
  nb.NA<-apply(snpsMat,1,countNA)
  i.remove<-which(nb.NA>0)
  snpsMat<-snpsMat[-i.remove,]
  
  select.entry.homozygous<-which((patho$Locus %in% rownames(snpsMat)) & (patho$zygosity=="Homozygous"))
  for (i in select.entry.homozygous) 
    snpsMat[patho$Locus[i],patho$patient[i]]<-2
  
  select.entry.heterozygous<-which((patho$Locus %in% rownames(snpsMat)) & (patho$zygosity=="Heterozygous"))
  for (i in select.entry.heterozygous) 
    snpsMat[patho$Locus[i],patho$patient[i]]<-1
  
  select.entry.homozygous<-which((control$Locus %in% rownames(snpsMat)) & (control$zygosity=="Homozygous"))
  for (i in select.entry.homozygous) 
    snpsMat[control$Locus[i],control$patient[i]]<-2
  
  select.entry.heterozygous<-which((control$Locus %in% rownames(snpsMat)) & (control$zygosity=="Heterozygous"))
  for (i in select.entry.heterozygous) 
    snpsMat[control$Locus[i],control$patient[i]]<-1
  
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
  
  for (chr in 1:1) {
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
  
  patient_data<-read.table("mappingZH_ISDBM_withFMC.txt")
  
  match.patient<-match(patient_data[,1],colnames(snpsMat))
  colnames(snpsMat)[match.patient]<-paste(patient_data[,1],patient_data[,3],sep="_")
  
  save(file="../../snpsMat.Rdata",snpsMat)
  write.table(file="../../snpsMat.txt",snpsMat)
  
}

dummy<-function() {
  
  #In highlander: Select patient, read_depth, patho, chr, pos columns
  #Filter with patient ZH135614, chr1
  #Sort by ascending order
  #
  #Note discrepancy in results. (e.g., locus 1:14653, or 1:17538)
  
  #Get chr1 info dor child in trio 1
  rs<-dbSendQuery(concoverage,"select Locus, Average_Depth_sample, Depth_for_ISDBM387751 from chr1")
  datacoverage<-fetch(rs, n=-1)
  
  rs<-dbSendQuery(concoverage,"select * from patho where patient='ZH135614' and chr=1")
  datapatho<-fetch(rs, n=-1)
  
  rs<-dbSendQuery(concoverage,"select patho.Locus, chr, pos, read_depth, Depth_for_ISDBM387751 from patho,chr1 
                  where patient='ZH135614' and chr=1 and patho.Locus=chr1.Locus")
  datamatch<-fetch(rs, n=-1)
  
  datamatch[1:10,]
  
  
  datamatch<-datamatch[sort(datamatch$pos,index.ret=T)$ix,]
  
  #Also, in Highlander, 
  #select patient chr, pos, num_genes, gene_ensembl, db_snp_id_137, filters
  # for chr4, pos 54969778
  #To see inconsistent info
  
  #3) Why is there no missing value in coverage file (i.e., at one locus, there should be missing data for some samples) 
  
}
