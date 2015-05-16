library(RMySQL)

createCoverageDB<-function() {
  
  
  con <- dbConnect(RSQLite::SQLite(), "../../coverage.db")
  
  #Add X
  for (num in 1:22) {
    if (num<10) add0<-"0" else add0<-""
    
    filename<-paste0("/home/yleborgn/bridge/data/erasme/coverage/ISDBMbatch_coverage_chr",add0,num,"_trios")
    
    #Takes about 2 minutes
    system.time(CoverageTab<-read.table(filename,header=T))
    
    #Retrieve chr and pos data in two different columns
    #Maybe not the most efficient way, takes about 3 minutes
    #system.time({
    #  chr_pos<-matrix(unlist(apply(CoverageTab[,1,drop=F],1,strsplit,":")),nrow(CoverageTab),2,byrow=T)
    #})
    
    #10s
    #system.time({
    #  CoverageTab<-cbind(chr=chr_pos[,1],pos=chr_pos[,2],CoverageTab)
    #})
    
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
  
  #2s
  system.time({
    rs<-dbSendQuery(con,"select Average_Depth_sample from chr17")
    avg_depth = fetch(rs, n=-1)
  })
  thresh<-median(avg_depth[,1])
  
  #30s
  system.time({
  rs<-dbSendQuery(concoverage,"select * from locus_snps,chr17 where locus_snps.Locus=chr17.Locus")
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
  snpsMat[patho$Locus[select.entry.homozygous],patho$patient[select.entry.homozygous]]<-2
  
  select.entry.heterozygous<-which((patho$Locus %in% coverage$Locus) & (patho$zygosity=="Heterozygous"))
  snpsMat[patho$Locus[select.entry.heterozygous],patho$patient[select.entry.heterozygous]]<-1
  
}

dummy<-function() {
  system.time({
    rs<-dbSendQuery(con,"select * from chr17")
    data = fetch(rs, n=-1)
  })
  
  #Locus Total_Depth Average_Depth_sample Depth_for_ISDBM379928 Depth_for_ISDBM379929
  #1 17:6011        1171                35.48                    28                    25
  #Depth_for_ISDBM379930 Depth_for_ISDBM379931 Depth_for_ISDBM379932 Depth_for_ISDBM379933
  
  system.time({
    rs<-dbSendQuery(con,"select Average_Depth_sample from chr17")
    avg_depth = fetch(rs, n=-1)
  })
  
  dbClearResult(rs)
  
  
}