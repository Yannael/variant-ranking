library(RMySQL)

createCoverageDB<-function() {
  
  
  con <- dbConnect(RSQLite::SQLite(), "coverage.db")
  
  #Add X
  for (num in 16:22) {
    if (num<10) add0<-"0" else add0<-""
    
    filename<-paste0("/home/yleborgn/bridge/data/erasme/coverage/ISDBMbatch_coverage_chr",add0,num,"_trios_seqCap")
    
    #Takes about 5 minutes
    st<-system.time(CoverageTab<-read.table(filename,header=T))
    print(st)
    #30s
    system.time({
      dbWriteTable(con,paste0("chr",num),CoverageTab,overwrite=T)
    })
  }
  for (num in c("X","Y")) {
    add0<-""
    
    filename<-paste0("/home/yleborgn/bridge/data/erasme/coverage/ISDBMbatch_coverage_chr",add0,num,"_trios_seqCap")
    
    #Takes about 5 minutes
    st<-system.time(CoverageTab<-read.table(filename,header=T))
    print(st)
    #30s
    system.time({
      dbWriteTable(con,paste0("chr",num),CoverageTab,overwrite=T)
    })
  }
  
  dbDisconnect(con)
  
}


countNA<-function(x) {length(which(is.na(x)))}

retrieveRefEntriesChr<-function(chr,patho,control) {
  
  thresh<-10
  
  #1,5min
  system.time({
    rs<-dbSendQuery(concoverage,paste0("select * from locus,",chr," where locus.Locus=",chr,".Locus"))
    coverage<-fetch(rs, n=-1)
  })
  
  varMat<-coverage[,7:39]
  ISDBMid<-as.vector(unlist(sapply(colnames(varMat),substr,11,21)))
  mapping<-read.table('mappingZH_ISDBM.txt',stringsAsFactors=F)
  i<-match(mapping[,2],ISDBMid)
  colnames(varMat)[i]<-mapping[,1]
  varMat<-varMat[,i]
  
  rownames(varMat)<-apply(coverage[,1:3],1,paste,collapse=":")
  
  varMat[varMat<thresh]<-NA
  varMat[varMat>=thresh]<-0
  
  nb.NA<-apply(varMat,1,countNA)
  i.remove<-which(nb.NA>0)
  varMat<-varMat[-i.remove,]
  
  pathoVars<-paste(patho$Locus,patho$reference,patho$alternative,sep=":")
  controlVars<-paste(control$Locus,control$reference,control$alternative,sep=":")
  idRow<-rownames(varMat)
  
  select.entry.homozygous<-which((pathoVars %in% idRow) & (patho$zygosity=="Homozygous"))
  for (i in select.entry.homozygous) 
    varMat[pathoVars[i],patho$patient[i]]<-2
  
  select.entry.heterozygous<-which((pathoVars %in% idRow) & (patho$zygosity=="Heterozygous"))
  for (i in select.entry.heterozygous) 
    varMat[pathoVars[i],patho$patient[i]]<-1
  
  select.entry.homozygous<-which((controlVars %in% idRow) & (control$zygosity=="Homozygous"))
  for (i in select.entry.homozygous) 
    varMat[controlVars[i],control$patient[i]]<-2
  
  select.entry.heterozygous<-which((controlVars %in% idRow) & (control$zygosity=="Heterozygous"))
  for (i in select.entry.heterozygous) 
    varMat[controlVars[i],control$patient[i]]<-1
  
  varMat<-cbind(Locus=coverage[-i.remove,1],reference=coverage[-i.remove,2],alternative=coverage[-i.remove,3],varMat)
  varMat
}

retrieveRefEntries<-function() {
  
  #Move patho and control groups to coverage DB
  congroups <- dbConnect(RSQLite::SQLite(), "groupsToComparePartial.db")
  concoverage<- dbConnect(RSQLite::SQLite(), "coverage.db")
  
  rs<-dbSendQuery(congroups,"select * from patho" )
  patho<-fetch(rs, n=-1)
  rs<-dbSendQuery(congroups,"select * from control" )
  control<-fetch(rs, n=-1)
  
  #dbWriteTable(concoverage,"patho",patho,overwrite=T)
  #dbWriteTable(concoverage,"control",control,overwrite=T)
  
  locusvariants<-data.frame(unique(patho[,c("Locus",'reference','alternative')]),stringsAsFactors=F)
  dbWriteTable(concoverage,"locus",locusvariants,overwrite=T)
  
  varMat<-NULL
  
  for (chr in c("X","Y")) {
    print(chr)
    st<-system.time(varMat<-rbind(varMat,retrieveRefEntriesChr(paste0("chr",chr),patho,control)))
    print(st)
                   
  }
  
  rs<-dbSendQuery(congroups,"select * from dbvariants" )
  dbvariants<-fetch(rs, n=-1)
  
  
  mapping_gene<-match(varMat$Locus,dbvariants$Locus)
  varMat<-cbind(geneid=dbvariants$gene_ensembl[mapping_gene],varMat)
  varMat<-varMat[sort(as.character(varMat$Locus),index.return=T)$ix,]
  
  #i.keep<-which((dbvariants[,'consensus_MAF']<0.01 | dbvariants[,'consensus_MAC']<1000) &(dbvariants[,'gene_symbol']!="NA"))
  i.keep<-which((dbvariants[,'consensus_MAF']<0.01 | dbvariants[,'consensus_MAC']<1000) )
  id.keep<-dbvariants[i.keep,'uniqueid']
  i.keep<-match(id.keep,rownames(varMat))
  varMat<-varMat[i.keep,]
  
  write.table(file="varMat2.txt",varMat,col.names=F,row.names=F,quote=F)
  
}

dummy<-function() {
  
  #In highlander: Select patient, read_depth, patho, chr, pos columns
  #Filter with patient ZH135614, chr1
  #Sort by ascending order
  #
  #Note discrepancy in results. (e.g., locus 1:14653, or 1:17538)
  dbWriteTable(concoverage,"patho",patho)
  
  
  #Get chr1 info dor child in trio 1
  rs<-dbSendQuery(concoverage,"select Locus, Average_Depth_sample, Depth_for_ISDBM387751 from chr1")
  datacoverage<-fetch(rs, n=-1)
  
  rs<-dbSendQuery(concoverage,"select * from patho where patient='ZH135614' and chr=1")
  datapatho<-fetch(rs, n=-1)
  
  rs<-dbSendQuery(concoverage,"select patho.Locus, chr, pos, read_depth, Depth_for_ISDBM387751 from patho,chr1 
                  where patient='ZH135614' and chr=1 and patho.Locus=chr1.Locus")
  datamatch<-fetch(rs, n=-1)
  
  
  datamatch<-datamatch[sort(datamatch$pos,index.ret=T)$ix,]
  datamatch[1:10,]
  
  i.diff<-which(abs(datamatch[,4]-datamatch[,5])>10)
  length(i.diff)
  datamatch[i.diff[1:5],]
  
  #Also, in Highlander, 
  #select patient chr, pos, num_genes, gene_ensembl, db_snp_id_137, filters
  # for chr4, pos 54969778
  #To see inconsistent info
  
  #3) Why is there no missing value in coverage file (i.e., at one locus, there should be missing data for some samples) 
  
}
