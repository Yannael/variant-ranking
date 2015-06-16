Sys.setenv(HADOOP_CMD="/usr/bin/hadoop",
           HADOOP_STREAMING="/opt/cloudera/parcels/CDH-5.1.3-1.cdh5.1.3.p0.12/lib/hadoop-0.20-mapreduce/contrib/streaming/hadoop-streaming-2.3.0-mr1-cdh5.1.3.jar")

library('rhdfs')
library('rmr2')
hdfs.init()

infile<-"in.variantMat.csv2"
#infile<-"snpsMat.tab"
outfile<-paste0(gsub(" ","_",Sys.time()),"_out.variantMat.csv")

local<-FALSE
if (local){
  rmr.options(backend = "local")
  hdfsfile<-paste0("/home/yleborgn/bridge/",infile)
  hdfsoutfile<-paste0("/home/yleborgn/bridge/",outfile)
} else {
  rmr.options(backend = "hadoop")
  bp = rmr.options("backend.parameters")
  #bp$hadoop<-c(list(D="mapreduce.map.memory.mb=8192"),bp$hadoop, list( D = "mapred.task.timeout=3600000"), list(D="mapreduce.map.memory.mb=8192"))
  bp$hadoop<-c(bp$hadoop, list( D = "mapred.task.timeout=3600000"))
  bp = rmr.options(backend.parameters=bp);
  
  print(rmr.options("backend.parameters"))
  
  hdfsfile<-paste0("/user/yleborgn/bridge/",infile)
  hdfsoutfile<-paste0("/user/yleborgn/bridge/",outfile)
}

createData<-T
if (createData) {
  prof.init<-system.time({
    load("../../snpsMat.Rdata")
    
    #alldata.i<-which(!is.na(apply(snpsMat,1,sum)))
    #snpsMat<-snpsMat[alldata.i,]
    
    #snpsMat<-snpsMat[1:10000,]
    
    snpsMat<-cbind(id=1:nrow(snpsMat), snpsMat)
    prof.savehdfs<-system.time(to.dfs(snpsMat,output=hdfsfile))
    
  })
}

Bcount <-  function(strng){
  chars = strsplit(strng,"")[[1]]
  length(which(chars == "B" ))
}

reduce =   function(k, v) { 
  #if (k=="ATP1A1") browser()
  tab<-v[[1]]$tab
  cnt<-v[[1]]$cnt
  gene.rs<-v[[1]]$gene.rs
  if (length(v)>1){
    for (i in 2:length(v)){
      tab<-tab+v[[i]]$tab
      cnt<-cbind(cnt,v[[i]]$cnt)
      gene.rs<-c(gene.rs,v[[i]]$gene.rs)
    }
  }
  
  pv<-1
  st<-10
  
  #if (sum(tab["case",])>(length(cases)/2) & sum(tab["ctrl",])>0){    
  #  ttest<- prop.trend.test(tab["case",], margin.table(tab,2))
  #  pv<-ttest$p.value ##pvalue(IT) ##c(L.case,L.control))
  #}
  pv<-0
  
  st<-entropy.tab(tab) #CHANGE was info.tab
  delta<--10
  if (all(apply(tab,1,sum)>0)){
    tab2<-tab/apply(tab,1,sum)
    delta<-sum(tab2["case",2:3]-tab2["ctrl",2:3])
  }
  keyval(k,list(list(tab=tab,pv=pv,st=st,delta=delta,gene.rs=gene.rs,cnt=cnt)))
}

map = function(k, X) {       
  listX<-NULL
  val<-NULL
  
  #Rows of X are variants. Row name is structure geneid_keysplit_rsid (e.g., OR4F5_keysplit_rs200676709)
  #We split around the token '_keysplit_' to retrieve geneid and rsid. 
  geneid_rsid<-rownames(X)
  geneid_rsid<-sapply(geneid_rsid,strsplit,'_')
  geneid_rsid<-unlist(geneid_rsid,use.name=F)
  geneid_rsid<-matrix(geneid_rsid,nrow(X),2,byrow=T)
  
  N<-nrow(X)
  browser()
  patient_class<-colnames(X)
  patient_class<-sapply(patient_class,strsplit,'\\.')
  patient_class<-unlist(patient_class,use.name=F)
  patient_class<-matrix(patient_class,ncol(X),2,byrow=T)
  
  controls<-which(patient_class[,2]=="Control")
  cases<-which(patient_class[,2]=="Patho")
  
  for (gene in unique(geneid_rsid[,1])) {#length(gene2rs)){
    
    gene.rs.i<-which(geneid_rsid[,1]==gene)
    tab<-array(0,c(2,3))
    colnames(tab)<-c(0,1,2)
    rownames(tab)<-c("ctrl","case")
    
    cnt<-NULL 
    
    if (length(gene.rs.i)>0 ){       
      for (n in 1:ncol(tab)){
        
        tab["ctrl",n]<-length(which(X[gene.rs.i,controls]==(n-1)))          
        tab["case",n]<-length(which(X[gene.rs.i,cases]==(n-1)))
      }
      
      listX<-c(listX,list(list(tab=tab,gene.rs=geneid_rsid[gene.rs.i,2])))
      val<-c(val,gene)
      
    }
  }
  
  keyval(val,listX)
}

countGeno<-function(x) {
  
  nb0<-sum(as.numeric(x==0))
  nb1<-sum(as.numeric(x==1))
  nb2<-sum(as.numeric(x==2))
  nbNA<-sum(as.numeric(is.na(x)))
  c(nb0,nb1,nb2,nbNA)
  
}

map.denovo = function(k, X) {       
  listX<-NULL
  val<-NULL
  N<-nrow(X)
  
  loci<-rownames(X)
  
  for (locus in 1:nrow(X)) {
    
    listX<-c(listX,list(list(snpsMat=X[locus,])))
    val<-c(val,rownames(X[locus,]))
    
  }
  
  keyval(val,listX)
}

reduce.denovo =   function(k, v) { 
  #browser()
  v<-v[[1]]$snpsMat[,2:34]
  patient_data<-colnames(v)
  patient_data<-sapply(patient_data,strsplit,'_')
  patient_data<-unlist(patient_data,use.name=F)
  patient_data<-matrix(patient_data,ncol(v),3,byrow=T)
  
  controls<-which(patient_data[,2]!="C")
  cases<-which(patient_data[,2]=="C")
  
  tab<-array(0,c(2,4))
  tab[1,]<-countGeno(v[,controls])
  tab[2,]<-countGeno(v[,cases])
  
  rownames(tab)<-c("ctrl","case")
  colnames(tab)<-c(0,1,2,"NA")
  
  score<-0
  for (i in 1:11) {
    genoC<-v[,which(patient_data[,2]=="C" & patient_data[,3]==i)]
    genoF<-v[,which(patient_data[,2]=="F" & patient_data[,3]==i)]
    genoM<-v[,which(patient_data[,2]=="M" & patient_data[,3]==i)]
    
    score<-score+as.numeric(genoC==1 & genoM==0 & genoF==0)
  }
  
  keyval(k,list(list(tab=tab,score=score)))
}

lpart<-1000
NAllpart<-6

map.pairs = function(k, X) {       
  val<-NULL
  
  N<-nrow(X)
  #browser()
  #In how many partitions must be split rows of X
  npart<-ceiling(N/lpart)+1
  #What is the starting partition in X index
  startPart<-floor(X[1,1]/lpart)
  
  listX<-NULL
  
  for (i in (startPart:(startPart+npart-1))) {
    indexpart<-which((X[,1]>=(i*lpart+1)) & (X[,1]<=((i+1)*lpart)))
    if (length(indexpart)>0) {
    listX<-c(listX,rep(list(list(snpsMat=X[indexpart,])),NAllpart))
    namePart<-NULL
    if (i>0) {
      namePart<-paste(0:(i-1),"_",i,sep="")
    }
    namePart<-c(namePart,paste(i,"_",i:(NAllpart-1),sep=""))
    val<-c(val,namePart)
    }
  }
  #browser()
  
  keyval(val,listX)
}

reduce.pairs =   function(k, v) { 
  
  browser()
  
#   v<-v[[1]]$snpsMat
#   patient_data<-colnames(v)
#   patient_data<-sapply(patient_data,strsplit,'_')
#   patient_data<-unlist(patient_data,use.name=F)
#   patient_data<-matrix(patient_data,ncol(v),3,byrow=T)
#   
#   controls<-which(patient_data[,2]!="C")
#   cases<-which(patient_data[,2]=="C")
#   
#   tab<-array(0,c(2,4))
#   tab[1,]<-countGeno(v[,controls])
#   tab[2,]<-countGeno(v[,cases])
#   
#   rownames(tab)<-c("ctrl","case")
#   colnames(tab)<-c(0,1,2,"NA")
#   
#   score<-0
#   for (i in 1:11) {
#     genoC<-v[,which(patient_data[,2]=="C" & patient_data[,3]==i)]
#     genoF<-v[,which(patient_data[,2]=="F" & patient_data[,3]==i)]
#     genoM<-v[,which(patient_data[,2]=="M" & patient_data[,3]==i)]
#     
#     score<-score+as.numeric(genoC==1 & genoM==0 & genoF==0)
#   }
  
  keyval(k,list(list(k)))
}



dummy<-function() {
  
  prof.run<-system.time({
    PVmr = mapreduce(
      input=hdfsfile,  
      # input.format=dat.in.format,
      map = map.denovo,
      reduce=reduce.denovo
    )
  })
  
  
  
  prof.run<-system.time({
    PVmr = mapreduce(
      input=hdfsfile,  
      # input.format=dat.in.format,
      map = map.pairs,
      reduce=reduce.pairs
    )
  })
  
  
  results<-from.dfs(PVmr)
  #to.dfs(results,output='/user/yleborgn/bridge/out16000.results')
  #to.dfs(results,output='/user/yleborgn/bridge/snpsMat10000.results')
  
  pv<-sapply(results[[2]],function(x) x$pv)
  st<-sapply(results[[2]],function(x) x$st)
  #genes<-sapply(results[[2]],function(x) x$gene.rs)
  
  #pv<-p.adjust(pv)
  print(summary(st))
  #st[which(is.na(st))]<-Inf
  
  rr<-apply(snpsMat,1,countNA)
  to.remove<-which(rr>0)
  
  st<-st[-to.remove]
  snpsMat2<-snpsMat[-to.remove,]
  names.locus<-results[[1]]
  names.locus<-names.locus[-to.remove]
  
  isel<-sort(st,decreasing=FALSE,index=TRUE)$ix
  
  score<-cbind(names.locus[isel[1:100]],st[isel[1:100]])
  match.snpsMat<-match(score[,1],rownames(snpsMat))
  snpsMatsub<-snpsMat[match.snpsMat,]
  save(file="../../res.Rdata",score,snpsMatsub)
  
  
  
  geneid_rsid<-results[[1]]
  geneid_rsid<-sapply(geneid_rsid,strsplit,'_')
  geneid_rsid<-unlist(geneid_rsid,use.name=F)
  geneid_rsid<-matrix(geneid_rsid,nrow(X),2,byrow=T)
  geneid<-unique(geneid_rsid[,1])
  
  #names.genes<-names(gene2rs)
  
  
  res<-cbind(chosen,st)
  
  res.rs<-NULL
  for (i in 1:100) {
    submat.i<-which(geneid_rsid[,1]==chosen[i])
    X.sub<-X[submat.i,]
    rownames(X.sub)<-geneid_rsid[submat.i,2]
    res.rs<-c(res.rs,list(X.sub))
  }
  
  names(res.rs)<-chosen[1:100]
  
  write.table(file="../../res.txt",res,res.rs)
  save(file="../../res.Rdata",res)
  
}

getCountsTab<-function(M) {
  apply(M,1,paste,collapse=";")
}

postProcessing<-function() {
  
  results<-from.dfs(PVmr)
  
  score<-sapply(results[[2]],function(x) x$score)
  tab<-lapply(results[[2]],function(x) x$tab)
  names.locus<-results[[1]]
  
  isel<-sort(score,decreasing=TRUE,index=TRUE)$ix
  
  names.locus<-names.locus[isel]
  score<-score[isel]
  tab<-tab[isel]
  
  congroups <- dbConnect(RSQLite::SQLite(), "../../groupsToComparePartial.db")
  rs<-dbSendQuery(congroups,"select * from dbSNPs" )
  dbSNPs<-fetch(rs, n=-1)
  
  names.locus.mat<-matrix(unlist(sapply(names.locus,strsplit,"_")),length(names.locus),2,byrow=T)
  
  snpMatch<-match(names.locus.mat[,2],dbSNPs$Locus)
  
  countTabCols<-t(sapply(tab,getCountsTab))
  
  matRes<-cbind(score,
                dbSNPs[snpMatch,c("Locus","dbsnp_id_137","gene_ensembl","gene_symbol","reference",
                                  "alternative")],
                countTabCols,
                dbSNPs[snpMatch,c("num_genes",'snpeff_effect', 'snpeff_impact',
                                  "cadd_phred","cadd_raw","vest_score",
                                  "pph2_hdiv_score","pph2_hdiv_pred","pph2_hvar_score","pph2_hvar_pred",
                                  "sift_score","sift_pred")])
  
  niceNames<-c("Ranking score","Locus","dbSNP ID (137)","Gene Ensembl ID","Gene Symbol","Ref",
               "Alt","Alleles (Control)","Alleles (Pathological)","#Genes",'snpeff_effect', 'snpeff_impact',
               "cadd_phred","cadd_raw","VEST score",
               "pph2_hdiv_score","pph2_hdiv_pred","pph2_hvar_score","pph2_hvar_pred",
               "sift_score","sift_pred")
  
  colnames(matRes)<-niceNames
  
  match.snpsMat<-match(names.locus,rownames(snpsMat))
  snpsMat<-snpsMat[match.snpsMat,]
  
  metadata<-list(timestart="2015-05-20 15:59:41 CEST",
                 timeend="2015-05-20 17:47:38 CEST",
                 controlgroup="Erasme_Control",
                 pathogroup="NEURODEV")
  
  name<-"NEURODEV_vs_ErasmeControl"
  
  res<-list(matRes=matRes,snpsMat=snpsMat,metadata=metadata,tab=tab,name=list(name))
  
  save(file="../../res.Rdata",res)
}
