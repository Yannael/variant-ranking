source("/home/jpoullet/Rscript/testing/config.R")
library('rhdfs')
library('rmr2')
hdfs.init()

infile<-"in.variantMat.csv"
outfile<-paste0(gsub(" ","_",Sys.time()),"_out.variantMat.csv")

local<-TRUE
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

entropy.tab<-function(tab){
  if (sum(tab)==0)
    return(10)
  S<-apply(tab,2,sum)
  S[which(S==0)]<-1
  px<-S/sum(S)
  px[which(is.na(px))]<-0
  p<-NULL
  H<-NULL
  
  w<-which(px>0)
  
  for (i in w){
    if (any(tab[,i]>0)){
      pp<-tab[,i]/sum(tab[,i])
      lp<-log2(pp)
      ##E<-as(entropy.shrink(tab[,i],lambda.freqs=0.1,verbose=FALSE),"numeric")
      wp<-which(pp>0)
      E<--sum(pp[wp]*lp[wp])
      ##E<-sum(pp)
      H<-c(H,E)
    }else{
      H<-c(H,NA)
    }
  }
  sum(px[w]*H,na.rm=TRUE)
}

createData<-F
if (createData) {
prof.init<-system.time({
  load("Web/variantMat.Rdata")
  
  rs<-colnames(variantMat1)
  rs2<-colnames(variantMat2)
  
  #We are only interested in variants in pathological population (pop1)
  variantMat2_aligned<-matrix(NA,nrow(variantMat2),ncol(variantMat1))
  subsetColMat2<-intersect(rs2,rs)

  matchMat2Mat1<-match(subsetColMat2,rs)
  
  variantMat2_aligned[,matchMat2Mat1]<-variantMat2[,subsetColMat2]
  
  rownames(variantMat1)<-rep("patho",nrow(variantMat1))
  rownames(variantMat2_aligned)<-rep("control",nrow(variantMat2))
  
  rs.geno.matrix.expe<-rbind(variantMat1,variantMat2_aligned)
  genes<-paste("gene",1:length(rs),sep="")
  idkeys<-paste0(genes,"_keysplit_",rs)
  
  colnames(rs.geno.matrix.expe)<-idkeys
  
  prof.savehdfs<-system.time(to.dfs(t(rs.geno.matrix.expe),output=hdfsfile))
  
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
  
  if (sum(tab["case",])>(length(cases)/2) & sum(tab["ctrl",])>0){    
    ttest<- prop.trend.test(tab["case",], margin.table(tab,2))
    pv<-ttest$p.value ##pvalue(IT) ##c(L.case,L.control))
  }
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
  geneid_rsid<-sapply(geneid_rsid,strsplit,'_keysplit_')
  geneid_rsid<-unlist(geneid_rsid,use.name=F)
  geneid_rsid<-matrix(geneid_rsid,nrow(X),2,byrow=T)
  
  N<-nrow(X)
  
  controls<-which(colnames(X)=="control")
  cases<-which(colnames(X)=="patho")
  
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


prof.run<-system.time({
  PVmr = mapreduce(
    input=hdfsfile,  
    # input.format=dat.in.format,
    map = map,
    reduce=reduce
  )
})

results<-from.dfs(PVmr)
#to.dfs(results,output='/user/yleborgn/bridge/out16000.results')

pv<-sapply(results[[2]],function(x) x$pv)
st<-sapply(results[[2]],function(x) x$st)

pv<-p.adjust(pv)
print(summary(st))
st[which(is.na(st))]<-Inf
isel<-sort(st,decreasing=FALSE,index=TRUE)$ix[1:(10)]
#names.genes<-names(gene2rs)
#chosen<-names.genes[isel]


