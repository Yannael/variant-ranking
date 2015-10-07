dummy<-function() {

data<-read.table("genoMatTrios.txt")
#40273 variants
nVarPerGene<-tapply(data[,1],data[,1],length)
length(nVarPerGene)
#15329 genes

dd<-data.frame(nVarPerGene=nVarPerGene)
ggplot(data=dd,aes(nVarPerGene)) + geom_histogram()

#De novo
time<-c(0.7,1.4,2.2,2.5,4.6,8.8)
nsamples<-as.factor(c(30,60,120,30,60,120))
cond<-as.factor(c("200000","200000","200000","880000","880000","880000"))
DF<-data.frame(time=time,nsamples=nsamples,cond=cond)

ggplot(data=DF, aes(x=nsamples, y=time, fill=cond)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black")+
xlab("Number of samples") +
  ylab("Computation time") +
  scale_fill_hue(name="Number of variants")

#Univariate gene
#1 core

times<-c(11,76,116,73,262,420,125,403,824)
nvar<-as.factor(rep(c(100000,500000,1000000),3))
nsamples<-as.factor(rep(c(100,500,1000),each=3))
DF<-data.frame(times=times,nsamples=nsamples,nvar=nvar)

ggplot(data=DF, aes(x=nvar, y=times, fill=nsamples)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black")+
  xlab("Number of variants") +
  ylab("Computation time (s)") +
  scale_fill_hue(name="Number of samples")+
  scale_y_continuous(limits = c(0, 1000))
 # theme(legend.position="bottom")

times<-c(2,8,13,10,34,61,16,57,113)
nvar<-as.factor(rep(c(100000,500000,1000000),3))
nsamples<-as.factor(rep(c(100,500,1000),each=3))
DF<-data.frame(times=times,nsamples=nsamples,nvar=nvar)

ggplot(data=DF, aes(x=nvar, y=times, fill=nsamples)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black")+
  xlab("Number of variants") +
  ylab("Computation time (s)") +
  scale_fill_hue(name="Number of samples")+
  scale_y_continuous(limits = c(0, 1000))
  

#Gene list
genes<-read.table("ID_genelist.txt",stringsAsFactors=F)[,1]
genelist<-paste0("(",paste0(genes,collapse=","),")")

}

createGenoMat<-function(N,G,S) {
  
  geno<-matrix(base::sample(c(0:2),N*S,replace=T),N,S) 
  nPerG<-N/G
  genes<-rep(1:G,each=nPerG)
  variants<-1:N
  data<-data.frame(genes=genes,variants=variants,geno)
  
  write.table(file=paste0("genoMat_",N,"_",G,"_",S,".txt"),data,col.names=F,row.names=F,quote=F)
  
  
}
