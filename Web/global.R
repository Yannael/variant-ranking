require("rvest")
require("magrittr")

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

genPhenotypeData<-function() {

  #require(VariantAnnotation)
  
Erasme.pheno.id<-read.table("web/Erasme.exomes.id.csv",head=TRUE,sep=";",stringsAsFactors=F)
phenotypes<-Erasme.pheno.id[order(Erasme.pheno.id[,1]),]
phenotypes<-Erasme.pheno.id[Erasme.pheno.id[,5]!="NON",]
phenotypes<-cbind(phenotypes,rep("EUR",nrow(phenotypes)))
colnames(phenotypes)<-c("Sample ID", "Pathology", "InSilico ID","Exome Nb","Control","Comment","Data source","Gender","Population","Super population")
phenotypesErasme<-phenotypes[,c("Data source","Sample ID","Pathology", "Gender","Super population","Population", "InSilico ID","Comment")]
phenotypesErasme[,"Data source"]<-"Erasme"

load("web/samplesID1000Gen.Rdata")
pedigre_data<-read.table("web/20130606_g1k.ped",sep="\t",header=TRUE,stringsAsFactors=F)

require(rvest)
require(magrittr)
pop_spec<-html("http://www.1000genomes.org/category/frequently-asked-questions/population")
pop_spec<-pop_spec %>% html_nodes("td")%>%html_text%>%matrix(ncol=6,byrow=T)
match.pop.spop<-pop_spec[-1,c(1,3)]

ID.select<-intersect(pedigre_data[,2],samplesID1000Gen)
phenotypes1000Gen<-pedigre_data[which(pedigre_data[,2] %in% ID.select),c("Individual.ID","Gender","Population")]
N<-nrow(phenotypes1000Gen)
phenotypes1000Gen<-cbind(phenotypes1000Gen,rep("",N),rep("1000 Genomes",N),rep("None",N),rep("",N),rep("",N))
colnames(phenotypes1000Gen)<-c("Sample ID","Gender","Population","Super population","Data source","Pathology","InSilico ID","Comment")

phenotypes1000Gen[,"Super population"]<-match.pop.spop[match(phenotypes1000Gen$Population,match.pop.spop),2]
phenotypes1000Gen<-phenotypes1000Gen[,c("Data source","Sample ID","Pathology", "Gender","Super population","Population", "InSilico ID","Comment")]

phenotypes1000Gen[phenotypes1000Gen[,"Gender"]==1,"Gender"]<-"Male"
phenotypes1000Gen[phenotypes1000Gen[,"Gender"]==2,"Gender"]<-"Female"

save(file="web/phenotypes.Rdata",phenotypesErasme,phenotypes1000Gen)

Hartsfield<-phenotypesErasme[which(phenotypesErasme[,"Pathology"]=="Hartsfield"),"Sample ID"]
Hydrocephalus<-phenotypesErasme[which(phenotypesErasme[,"Pathology"]=="hydrocephalus"),"Sample ID"]
IGAN<-phenotypesErasme[which(phenotypesErasme[,"Pathology"]=="IGAN"),"Sample ID"]
TransplantTolerance<-phenotypesErasme[which(phenotypesErasme[,"Pathology"]=="TransplantTolerance"),"Sample ID"]

genomes1000EUR<-phenotypes1000Gen[which(phenotypes1000Gen[,'Super population']=="EUR"),"Sample ID"]

mySampleGroups<-list(c("Erasme_Hartsfield","Erasme_Hydrocephalus","Erasme_IGAN","Erasme_TransplantTolerance","1000genomes_EUR"),
                     list(Hartsfield,Hydrocephalus,IGAN,TransplantTolerance,genomes1000EUR))


save(file="web/mySampleGroups.Rdata",mySampleGroups)
}

load(file="phenotypes.Rdata")
phenotypesAll<-rbind(phenotypesErasme,phenotypes1000Gen)

username<-'yleborgn'

load("mySampleGroups.Rdata")  
  
load("myResults.Rdata")

#disGenNet<-read.csv("all_gene_disease_associations.txt",sep="\t")
