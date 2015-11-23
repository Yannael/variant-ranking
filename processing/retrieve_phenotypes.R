getPheno1000Genome<-function() {
  
  load("../samplesID1000Gen.Rdata")
  pedigre_data<-read.table("../20130606_g1k.ped",sep="\t",header=TRUE,stringsAsFactors=F)
  
  require(rvest)
  require(magrittr)
  pop_spec<-read_html("http://www.1000genomes.org/category/frequently-asked-questions/population")
  pop_spec<-pop_spec %>% html_nodes("td")%>%html_text%>%matrix(ncol=6,byrow=T)
  match.pop.spop<-pop_spec[-1,c(1,3)]
  
  ID.select<-intersect(pedigre_data[,2],samplesID1000Gen)
  phenotypes1000Gen<-pedigre_data[which(pedigre_data[,2] %in% ID.select),c("Individual.ID","Gender","Population")]
  N<-nrow(phenotypes1000Gen)
  phenotypes1000Gen<-cbind(phenotypes1000Gen,rep("",N),rep("1000 Genomes",N),rep("Control",N))
  colnames(phenotypes1000Gen)<-c("Sample ID","Gender","Population","Super population","Data source","Pathology")
  
  phenotypes1000Gen[,"Super population"]<-match.pop.spop[match(phenotypes1000Gen$Population,match.pop.spop),2]
  phenotypes1000Gen<-phenotypes1000Gen[,c("Data source","Sample ID","Pathology", "Gender","Super population","Population")]
  
  phenotypes1000Gen[phenotypes1000Gen[,"Gender"]==1,"Gender"]<-"Male"
  phenotypes1000Gen[phenotypes1000Gen[,"Gender"]==2,"Gender"]<-"Female"
  phenotypes1000Gen
}

getPhenoHighlander<-function() {
  #Connect to Highlander DB
  require(RMySQL)
  source("../../connectHighlander.R")
  
  #print(dbListTables(highlanderdb))
  #dbListFields(highlanderdb, 'exomes_ug')
  
  rs = dbSendQuery(highlanderdb, "select distinct patient,pathology from exomes_hc")
  data = fetch(rs, n=-1)
  
  N<-nrow(data)
  phenotypesULB<-cbind(rep("ULB",N),data[,1],data[,2],rep("",N),rep("EUR",N),rep("",N))
  colnames(phenotypesULB)<-c("Data source","Sample ID","Pathology", "Gender","Super population","Population")
  
  phenotypesULB[phenotypesULB[,"Pathology"]=="CTRL","Pathology"]<-"Control"
  phenotypesULB
}

getPhenoData<-function() {
  
  phenotypes1000Gen<-getPheno1000Genome()
  phenotypesULB<-getPhenoHighlander()
  
  phenotypes<-rbind(phenotypesULB,phenotypes1000Gen)
  
  phenotypesdb <- dbConnect(RSQLite::SQLite(), "phenotypes.db")
  colnames(phenotypes)<-c("Data_Source","Sample_ID","Pathology","Gender","Super_Population","Population")
  dbWriteTable(phenotypesdb,"phenotypes",phenotypes,overwrite=T,row.names=F)
  dbDisconnect(phenotypesdb)
  
  data<-phenotypes
  data[,1]<-as.factor(data[,1])
  data[,2]<-as.charater(data[,2])
  data[,3]<-as.factor(data[,3])
  data[,4]<-as.factor(data[,4])
  data[,5]<-as.factor(data[,5])
  data[,6]<-as.factor(data[,6])
  filters<-getFiltersFromTable(data)
  save(file="filterPhenotypeSpec.Rdata",filters)
  
  #condb<-dbConnect(RSQLite::SQLite(), paste0("data.db"))
  #dbWriteTable(condb,"phenotypes",phenotypes,row.names=F,overwrite=T) #228MB
  #dbDisconnect(condb)
  
}

createPhenoGroups<-function() {
  
  NEURODEV<-phenotypesULB[which(phenotypesULB[,"Pathology"]=="NEURODEV"),"Sample ID"]
  Control<-phenotypesULB[which(phenotypesULB[,"Pathology"]=="Control"),"Sample ID"]
  
  genomes1000EUR<-phenotypes1000Gen[which(phenotypes1000Gen[,'Super population']=="EUR"),"Sample ID"]
  
  
  mySampleGroups<-list(c("ULB_NEURODEV","ULB_Control","1000genomes_EUR"),
                       list(NEURODEV,Control,genomes1000EUR))
  
  
  save(file="Web/mySampleGroups.Rdata",mySampleGroups)
}


