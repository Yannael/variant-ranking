getPheno1000Genome<-function() {
  
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
  phenotypes1000Gen
}

getPhenoHighlander<-function() {
  #Connect to Highlander DB
  require(RMySQL)
  source("../../connectHighlander.R")
  
  #print(dbListTables(highlanderdb))
  #dbListFields(highlanderdb, 'exomes_ug')
  
  rs = dbSendQuery(highlanderdb, "select distinct patient,pathology from exomes_ug")
  data = fetch(rs, n=-1)
  
  N<-nrow(data)
  phenotypesErasme<-cbind(rep("Erasme",N),data[,1],data[,2],rep("",N),rep("EUR",N),rep("",N))
  colnames(phenotypesErasme)<-c("Data source","Sample ID","Pathology", "Gender","Super population","Population")
  
  phenotypesErasme[phenotypesErasme[,"Pathology"]=="CTRL","Pathology"]<-"None"
  phenotypesErasme
}

getPhenoData<-function() {
  
  phenotypes1000Gen<-getPheno1000Genome()
  phenotypesErasme<-getPhenoHighlander()
  
  save(file="web/phenotypes.Rdata",phenotypesErasme,phenotypes1000Gen)
  
}

createPhenoGroups<-function() {
  
  NEURODEV<-phenotypesErasme[which(phenotypesErasme[,"Pathology"]=="NEURODEV"),"Sample ID"]
  Control<-phenotypesErasme[which(phenotypesErasme[,"Pathology"]=="None"),"Sample ID"]
  
  genomes1000EUR<-phenotypes1000Gen[which(phenotypes1000Gen[,'Super population']=="EUR"),"Sample ID"]
  
  
  mySampleGroups<-list(c("Erasme_NEURODEV","Erasme_Control","1000genomes_EUR"),
                       list(NEURODEV,Control,genomes1000EUR))
  
  
  save(file="web/mySampleGroups.Rdata",mySampleGroups)
}


