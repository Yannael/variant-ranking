require("rvest")
require("magrittr")


load(file="phenotypes.Rdata")
phenotypesAll<-rbind(phenotypesErasme,phenotypes1000Gen)

username<-'yleborgn'

load("mySampleGroups.Rdata")  

load("myResults.Rdata")

#disGenNet<-read.csv("all_gene_disease_associations.txt",sep="\t")
