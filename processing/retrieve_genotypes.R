#require(bigrquery)
require(VariantAnnotation)

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

geneIDs<-c("TUBB2A","TCL1A")

getPosFromGene<-function(geneID) {
  
  
  
}

#From https://www.ebi.ac.uk/gwas/search?query=iga%20nephropathy
chromosome<-c("17")
base_range<-c(7460000,7490000)

list_variants<-c("rs3803800")


getSNPsLocalVCF1000Genomes<-function(chr, range, patients) {
  #To update for selecting chr range patients fields
  
  path_vcf<-"~/bridge/data/1000genomes/originalVCF/"
  
  rng <- GRanges(seqnames="17", ranges=IRanges(
    start=c(base_range[1]),
    end=c(base_range[2]),
    names=c("range_analysis")))
  tab <- TabixFile(paste0(path_vcf,"mod_ALL.chr17.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"))
  vcf_rng <- readVcf(tab, "hg19", param=rng)
  
}


getSNPsInfosHighlander<-function(chromosome,base_range, samples) {
  
  require(RMySQL)
  source("../connectHighlander.R")
  
  samples_select<-paste(samples,collapse="' OR patient='")
  
  sql<-paste0("
  select 
  patient,chr,pos,read_depth,
  allelic_depth_ref,allelic_depth_alt,
  zygosity,genotype_quality,dbsnp_id_137,
  dbsnp_id_141,filters,cadd_phred,cadd_raw,vest_score
  from exomes_ug
  where (patient='", samples_select,"') 
  and chr=",chromosome,"
  and pos>",base_range[1],"
and pos<",base_range[2])
  
  rs = dbSendQuery(highlanderdb,sql )
  data = fetch(rs, n=-1)
  
  id<-paste(data[,"chr"], data[,"pos"],sep="_")
  data<-cbind(id=id,data)
  data[,'id']<-as.character(data[,'id'])
  data<-data[sort(data[,'id'],index.r=T)$ix,]
  
  data
}

getGenoFromSNPsInfos<-function(SNPsInfos) {
  
  geno<-table(SNPsInfos[,"patient"],SNPsInfos[,'id'])
  hetero.id<-which(SNPsInfos[,'zygosity']=="Heterozygous")
  if (length(hetero.id)>0) {
    for (i in 1:length(hetero.id)) {
      geno[SNPsInfos[hetero.id[i],"patient"],SNPsInfos[hetero.id[i],"id"]]<-2
    }
  }
  
  geno
  
}

dummy<-function() {
  
  load(file="web/mySampleGroups.Rdata")
  
  SNPsInfos1<-getSNPsInfosHighlander(chromosome, base_range,mySampleGroups[[2]][[1]])#Patho
  SNPsInfos2<-getSNPsInfosHighlander(chromosome, base_range,mySampleGroups[[2]][[2]])#Control
  
  variantMat1<-getGenoFromSNPsInfos(SNPsInfos1)
  variantMat2<-getGenoFromSNPsInfos(SNPsInfos2)
  
  save(file="Web/variantMat.Rdata",variantMat1,variantMat2)
  
}

