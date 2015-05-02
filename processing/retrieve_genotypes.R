require(bigrquery)
require(VariantAnnotation)

project <- "caramel-brook-93006"

#From https://www.ebi.ac.uk/gwas/search?query=iga%20nephropathy
chromosome<-c("17")
base_range<-c(7460000,7480000)

list_variants<-c("rs3803800")

getSNPsGoogleGenomics<-function() {
  DB<-"[genomics-public-data:1000_genomes.variants]"
    sql<-paste0("
  SELECT
  reference_name,
  start,
  end,
  reference_bases,
  GROUP_CONCAT(alternate_bases) WITHIN RECORD AS alternate_bases,
  call.call_set_name,
  GROUP_CONCAT(STRING(call.genotype)) WITHIN call AS genotype
  FROM ", DB," 
  WHERE
  reference_name = '",chromosome,"'
  AND start BETWEEN ",base_range[1],"
  AND ",base_range[2],"
  HAVING
  alternate_bases IS NOT NULL
  and call.call_set_name='HG00096'
  
ORDER BY
  start,
  alternate_bases,
  call.call_set_name
  ")
  data<-query_exec(sql,project)
data
  
}

getSNPsLocalVCF<-function() {
  
  path_vcf<-"~/bridge/data/1000genomes/originalVCF/"
  
  rng <- GRanges(seqnames="17", ranges=IRanges(
    start=c(base_range[1]),
    end=c(base_range[2]),
    names=c("range_analysis")))
  tab <- TabixFile(paste0(path_vcf,"mod_ALL.chr17.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"))
  vcf_rng <- readVcf(tab, "hg19", param=rng)
  
}

getSNPsErasmeVCF<-function() {
  
  rng <- GRanges(seqnames="17", ranges=IRanges(
         start=c(base_range[1]),
         end=c(base_range[2]),
         names=c("range_analysis")))
  tab <- TabixFile(paste0(path_vcf,"mod_ALL.chr17.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"))
  vcf_rng <- readVcf(tab, "hg19", param=rng)
  
  pheno<-read.table("~/bridge/data/erasme/dataset1/Erasme.exomes.id.csv",sep=";",header=T,stringsAsFactors=F)
  i_igan<-which(pheno$Disease=="IGAN")
  names_vcf<-paste0("~/bridge/data/erasme/dataset1/F",pheno$InSilicoDB.id[i_igan],"_hc_all_recalibrated_variants.vcf.gz")

  name_vcf<-names_vcf[1]
  
  rng <- GRanges(seqnames="17", ranges=IRanges(
    start=c(base_range[1]),
    end=c(base_range[2]),
    names=c("range_analysis")))
  tab <- TabixFile(name_vcf)
  vcf <- readVcf(name_vcf, "hg19")
  
  load("snps.erasme.Rdata")
  
  
  for (i in 1:length(list.snps)) {
    
    variants.patho<-colnames(list.snps[[i]])
    GT <-geno(vcf_rng)$GT
    variants.control<-rownames(GT)
    vv<-intersect(variants.patho,variants.control)
  }
  
  
  
}

