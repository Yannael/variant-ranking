project <- "caramel-brook-93006"


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
  #BAJA
  name_vcf<-"~/bridge/data/erasme/dataset1/FISDBM261950_hc_all_recalibrated_variants.vcf.gz"
  
  rng <- GRanges(seqnames="17", ranges=IRanges(
    start=c(base_range[1]),
    end=c(base_range[2]),
    names=c("range_analysis")))
  tab <- TabixFile(name_vcf)
  vcf <- readVcf(name_vcf, "hg19")
  
  name_vcf<-"~/bridge/data/erasme/dataset1/FISDBM261950_hc_all_recalibrated_variants.vcf.gz"
  vcf <- readVcf(name_vcf, "hg19")
  
  data<-rowData(vcf_rng)
  #i.keep<-which((seqnames(data)=="chr17") & (ranges(data)@start>7460000) & (ranges(data)@start <7480000))
  alt1<-which(sapply(info(vcf_rng)$AF,length)==1)
  i.keep<-which((seqnames(data)=="17") & (ranges(data)@start>7460000) & (ranges(data)@start <7480000) )
  i.keep<-intersect(i.keep,alt1)
  i.keep<-which(unlist(info(vcf_rng)$AF[i.keep])>0.01)
  
  dataRanges<-data[i.keep,]
  start<-ranges(data[i.keep,])@start
  wdth<-ranges(data[i.keep,])@width
  rsnames<-names(data[i.keep,])
  refchr<-as.character(ref(vcf_rng)[i.keep])
  altchr<-as.character(unlist(alt(vcf_rng)[i.keep]))
  genoGT<-(geno(vcf_rng)$GT[i.keep,])
  #genoAD<-(geno(vcf)$AD[i.keep,1])
  #genoDP<-(geno(vcf)$DP[i.keep,1])
  #genoGQ<-(geno(vcf)$GQ[i.keep,1])
  #genoPL<-(geno(vcf)$PL[i.keep,1])
  
  #DF<-cbind(start,wdth,rsnames,refchr,altchr,genoGT,genoAD,genoDP,genoGQ,genoPL)
  DF<-cbind(start,wdth,rsnames,refchr,genoGT)
  
  
  
  load("snps.erasme.Rdata")
  
  
  for (i in 1:length(list.snps)) {
    
    variants.patho<-colnames(list.snps[[i]])
    GT <-geno(vcf_rng)$GT
    variants.control<-rownames(GT)
    vv<-intersect(variants.patho,variants.control)
  }
}

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


getPheno1000GenomesGoogle<-function() {
  project<-"caramel-brook-93006"
  project<-"lyrical-line-92909"
  require(bigrquery)
  sql<-"
  select call.call_set_name as callname
from [genomics-public-data:1000_genomes_phase_3.variants]
  where reference_name='17'
and start=16049958
  "
  data<-query_exec(sql,project)  
  
  sql<-"
  SELECT
  COUNT(sample) AS all_samples,
  SUM(IF(In_Phase1_Integrated_Variant_Set = TRUE, 1, 0)) AS samples_in_variants_table
FROM
  [genomics-public-data:1000_genomes.sample_info]"
  data<-query_exec(sql,project)  
  
}


