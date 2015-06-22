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

#From https://www.ebi.ac.uk/gwas/search?query=iga%20nephropathy
chromosome<-c("17")
base_range<-c(7440000,7490000)

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


getSNPsInfosHighlander<-function(samples,chromosome=NULL,base_range=NULL) {
  
  if (is.null(chromosome)) chr_select<-""
  else chr_select<-paste0("and chr=",chromosome)
  
  if (is.null(base_range)) base_range_select<-""
  else base_range_select<-paste0("and pos>=",base_range[1]," and pos<=",base_range[2])
  
  require(RMySQL)
  source("../connectHighlander.R")
  
  samples_select<-paste(samples,collapse="' OR patient='")
  
  sql<-paste0("
  select 
  patient,chr,pos,read_depth,
  zygosity,reference, alternative, gene_symbol,
  gene_ensembl, num_genes,clinvar_rs, 
  dbsnp_id_137,
  allelic_depth_ref,
  allelic_depth_alt,
  snpeff_effect, snpeff_impact,
  consensus_maf,
  cadd_phred,cadd_raw,vest_score,pph2_hdiv_score, pph2_hdiv_pred, 
  pph2_hvar_score, pph2_hvar_pred, sift_score, sift_pred, short_tandem_repeat
  from exomes_ug
  where 
        filters='PASS' and
        read_depth>15 and
        mapping_quality_zero_reads=0 and
        downsampled='false' and
        allele_num=2 and
        genotype_quality>70 and
        
              (patient='", samples_select,"') ",
              chr_select,
              base_range_select
  )
  #change_type='SNP' and 
  rs = dbSendQuery(highlanderdb,sql )
  data = fetch(rs, n=-1)
  
  dbClearResult(rs)
  dbDisconnect(highlanderdb)
  
  data
}


createSNPsHighlanderDB<-function() {
  #Load all data from Erasme trios. 
  
  load(file="Web/mySampleGroups.Rdata")
  
  system.time(
    SNPs_patho<-getSNPsInfosHighlander(mySampleGroups[[2]][[1]],chr="1") #Patho
  )
  system.time(
    SNPs_control<-getSNPsInfosHighlander(mySampleGroups[[2]][[2]],chr="1") #Control
  )
  
  chr<-SNPs_patho[,"chr"]
  pos<-SNPs_patho[,"pos"]
  id_SNPs_patho<-paste(chr,pos,sep=":")
  SNPs_patho$Locus<-id_SNPs_patho
  
  chr<-SNPs_control[,"chr"]
  pos<-SNPs_control[,"pos"]
  id_SNPs_control<-paste(chr,pos,sep=":")
  SNPs_control$Locus<-id_SNPs_control
  
  id_control_to_remove<-setdiff(id_SNPs_control,id_SNPs_patho)
  i<-which(id_SNPs_control %in% id_control_to_remove)
  if (length(i)>0) SNPs_control<-SNPs_control[-i,]
  
  
  dbSNPs<-SNPs_patho[,c('Locus','chr', 'pos', 'dbsnp_id_137','gene_ensembl', 'gene_symbol',
                        'reference','alternative','num_genes','clinvar_rs',
                        'snpeff_effect', 'snpeff_impact',
                        'consensus_maf',
                        'cadd_phred','cadd_raw','vest_score','pph2_hdiv_score','pph2_hdiv_pred', 
                        'pph2_hvar_score', 'pph2_hvar_pred', 'sift_score', 'sift_pred', 'short_tandem_repeat'
                     )]
  
  SNPs_patho<-SNPs_patho[,c("patient","chr","pos","Locus","read_depth","zygosity",'alternative')]
  SNPs_control<-SNPs_control[,c("patient","chr","pos","Locus","read_depth","zygosity",'alternative')]
  
  system.time({
    con <- dbConnect(RSQLite::SQLite(), "groupsToComparePartial.db")
    dbWriteTable(con,"patho",SNPs_patho,overwrite=T)
    dbWriteTable(con,"control",SNPs_control,overwrite=T)
    dbWriteTable(con,"dbSNPs",dbSNPs,overwrite=T)
    dbDisconnect(con)
  })
  
}



