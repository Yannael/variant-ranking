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

getvariantsLocalVCF1000Genomes<-function(chr, range, patients) {
  #To update for selecting chr range patients fields
  
  path_vcf<-"~/bridge/data/1000genomes/originalVCF/"
  
  rng <- GRanges(seqnames="17", ranges=IRanges(
    start=c(base_range[1]),
    end=c(base_range[2]),
    names=c("range_analysis")))
  tab <- TabixFile(paste0(path_vcf,"mod_ALL.chr17.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"))
  vcf_rng <- readVcf(tab, "hg19", param=rng)
  
}


getvariantsInfosHighlander<-function(samples,chromosome=NULL,base_range=NULL) {
  
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
  zygosity,reference, alternative, change_type, gene_symbol,
  gene_ensembl, num_genes,clinvar_rs, 
  consensus_MAF, consensus_MAC,
  dbsnp_id_137,
  allelic_depth_ref,
  allelic_depth_alt,
  snpeff_effect, snpeff_impact,
  cadd_phred,cadd_raw,vest_score,pph2_hdiv_score, pph2_hdiv_pred, 
  pph2_hvar_score, pph2_hvar_pred, sift_score, sift_pred, short_tandem_repeat
  from exomes_ug
  where 
        filters='PASS' and
        read_depth>10 and
        mapping_quality_zero_reads=0 and
        downsampled='false' and
        allele_num=2 and
        genotype_quality>70 and
        snpeff_effect<>'SYNONYMOUS_CODING' and 
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


createvariantsHighlanderDB<-function() {
  #Load all data from Erasme trios. 
  
  load(file="Web/mySampleGroups.Rdata")
  
#   system.time(
#     variants_patho<-getvariantsInfosHighlander(mySampleGroups[[2]][[1]],chr="1") #Patho
#   )
#   system.time(
#     variants_control<-getvariantsInfosHighlander(mySampleGroups[[2]][[2]],chr="1") #Control
#   )
  system.time(
    variants_patho<-getvariantsInfosHighlander(mySampleGroups[[2]][[1]]) #Patho
  )
  system.time(
    variants_control<-getvariantsInfosHighlander(mySampleGroups[[2]][[2]]) #Control
  )
  
  chr<-variants_patho[,"chr"]
  pos<-variants_patho[,"pos"]
  id_variants_patho<-paste(chr,pos,sep=":")
  variants_patho$Locus<-id_variants_patho
  
  chr<-variants_control[,"chr"]
  pos<-variants_control[,"pos"]
  id_variants_control<-paste(chr,pos,sep=":")
  variants_control$Locus<-id_variants_control
  
  id_control_to_remove<-setdiff(id_variants_control,id_variants_patho)
  i<-which(id_variants_control %in% id_control_to_remove)
  if (length(i)>0) variants_control<-variants_control[-i,]
  
  
  dbvariants<-variants_patho[,c('Locus','chr', 'pos', 'dbsnp_id_137','gene_ensembl', 'gene_symbol',
                        'reference','alternative','change_type','num_genes','clinvar_rs',
                        'snpeff_effect', 'snpeff_impact',
                        'consensus_MAC','consensus_MAF',
                        'cadd_phred','cadd_raw','vest_score','pph2_hdiv_score','pph2_hdiv_pred', 
                        'pph2_hvar_score', 'pph2_hvar_pred', 'sift_score', 'sift_pred', 'short_tandem_repeat'
                     )]
  
  variants_patho<-variants_patho[,c("patient","chr","pos","Locus","read_depth","zygosity",'reference','alternative')]
  variants_control<-variants_control[,c("patient","chr","pos","Locus","read_depth","zygosity",'reference','alternative')]
  
  system.time({
    con <- dbConnect(RSQLite::SQLite(), "groupsToComparePartial.db")
    dbWriteTable(con,"patho",variants_patho,overwrite=T)
    dbWriteTable(con,"control",variants_control,overwrite=T)
    dbWriteTable(con,"dbvariants",dbvariants,overwrite=T)
    dbDisconnect(con)
  })
  
}



