require(bigrquery)
project <- "caramel-brook-93006"

#From https://www.ebi.ac.uk/gwas/search?query=iga%20nephropathy
chromosome<-c("17")
base_range<-c(7550000,7590000)

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
  
  
  
  
}

