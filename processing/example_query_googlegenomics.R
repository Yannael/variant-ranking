require(bigrquery)
project <- "caramel-brook-93006"

getBRCA1<-function() {
  
  # Retrieve sample-level information for BRCA1 variants.
  sql<-"
  SELECT
  reference_name,
  start,
  end,
  reference_bases,
  GROUP_CONCAT(alternate_bases) WITHIN RECORD AS alternate_bases,
  call.call_set_name,
  GROUP_CONCAT(STRING(call.genotype)) WITHIN call AS genotype,
  call.phaseset,
  call.genotype_likelihood,
FROM
  [genomics-public-data:1000_genomes.variants]
WHERE
  reference_name = '17'
  AND start BETWEEN 41196311
  AND 41277499
HAVING
  alternate_bases IS NOT NULL
ORDER BY
  start,
  alternate_bases,
  call.call_set_name
"
  data<-query_exec(sql,project)
}

getBRCA2<-function() {
  
  # Sample variant counts within BRCA1.
  sql<-"
SELECT
call_set_name,
COUNT(call_set_name) AS variant_count,
FROM (
  SELECT
  reference_name,
  start,
  END,
  reference_bases,
  GROUP_CONCAT(alternate_bases) WITHIN RECORD AS alternate_bases,
  call.call_set_name AS call_set_name,
  NTH(1,
      call.genotype) WITHIN call AS first_allele,
  NTH(2,
      call.genotype) WITHIN call AS second_allele,
  FROM
  [genomics-public-data:1000_genomes_phase_3.variants]
  WHERE
  reference_name = '17'
  AND start BETWEEN 41196311
  AND 41277499
  HAVING
  first_allele > 0
  OR second_allele > 0
)
GROUP BY
call_set_name
ORDER BY
call_set_name
"
  data<-query_exec(sql,project)  
}


getFirstVariant<-function() {
sql<-'SELECT
reference_name,
call.call_set_name
FROM
[genomics-public-data:1000_genomes_phase_3.variants]
limit 1
'
data<-query_exec(sql,project)
#2535 call sets.
}

test<-function() {

  
  sql<-'
  select reference_name,
start,
end,
reference_bases,
  alternate_bases,
count(call.call_set_name) WITHIN RECORD AS num_samples,
VT

from [genomics-public-data:1000_genomes.variants]

where reference_name="17"
AND start BETWEEN 41196311
  AND 41277499
  '
  data<-query_exec(sql,project)
  
  sql<-"
   SELECT
    reference_name,
    start,
  END,
  reference_bases,
  GROUP_CONCAT(alternate_bases) WITHIN RECORD AS alternate_bases,
  call.call_set_name AS call_set_name,
  NTH(1,
  call.genotype) WITHIN call AS first_allele,
  NTH(2,
  call.genotype) WITHIN call AS second_allele,
NTH(1,
  call.call_set_name) WITHIN record AS csn1,
  GROUP_CONCAT(STRING(call.genotype)) WITHIN call AS callgeno
  FROM
  [genomics-public-data:1000_genomes_phase_3.variants]
  WHERE
  reference_name = '17'
  AND start BETWEEN 41196311
  AND 41277499
  HAVING
call_set_name='HG00657'
and (
  first_allele > 0
  OR second_allele > 0)
order by start
"
  data<-query_exec(sql,project)


  
  
  
}

