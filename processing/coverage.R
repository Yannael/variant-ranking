loadCoverage<-function(num) {
  num<-1
  if (num<10) add0<-"0" else add0<-""
  filename<-paste0("/home/yleborgn/bridge/data/erasme/coverage/ISDBMbatch_coverage_chr",add0,num,"_trios")
  
  system.time(CoverageTab<-read.table(filename,nrow=10,header=T))
  
}