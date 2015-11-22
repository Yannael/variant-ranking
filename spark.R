#Sys.setenv(SPARK_CLASSPATH="/Users/yalb/Projects/Github/variant-ranking/sqlite-jdbc-3.8.11.1.jar")
#Sys.setenv(SPARK_HOME="/home/yleborgn/spark")
Sys.setenv(DYLD_FALLBACK_LIBRARY_PATH="/Library/Java/JavaVirtualMachines/jdk1.7.0_80.jdk/Contents/Home/jre/lib/server/")
Sys.setenv(HADOOP_JAR="/Users/yalb/spark/assembly/target/scala-2.10/spark-assembly-1.5.0-hadoop2.4.0.jar")
# This line loads SparkR from the installed directory
Sys.setenv(SPARK_HOME="/Users/yalb/spark")
.libPaths(c(file.path(Sys.getenv("SPARK_HOME"), "R", "lib"), .libPaths()))
library(SparkR)

#Sys.setenv('SPARKR_SUBMIT_ARGS'='"--packages" "com.databricks:spark-csv_2.10:1.2.0" "sparkr-shell"')

sparkEnvir <- list('spark.sql.parquet.binaryAsString'='true')
sc<-sparkR.init(master="local[2]",sparkEnvir=sparkEnvir)
sqlContext <- sparkRSQL.init(sc)
sc <- sparkR.init(master="local")
hiveContext <- sparkRHive.init(sc)

#sc<-sparkR.init()
#sc <- sparkR.init(master="local",sparkJars="sqlite-jdbc-3.8.11.1.jar",sparkPackages ="com.databricks:spark-csv_2.11:1.2.0")

DF2 <- loadDF(sqlContext, "file:///Users/yalb/Projects/Github/variant-ranking/variants")  
DF2 <- loadDF(sqlContext, "file:///Users/yalb/Projects/Github/variant-ranking/genotypes.db/genotypes")  
dim(DF2)
registerTempTable(DF2, "DF")
pos<-sql(sqlContext, "SELECT dbsnp_id_137 FROM DF WHERE position >= 13 AND position <= 100000 limit 2")
pos<-sql(sqlContext, "SELECT * FROM DF where position<1000000 limit 10000")
pos<-sql(sqlContext, "SELECT * FROM DF where position<10000 limit 10000")
pos<-sql(sqlContext, "SELECT * FROM DF where v1=2 and v2=2 and v3=2 and v4=2 limit 1000000")
system.time(data<-collect(pos))

library(RJDBC)
drv <- JDBC(driverClass = "org.apache.hive.jdbc.HiveDriver",
            classPath = list.files("/Users/yalb/Downloads/impala-jdbc-0.5-2",pattern="jar$",full.names=T),
            identifier.quote="`")
conn <- dbConnect(drv, "jdbc:hive2://192.168.99.101:21050/;auth=noSasl")
conn <- dbConnect(drv, "jdbc:hive2://127.0.0.1:21050/;auth=noSasl")
system.time(data <- dbGetQuery(conn, "select * from gvr4.variants where position<10000 limit 100000"))
system.time(data <- dbGetQuery(conn, "select distinct sample_id from gvr4.variants "))


system.time(df <- read.df(sqlContext, "genoMat_1000000_1000_100.txt", source = "com.databricks.spark.csv"))
system.time(dim(df))

parseFields <- function(record) {
  Sys.setlocale("LC_ALL", "C") # necessary for strsplit() to work correctly
  parts <- strsplit(record, " ")[[1]]
  list(id=parts[1], title=parts[2], modified=parts[3], text=parts[4], username=parts[5])
}
parsedRDD <- lapply(df, parseFields)

condb<-dbConnect(RSQLite::SQLite(), paste0("data.db"))
data<-dbGetQuery(condb,paste0("select * from variants "))
dbDisconnect(condb)

save(file="data.Rdata",data)

assign("classpath","/Users/yalb/Projects/Github/variant-ranking/sqlite-jdbc-3.8.11.1.jar")

df <- loadDF(sqlContext, source="jdbc", url="jdbc:sqlite:samplesSets.db", dbtable="schema.groups")

system.time(DF <- createDataFrame(sqlContext, data[1:4000,]))


saveDF((DF, "namesAndAges.parquet", "parquet"))

DF2 <- loadDF(sqlContext, "namesAndAges.parquet")  


destDF <- filter(DF,DF$Data_Source== "ULB")


library(dplyr)
library(RJDBC)
library(dplyr.spark)

library(rJava)
.jinit()
spark.src = src_SparkSQL("localhost", "10000")
my_db<-src_SparkSQL("localhost", "10000")

DF<-data[1:1000,1:30]
system.time(copy_to(my_db, DF, temporary = FALSE))

DF2<-tbl(my_db,"DF")
cache(DF2)

c1 = filter(DF2, chr=2)

#install.packages("nycflights13")

my_db = src_SparkSQL()

genoMat<-read.table("geno100m30.txt")
DF<-genoMat[1:100000,1:30]

#data 41 col
system.time(DF <- createDataFrame(sqlContext, data[1:10000,]))
#12s
system.time(destDF <- select(DF, "Data_Source"))
#0.034s
system.time(aa<-collect(destDF))
#10s


system.time(df <- createDataFrame(sqlContext, DF52) )


my_db <- src_sqlite("data.db")
system.time(db<-tbl(my_db,"variants"))
system.time(dd<-filter(db, position < "2005"))
system.time(dim(dd))

system.time(dd<- sql(my_db,"SELECT * FROM variants where pos<1000"))

