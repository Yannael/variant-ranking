Sys.setenv(SPARK_HOME="/Users/yalb/spark")
#Sys.setenv(SPARK_HOME="/home/yleborgn/spark")
Sys.setenv(DYLD_FALLBACK_LIBRARY_PATH="/Library/Java/JavaVirtualMachines/jdk1.7.0_79.jdk/Contents/Home/jre/lib/server/")
Sys.setenv(HADOOP_JAR="/Users/yalb/spark/assembly/target/scala-2.10/spark-assembly-1.4.1-hadoop2.4.0.jar")
# This line loads SparkR from the installed directory
.libPaths(c(file.path(Sys.getenv("SPARK_HOME"), "R", "lib"), .libPaths()))
library(SparkR)
sc <- sparkR.init(master="local")

sqlContext <- sparkRSQL.init(sc)

condb<-dbConnect(RSQLite::SQLite(), paste0("data.db"))
data<-dbGetQuery(condb,paste0("select * from variants "))
dbDisconnect(condb)

save(file="data.Rdata",data)

system.time(DF <- createDataFrame(sqlContext, data[1:1000,]))

system.time(flightsDF <- read.df(sqlContext, flightsCsvPath, source = "com.databricks.spark.csv", header = "true"))

destDF <- select(DF, "Data_Source", "ULB")


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

dd<- sql(my_db,"SELECT * FROM variants")

