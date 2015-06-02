# Set this to where Spark is installed
Sys.setenv(SPARK_HOME="/home/yleborgn/spark")
# This line loads SparkR from the installed directory
.libPaths(c(file.path(Sys.getenv("SPARK_HOME"), "R", "lib"), .libPaths()))
library(SparkR)

# Initialize SparkContext and SQLContext
sc <- sparkR.init(appName="SparkR-DataFrame-example")
sqlContext <- sparkRSQL.init(sc)

# Create a simple local data.frame
localDF <- data.frame(name=c("John", "Smith", "Sarah"), age=c(19, 23, 18))

# Convert local data frame to a SparkR DataFrame
df <- createDataFrame(sqlContext, localDF)

# Print its schema
printSchema(df)
# root
#  |-- name: string (nullable = true)
#  |-- age: double (nullable = true)

# Create a DataFrame from a JSON file
path <- file.path(Sys.getenv("SPARK_HOME"), "examples/src/main/resources/people.json")
peopleDF <- jsonFile(sqlContext, path)
printSchema(peopleDF)

# Register this DataFrame as a table.
registerTempTable(peopleDF, "people")

# SQL statements can be run by using the sql methods provided by sqlContext
teenagers <- sql(sqlContext, "SELECT name FROM people WHERE age >= 13 AND age <= 19")

# Call collect to get a local data.frame
teenagersLocalDF <- collect(teenagers)

# Print the teenagers in our dataset 
print(teenagersLocalDF)

# Stop the SparkContext now
sparkR.stop()

#######################
sc <- sparkR.init("local")
lines <- textFile(sc, "hdfs://Users/yleborgn/test.txt")
wordsPerLine <- lapply(lines, function(line) { length(unlist(strsplit(line, " "))) })



