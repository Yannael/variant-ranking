
# coding: utf-8

# In[2]:

#from pyspark import SparkContext, SparkConf
#import json
#import time
#import sys


# In[11]:

from pyspark import SparkContext, SparkConf
import json
import time
import sys

analysisName=sys.argv[1]
indexControl=int(sys.argv[2])
scope=sys.argv[3]
scale=sys.argv[4]

filename="analyses/"+analysisName+"_genoMat.txt"

nPartitions=4
conf = (SparkConf()
         .setMaster("local["+str(nPartitions)+"]")
         .setAppName("Ranking")
#         .set("spark.executor.memory", "5g")
#         .set("spark.driver.memory", "5g")
#         .set("spark.python.worker.memory", "5g")
       )
sc = SparkContext(conf=conf)



# In[18]:

#analysisName="Control_vs_Transplant_Rare"
#indexControl=27
#scope="monogenic"
#scale="variant"
#filename="analyses/"+analysisName+"_genoMat.txt"
#nPartitions=4


# In[19]:

def splitValues(variantData):
    variantData=variantData.split(" ")
    variantGene=variantData.pop(0)
    variantLocus=variantData.pop(0)
    
    return (variantGene,(variantLocus,variantData))

def makePairParts(k,v,nbPart):
    result=[]
    for i in range(0,nbPart):
        result.append(((k,i),v))
        
    return [(str(sorted([k,i])),(v)) for i in range(0,nbPart)]

def f(splitIndex ,v): 
    return [(splitIndex,list(v))]


# In[20]:

def scoreVariantUnivariate(variantData):
    variantData=variantData.split(" ")
    variantGene=variantData.pop(0)
    variantLocus=variantData.pop(0)
    
    score=0
    
    sumCase=sum([int(int(x)>0) for x in variantData[0:indexControl]])
    sumControl=sum([int(int(x)>0) for x in variantData[indexControl:len(variantData)]])
    score=sumCase-sumControl
    #if sumControl>0:
    #    score=0
        
    return ((variantGene,variantLocus),(score,sumCase,sumControl))


# In[21]:

def scoreGeneUnivariate(k,variantList):
    variantList=list(variantList)
    result=[k,(0,0,0)]
    sumCase=0
    sumControl=0
    score=0
    genosum=[]
    if len(variantList)>0:
        for i in range(0,len(variantList)):
            (locus,geno)=variantList[i]
            if genosum==[]:
                genosum=[int(x) for x in geno]
            else:
                genosum=[int(x)+int(y) for x,y in zip(genosum,geno)]
                
        sumCase=sum([int(x>0) for x in genosum[0:indexControl]])
        sumControl=sum([int(x>0) for x in genosum[indexControl:len(genosum)]])
        score=sumCase-sumControl
        #if sumControl>0:
        #    score=0

    #if score>0:
    result=[k,(score,sumCase,sumControl)]
    
    return result


# In[22]:

def scoreDigenicGene(k,variantLists):
    variantLists=list(variantLists)
    result=[]
    geno1sum=[]
    geno2sum=[]
    score=0
    gene1=""
    gene2=""
    sumCase=-1
    sumControl=-1
    if len(variantLists)==2:
        (genes,variantList1)=list(variantLists[0])
        (genes,variantList2)=list(variantLists[1])
        gene1=genes[0]
        gene2=genes[1]
        variantList1=list(variantList1)
        variantList2=list(variantList2)
        for i in range(0,len(variantList1)):
            (locus1,geno1)=variantList1[i]
            if geno1sum==[]:
                geno1sum=[int(x) for x in geno1]
            else:
                geno1sum=[int(x)+int(y) for x,y in zip(geno1sum,geno1)]
                
        for i in range(0,len(variantList2)):
            (locus2,geno2)=variantList2[i]
            if geno2sum==[]:
                geno2sum=[int(x) for x in geno2]
            else:
                geno2sum=[int(x)+int(y) for x,y in zip(geno2sum,geno2)]
                
        genosum=[x+y for x,y in zip(geno1sum,geno2sum)]
        sumCase=sum([int(x>0) for x in genosum[0:indexControl]])
        sumControl=sum([int(x>0) for x in genosum[indexControl:len(genosum)]])
        score=sumCase-sumControl
        #if sumControl>0:
        #    score=0
        
    return (k,((gene1,gene2),score,sumCase,sumControl))

def getGene(variantData):
    variantData=variantData.split(" ")
    variantGene=variantData.pop(0)
    
    return (variantGene)

def createPairsGenes(k,v,genes):
    return [(str(sorted([k,gene])),(sorted([k,gene]),v)) for gene in genes]


# In[25]:

start_time = time.time()


if scope=='monogenic':
    if scale=='variant':
        scores=sc.textFile(filename).repartition(nPartitions).map(scoreVariantUnivariate).filter(lambda (k,(v1,v2,v3)): v1>=0).takeOrdered(1000, key=lambda (k,(v1,v2,v3)): -v1)

    if scale=='gene':
        scores=sc.textFile(filename,nPartitions).map(splitValues).groupByKey().map(lambda (k,v):scoreGeneUnivariate(k,v)).filter(lambda (k,(v1,v2,v3)): v1>0).takeOrdered(1000, key=lambda (k,(v1,v2,v3)): -v1)
    
if scope=='digenic':
    genes=sc.textFile(filename,nPartitions).map(getGene).distinct().takeOrdered(1000000)#.flatMap(lambda (k,v):scoreCompound(k,v)).takeOrdered(100000, key=lambda (k,v1,v2,v3): -v1)
    scores=sc.textFile(filename,nPartitions).map(splitValues).groupByKey().flatMap(lambda (k,v):createPairsGenes(k,v,genes)).groupByKey().map(lambda (k,v):scoreDigenicGene(k,v)).takeOrdered(1000, key=lambda (k,(genes,v1,v2,v3)): -v1)

end_time=time.time()
runtime=end_time - start_time
print(runtime)


# In[26]:

#genes=sc.textFile(filename,nPartitions).map(getGene).distinct().takeOrdered(1000000)
#scores[0:10]


# In[17]:

scores=[analysisName,scale,scope,start_time,end_time,runtime,scores]

with open("analyses/"+analysisName+'.txt', 'w') as outfile:
    json.dump(scores, outfile)
    


# In[ ]:

sc.stop()

