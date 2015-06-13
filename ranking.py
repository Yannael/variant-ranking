
# coding: utf-8

# In[2]:

from pyspark import SparkContext
sc = SparkContext("local", "Simple App")


# In[101]:

snpsMat="/Users/yalb/snpsMat.txt"
snpsMat = sc.textFile(snpsMat)
snpsMat.count()


# In[9]:

pedigree="/Users/yalb/pedigree.txt"
pedigree = sc.textFile(pedigree)
pedigree.count()
pedigree.first()


# In[182]:

def scoreDeNovo(variantData):
    #variantData=snpsMat.first()
    variantData=variantData.split(" ")
    variantLocus=variantData.pop(0)
    
    score=0
    for i in range(0,11):
        score=score+int(variantData[(i*3)]=="1" and variantData[i*3+1]=="0" and variantData[i*3+2]=="0")
    return (variantLocus,score)


# In[283]:

snpsMat="/Users/yalb/snpsMat.txt"
snpsMat = sc.textFile(snpsMat)
snpsMat=snpsMat.collect()
snpsMat=sc.parallelize(snpsMat,8)
snpsMat=snpsMat.repartition(8)


# In[226]:

scores=snpsMat.map(scoreDeNovo).filter(lambda x: x[1]>0)


# In[282]:

snpsMat.take(1)


# In[228]:

def getKeyValue(variantData):
    variantData=variantData.split(" ")
    variantLocus=variantData.pop(0)
    (gene,locus)=variantLocus.split("_")
    return(gene,(locus,variantData))


# In[274]:

def getScorePair(geno1,geno2):
    score=0
    listTrios=[]
    for i in range(0,11):
        scoreTrio=int(geno1[i*3]=="1" and geno2[i*3]=="1" and 
                        ((geno1[i*3+1]=="1" and geno2[i*3+1]=="0" and geno1[i*3+2]=="0" and geno2[i*3+2]=="1") or
                        (geno1[i*3+1]=="0" and geno2[i*3+1]=="1" and geno1[i*3+2]=="1" and geno2[i*3+2]=="0")))
        if scoreTrio>0:
            listTrios.append(i)
            score=score+scoreTrio
    return (score,listTrios)


# In[277]:

def scoreCompound(variantList):
    variantList=list(variantList)
    result=[]
    for i in range(0,len(variantList)-1):
        (locus1,geno1)=variantList[i]
        for j in range(i+1,len(variantList)):
            (locus2,geno2)=variantList[j]
            (score,listTrios)=getScorePair(geno1,geno2)
            if score>0:
                result.append((locus1,locus2,score,listTrios))
        
    return result


# In[284]:

snpsMat.map(getKeyValue).groupByKey().mapValues(scoreCompound).filter(lambda (k,v):len(v)>0).collect()


# In[10]:

hbase_rdd.first()


# In[ ]:



