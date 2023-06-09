---
layout: page
title: RNA-seq Data preparation
order: 101
session: 1
length: 10
toc: true
adapted: true
---
## Preparation for the data

For this practical we will start with raw unnormalised matrix counts. You can seek advice from us or others regarding the on the aligning of transcripts and the quantification of them. The quantification data required to run differential expression analysis using DEseq2 are raw readcounts for either genes or transcripts. We will use the output from HTseq as a starting point.

We will begin by downloading the starting the data and assigning it to a variable.Variables are objects in R that have values and in this course we will create many variables with objects such as characters or even matrices. You can create your own varibles just make sure that they don't start with a number and that they are informative.

```
rawCounts <- read.delim("http://genomedata.org/gen-viz-workshop/intro_to_deseq2/tutorial/E-GEOD-50760-raw-counts.tsv")


# Read in the sample mappings
sampleData <- read.delim("http://genomedata.org/gen-viz-workshop/intro_to_deseq2/tutorial/E-GEOD-50760-experiment-design.tsv")

```

```
head(rawCounts)
head(sampleData)
```


#got data from https://genviz.org/module-04-expression/0004/02/01/DifferentialExpression/

#Now we have our raw data but we need to ensure that it is in the right format fro Deseq2.
#For Deseq2 the rows much be gene names and the columns must be sample_ids. It is important to ensure that Deseq2 doesn't allow sample names that start with numerical values.

```
#extract the gene symbol from the raw read counts we imported.
geneID <- rawCounts$Gene.ID
#it is never a bad idea to have a look at the code you wrote to ensure that teh output is what you expected. You can have a lok at the fist few entries of a datatype by using head or the last few rows of a data type by using tail.
head(geneID)
#This creates a boonlean table of columns that start with SRR which are our sample names.
sampleIndex <- grepl("SRR\\d+", colnames(rawCounts))
#Sanity check
head(sampleIndex)

#creates a matrix that is a type of data with the columns matching sampleIndex=TRUE
rawCounts <- as.matrix(rawCounts[,sampleIndex])
head(rawCounts)
#assigns the gene symbol as the rownames of the data.
rownames(rawCounts) <- geneID
head(rawCounts)
#now your data read counts are ready for deseq2.



head(sampleData)
#creates the rownames as the sample_name which are under the Run folder.
rownames(sampleData) <- sampleData$Run
keep <- c("Sample.Characteristic.biopsy.site.", "Sample.Characteristic.individual.")
#creating a table with the columns that we want such as Sample.Characteristic.biopsy.site.", "Sample.Characteristic.individual.
sampleData <- sampleData[,keep]
#we rename the column names as tissueType adn IndividualId so that they are more informative.
colnames(sampleData) <- c("tissueType", "individualID")
head(sampleData)
#We turned the IndividualsID into factors. Factors is a data structure used to for the storage of categorical data. Which is required for Deseq2.
sampleData$individualID <- factor(sampleData$individualID)
head(sampleData)

#rawCounts <- rawCounts[,unique(rownames(sampleData))] is it necessary?

#this checks that the column names of the rawConts match with the rownames of the SampleData. if they are not in the right order Deseq2 won't work.
all(colnames(rawCounts) == rownames(sampleData))

# rename the tissue types

#this is  created a function called rename_tissue that thanks in x and changes words such as normal-looking surrounding colonic epithelium to normal and maks it easier for us to deal with.  

rename_tissues <- function(x){
  x <- switch(as.character(x), "normal"="normal_looking_surrounding_colonic_epithelium", "primary tumor"="primary_colorectal_cancer",  "colorectal cancer metastatic in the liver"="metastatic_colorectal_cancer_to_the_liver")
  return(x)
}
sampleData$tissueType <- unlist(lapply(sampleData$tissueType, rename_tissues))

# Order the tissue types so that it is sensible and make sure the control sample is first: normal sample -> primary tumor -> metastatic tumor
sampleData$tissueType <- factor(sampleData$tissueType, levels=c("normal_looking_surrounding_colonic_epithelium", "primary_colorectal_cancer", "metastatic_colorectal_cancer_to_the_liver"))

# Create the DEseq2DataSet object
deseq2Data <- DESeqDataSetFromMatrix(countData=rawCounts, colData=sampleData, design= ~ individualID + tissueType)
#Sanity check to understand the number of rows.
dim(deseq2Data)

```
