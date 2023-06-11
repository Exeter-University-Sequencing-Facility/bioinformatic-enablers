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

We will begin by downloading the starting the data and assigning it to a variable.Variables are objects in R that have values and in this course we will create many variables with objects such as characters or even matrices. You can create your own varibles just make sure that they don't start with a number and that they are informative.We reccomnd it is best to copy the code as you would learn and understand it better if you take the time to type it out.

```
rawCounts <- read.delim("http://genomedata.org/gen-viz-workshop/intro_to_deseq2/tutorial/E-GEOD-50760-raw-counts.tsv")


# Read in the sample mappings
sampleData <- read.delim("http://genomedata.org/gen-viz-workshop/intro_to_deseq2/tutorial/E-GEOD-50760-experiment-design.tsv")

```
Throughout this tutorial we will use the command ```head``` that displays the first few rows of a table. It is always good to check and look what is going on in between steps and that a command has performed in the way that you expected it to perform.

```
head(rawCounts)
head(sampleData)
```
As you can see with the rawCounts there is a column corresponsind to ENSEMBL gene names another one corresponding to gene symbols while the rest are of the columns are the sample names. For the sample Data we have a lot of information regarding links and the type of tissue in this exercise we are interested in just a few columns. It is important to ensure that we extract the needed information.


Now we have our raw data but we need to ensure that it is in the right format for Deseq2. Deseq2 is a method to detect differentially expressed genes it uses various algorithms to calculate it and it takes into consideration the sequencing depths. [Here](http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#standard-workflow) is a link with an example satndard workflow for Deseq2. For Deseq2 the rows much be gene names and the columns must be sample names. It is important to ensure that for Deseq2  that the sample names do not start with numbers.

The gene ids were extracted and assigned to the variable called geneID.
```
#extract the gene symbol from the raw read counts we imported.
geneID <- rawCounts$Gene.ID
head(geneID)
```

Now we are extracting the sample names, not their content just their names by using the command ```grepl``` that find matching patterns.

```
sampleIndex <- grepl("SRR\\d+", colnames(rawCounts))
#Sanity check
head(sampleIndex)
```
Once we have checked we now create a matrix extracting the columns corresponding to the sampleIndex.
```
rawCounts <- as.matrix(rawCounts[,sampleIndex])
head(rawCounts)
```
We are almost done we now generate rownames for our new table using the geneIDs.
```
rownames(rawCounts) <- geneID
head(rawCounts)
```
Now your data RAW counts are ready for Deseq2.We can use the SampleData variable to create the suitable metadata needed for Deseq2. The metadata highlights what group each data is assigned to.

```
head(sampleData)
```
You can see from the SampleData that the Run columns corresponds to the sample name. What we want are the columns Run,Sample.Characteristic.biopsy.site. and Sample.Characteristic.individual. In this particular example we are going to use both the biopsy site and the patient details. In theory depending on what you are looking for you could just use the biopsy type but by including the Sample.Characteristic.individual. you can use this information for other experiments when there might be a batch effect going on.

The code below will create the sample names as rownames.
```
rownames(sampleData) <- sampleData$Run

```
Now we are trying to retaain the columns that we want such as Sample.Characteristic.biopsy.site.", "Sample.Characteristic.individual.
```
 keep <-c("Sample.Characteristic.biopsy.site.","Sample.Characteristic.individual.")
```
We are now adding those two columns to the sampleData.
```
sampleData <- sampleData[,keep]
```
Perfect now we have a table with all the information we need. A few more touches that need to be done are to rename the column names as tissueType adn IndividualId so that they are more informative for our benefits.
```
colnames(sampleData) <- c("tissueType", "individualID")
head(sampleData)
```
We turned the IndividualID and the Tissuetype into factors. Factors is a data structure used to for the storage of categorical data. Which is required for Deseq2.
```
sampleData$individualID <- factor(sampleData$individualID)
sampleData$tissueType <- factor(sampleData$tissueType)
head(sampleData)
```

This checks that the column names of the rawCounts match with the rownames of the SampleData. if they are not in the right order Deseq2 won't work.
```
all(colnames(rawCounts) == rownames(sampleData))
```
Now we plan to rename the tissue types
by creating a function called rename_tissue that takes in x and changes words such as normal into normal_looking_surrounding_colonic_epithelium and the description has more information. It is also important to not have spaces in names.
```
rename_tissues <- function(x){
  x <- switch(as.character(x), "normal"="normal_looking_surrounding_colonic_epithelium", "primary tumor"="primary_colorectal_cancer",  "colorectal cancer metastatic in the liver"="metastatic_colorectal_cancer_to_the_liver")
  return(x)
}
sampleData$tissueType <- unlist(lapply(sampleData$tissueType, rename_tissues))
```

Order the tissue types so that it is sensible and make sure the control sample is first: normal sample -> primary tumor -> metastatic tumor
```
sampleData$tissueType <- factor(sampleData$tissueType, levels=c("normal_looking_surrounding_colonic_epithelium", "primary_colorectal_cancer", "metastatic_colorectal_cancer_to_the_liver"))
```
Once that is sorted we create the DEseq2DataSet object
```
deseq2Data <- DESeqDataSetFromMatrix(countData=rawCounts, colData=sampleData, design= ~ individualID + tissueType)
```
Now we are doing a sanity check to check the number of rows.
```
dim(deseq2Data)

```
Some bits of this workshop were adapted from https://genviz.org/module-04-expression/0004/02/01/DifferentialExpression/
