---
layout: page
title: RNA-seq running Deseq2
order: 102
session: 1
length: 10
toc: true
adapted: true
---
# Pre-filtering

In certain cases it might be a good idea to pre-filter your data. As there will be a few genes with very low read counts across the board in all groups.


 With filtering you want to try a balance of removing genes that have very low counts for all groups but you want a balance and not be so agressive that you remove useful genes. Typically there are two main ways of filtering at this step you can either, do it via 1 CPM(count per million) value or via a threshold on your data. The aim of pre-filtering is to reduce the reduce the number of statistical test being performed and to reduce the noise that can occur when there are many low genes with low reads.

We will begin by getting a summary of the rawCounts. Here we will be able to observe general figures such as means and interquartile ranges for our data.  
```
summary(rawCounts)
```
By looking at the summary of the rawCounts we can choose 5 as a create number. We can choose how many genes would be left if we filtered by 5 counts, by using the ```dim``` command. This command checks the dimension of matrix or table.
```
dim(deseq2Data[rowSums(counts(deseq2Data)) > 5, ])
```
We can observe that 5 isn't that aggresive of a prunning so so we apply the filtering to the deseq2Data object we previously created.
```
deseq2Data <- deseq2Data[rowSums(counts(deseq2Data)) > 5, ]
```
Now we run DESeq this should take a minute to perform but go ahead and use the ```head``` command to see the new object. 
```
deseq2Data <- DESeq(deseq2Data)
```
