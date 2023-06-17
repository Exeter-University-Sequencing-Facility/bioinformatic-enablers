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


 With filtering you want to try a balance of removing genes that have very low counts for all groups but not be so aggressive that you remove useful genes. Typically there are many ways of filtering at this step you can either, do it via 1 CPM(count per million) value or via a threshold on your data. The aim of pre-filtering is to reduce the reduce the number of statistical test being performed and to reduce the noise that can occur when there are many low genes with low counts. There is not a definitely way of filtering data at this stage and it depends on your experiment.

We know that the lowest number of reads per sample was 5 million pair end reads thus we know that 1 count per million reds is 5. We use five as our cut off. We can use dim again to check the number of rows in our new table when we remove anything that has less than 5 counts.

```
dim(deseq2Data[rowSums(counts(deseq2Data)) > 5, ])
```
How does this compare to the previous total?

<details>
           <summary>What percentage of rows have been retained?</summary>
           <p>53.9%</p>
         </details>

The important thing with this step is to be aware that we have lost around 47% of our genes and if things don't look right later one this is where we can try to finetune. As we take this aspects into consideration we are happy to apply the filtering to the deseq2Data object we previously created.
```
deseq2Data <- deseq2Data[rowSums(counts(deseq2Data)) > 5, ]
dim(deseq2Data)
```
Now we run DESeq this should take a minute to perform but go head and use the ```head``` command to see the new object.
```
deseq2Data <- DESeq(deseq2Data)
```

Now before we go to the next step it is important to
