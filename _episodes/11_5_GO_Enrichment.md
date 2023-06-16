---
layout: page
title: RNA-seq Functional Analysis
order: 104
session: 1
length: 10
toc: true
adapted: true
---

# Functional Analysis

At this stage you will have a list of differentially expressed genes and you will know if they are significant or not but if you want to know if certain processes are suppressed or activated. For that you would to know the processes that those genes are involved in. For this we are going to use to use Clusterprofiler. Clusterprofiler is a tool that is used for enrichment analysis. One can perform gene set enrichement analysis(gsea)or overall representation analysis.(ora). You can read in more detail the different between the two methods over [here](https://yulab-smu.top/biomedical-knowledge-mining-book/enrichment-overview.html). We are going to use gsea in this example as it is tends to be more precise when looking at subtler and smaller differences. The first step for gene set enrichment analysis is to rank our data. There are multiple ways to rank it, some sort it by descending order of log2foldChange. In this case we are sorting by log2foldchange and by p-adjusted value.

#As an exercise you can try other ways of sorting this data and see if it yields other results.

```
sorted_result_primary_vs_normal<-sorted_result_primary_vs_normal %>% group_by(round(desc(log2FoldChange))) %>% arrange(padj, .by_group = TRUE)
head(sorted_result_primary_vs_normal)
na_omited_sorted_result_primary_vs_normal<-na.omit(sorted_result_primary_vs_normal) #na values have to be removed as clusterprofiler won't account for them.
dim(na_omited_sorted_result_primary_vs_normal) #check demensions
dim(sorted_result_primary_vs_normal) #check demensions

original_gene_list <-na_omited_sorted_result_primary_vs_normal$`round(desc(log2FoldChange))`
```
We are using the rounded log2foldchanges so that it is still descending clusterprofiler doesn't like it if it is not in a descending value.

```
original_gene_list<-(-original_gene_list) #change the negative signs to positive and vice-versa.
names(original_gene_list) <-na_omited_sorted_result_primary_vs_normal$genes
head(original_gene_list)
```
If we look at the original_gene_list we should see it list of genes with its rounded log fold change. For the next step we are going to use the org.Hs.eg.db. That is a complex R object that contains genome wide annotation for humans. For very well annotated and model organisms there should be libraries like these that can you install and use such as org.Dm.eg.db.
```
library(org.Hs.eg.db)
gbp <- gseGO(original_gene_list, minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.25, pAdjustMethod = "BH", OrgDb=org.Hs.eg.db,seed = 123,ont="BP",keyType= 'ENSEMBL')

p_dot <- dotplot(gbp, showCategory=200, split=".sign", title = "Dotplot - all") + facet_grid(.~.sign)

p_dpt
```
If we want to save the plot you can save using Export > Save as Image >
remember to give it an informative name. We recommend to save it with a width of 1000 px and a height of 5000 px so that you can visualise the dot plot.

When comparing colorectal cancer tumour cells and their normal counterparts are the GO terms represented what you expected? Are there any surprises?  
```
edox <- setReadable(gbp, 'org.Hs.eg.db', 'ENSEMBL') # reads the ENSEMBL values into edox from the gseGo we developed.

p3 <- cnetplot(edox, foldChange=original_gene_list, circular = TRUE, colorEdge = TRUE)

p3
```

# Finishing touches

It is important to save what you did. To save anything important such as figure and etc but also to save files and scripts. 

#reference

https://www.bioinformatics.babraham.ac.uk/projects/seqmonk/Help/3%20Visualisation/3.2%20Figures%20and%20Graphs/3.2.25%20The%20PCA%20Plot.html#:~:text=The%20PCA%20plot%20is%20a,which%20best%20separate%20your%20data.
