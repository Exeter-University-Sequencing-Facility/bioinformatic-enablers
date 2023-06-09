---
layout: page
title: RNA-seq Introduction
order: 100
session: 1
length: 10
toc: true
adapted: true
---

## Welcome to the RNA-seq portion of this course

By the end of this session you should be able to:
  * Manipulate tables into suitable formats for Deseq2.
  * Perform differential gene expression with Deseq2.
  * Generate some visual plots.
  * Perform Gene Set enrichment Analysis with Clusterprofiler4.


But before we start with any of the fun stuff it is important to understand the experimental design especially if you didn't perform the experiment ourselves. The paper for the dataset we will be working is [paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5528589/ "paper")

The primary objective of the study was to identify key genes associated with colorectal cancer agressiveness. The RNA-seq has 54 samples from 3 tissues:primary colorectal, liver metastases and normal colon belonging to 18 patients. In this practical we are going to look at the differentially expressed genes between the primary colorectal tissue and the normal colon tissues.
