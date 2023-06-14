---
layout: page
title: Alignment against reference sequence.
order: 64
session: 2
length: 15
toc: true
---

## Genomic Alignment

Here we will use a program called [BWA](https://github.com/lh3/bwa) to align your reads against the reference genome.  
The program will examine each read in turn and work out from which part of the genome it originates. 

The result is a file called a 'Sequence Alignment Map' or SAM file. More frequenctly a Binary (BAM) version is used see the [wiki page](https://en.wikipedia.org/wiki/SAM_(file_format)) for details of the format.



```
srun  --export=ALL -D . -p bioseq  --time=12:00:00 -A Research_Project-BioTraining --nodes=1  --ntasks-per-node=8 --pty bash -i
```

Make sure you are in the correct folder
```
cd /lustre/projects/Research_Project-BioTraining/ecr2023/${USER}
```

```
. "/gpfs/ts0/shared/software/Miniconda3/4.9.2/etc/profile.d/conda.sh"
conda activate /lustre/projects/Research_Project-BioTraining/ecr2023/bioconda-envs/bwa
```

```
sample=wildtype
input_folder=11_trimmed_reads
output_folder=31_bwa_aligned
reference_fasta=ncbi_dataset/data/GCF_000011545.1/GCF_000011545.1_ASM1154v1_genomic.fna
mkdir -p ${output_folder}
```


```
time bwa mem \
    -t 8 \
    -R "@RG\tID:${sample}\tSM:${sample}\tLB:l1\tPL:ILLUMINA\tDS:MiSeq" \
    ${reference_fasta} \
    ${input_folder}/${sample}_R1_fastp.fastq.gz \
    ${input_folder}/${sample}_R2_fastp.fastq.gz > \
    ${output_folder}/${sample}_temp.sam \
```

### Quick check
```
samtools flagstat ${output_folder}/${sample}_temp.sam
```
Look for a high mapping rate.


### Convert to BAM and sort
```
samtools view -bSh \
    ${output_folder}/${sample}_temp.sam > \
    ${output_folder}/${sample}_temp.bam 
```
```
time samtools sort ${output_folder}/${sample}_temp.bam  \
    --threads 8 \
    -T ${output_folder}/${sample}_unsorted.bam  \
    -o ${output_folder}/${sample}_sorted.bam
```
```
samtools flagstat ${output_folder}/${sample}_sorted.bam
```
This should be the same as before. Notice that the BAM file is much smaller than the SAM formatted one.

### Run Qualimap

Qualimap gives us some useful data on alignment rates ond coverage of the genome

```
qualimap_folder=${QC_FOLDER}/qualimap
reference_gff=ncbi_dataset/data/GCF_000011545.1/genomic.gff
reference_fasta=ncbi_dataset/data/GCF_000011545.1/GCF_000011545.1_ASM1154v1_genomic.fna

qualimap bamqc \
    --java-mem-size=4G \
    -bam ${output_folder}/${sample}_sorted.bam \
    -gff ${reference_gff} \
    -outdir ${qualimap_folder}/${sample}
```


## Now execute in a loop
```
while read sample; do 
    echo $sample
    bwa mem \
        -t 8 \
        -R "@RG\tID:${sample}\tSM:${sample}\tLB:l1\tPL:ILLUMINA\tDS:MiSeq" \
        ${reference_fasta} \
        ${input_folder}/${sample}_R1_fastp.fastq.gz \
        ${input_folder}/${sample}_R2_fastp.fastq.gz > \
        ${output_folder}/${sample}_temp.sam \
        
    samtools view -bSh \
        ${output_folder}/${sample}_temp.sam > \
        ${output_folder}/${sample}_temp.bam 
    samtools sort ${output_folder}/${sample}_temp.bam  \
        --threads 8 \
        -T ${output_folder}/${sample}_unsorted.bam  \
        -o ${output_folder}/${sample}_sorted.bam
    samtools flagstat ${output_folder}/${sample}_sorted.bam > ${output_folder}/${sample}_aligment_stats.txt
    qualimap bamqc \
        --java-mem-size=4G \
        -bam ${output_folder}/${sample}_sorted.bam \
        -gff ${reference_gff} \
        -outdir ${qualimap_folder}/${sample}
done < samples.txt
```




Run multiqc again - changing the output filename and have a look at the results.

## Summary

We now have aligned and sorted files for ongoing analysis.

