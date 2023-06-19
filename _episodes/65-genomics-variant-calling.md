---
layout: page
title: Variant Detection.
order: 65
session: 2
length: 15
toc: true
---

# Variant Detection

Remember we are following the Broad Institute's best practice for [Variant Calling](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-)

## Preparing the BAM file

There are a few steps that 'should be' performed before calling variants.

- Mark Duplicates
- Left Align Indels
- Recalibrate QScores


Make sure you are in a compute node and not the login node.
Make sure you are in the correct folder.
```
cd /lustre/projects/Research_Project-BioTraining/ecr2023/${USER}
```

Activate the gatk conda environment

```
. "/gpfs/ts0/shared/software/Miniconda3/4.9.2/etc/profile.d/conda.sh"
conda activate /lustre/projects/Research_Project-BioTraining/ecr2023/bioconda-envs/gatk4
```

### Mark Duplicates

During the preparation for the sample (library preparation) and the sequence itself, duplication of molecules may occur. This process finds such duplicates in the BAM file and marks them as such. Later steps in can then use this information to give the read less weight than a unique molecule when identifying variants.

Run the gatk module `MarkDuplicates`

```
input_folder=31_bwa_aligned
output_folder=32_mark_duplicates
sample=wildtype
mkdir -p ${output_folder}

gatk MarkDuplicates \
    --INPUT=${input_folder}/${sample}_sorted.bam  \
    --OUTPUT=${output_folder}/${sample}_sorted_md.bam  \
    --METRICS_FILE=${output_folder}/${sample}_md.metrics;
```

### Left Align Indels

This is strictly not necessary as latter steps have this build in. I have left

This is the process of ensuring that alignments with more than one possible solution are treated consistently. I will explain with an example.

If xyour reference sequence is GACACACACG
and the sequence in your sample is GACACACG, then there are several equally valid alignments

```{txt}
GACACACACG  
G--ACACACG

GACACACACG
GACACAC--G
```

etc.

The nature of the aligner means that it will not always pick the same alignment solution - previous versions of the variant caller might generate mutliple calls for the same variant if this was not done.

This is strictly not necessary as latter steps have this build in. I have left this is as it is still useful for visualisation.

```
input_folder=32_mark_duplicates
output_folder=33_align_indels
mkdir -p ${output_folder}
reference_fasta=ncbi_dataset/data/GCF_000011545.1/GCF_000011545.1_ASM1154v1_genomic.fna
sample=wildtype

gatk LeftAlignIndels \
    --reference=${reference_fasta} \
    --input=${input_folder}/${sample}_sorted_md.bam  \
    --OUTPUT=${output_folder}/${sample}_sorted_md_lii.bam

gatk BuildBamIndex \
    --INPUT=${output_folder}/${sample}_sorted_md_lii.bam
```

## Recalibrate Basecalling

This step adjusts the Illumina derived basecalls based on real error rates rather than confidence levels.

When the Illumina software generates a basecall, it is out of any context - it just reports how confident it is that the basecall is correct (relative to the second best basecall).

After aligment it is possible to refine this estimate based on real error rates. e.g. when Illumina give a Qscore of 30 expects to be wrong 0.1% of the time. After alignment that might turn out to be 0.15% on average and the QScore can be adjusted accordingly.

To perform this simply a database of expected SNP is required.
We will skip this step. More details [here](https://gatk.broadinstitute.org/hc/en-us/articles/360035890531)

### Now the loop for the rest of the sample

```
while read sample; do 
    echo ${sample}
    input_folder=31_bwa_aligned
    output_folder=32_mark_duplicates
    mkdir -p ${output_folder}

    gatk MarkDuplicates \
        --INPUT=${input_folder}/${sample}_sorted.bam  \
        --OUTPUT=${output_folder}/${sample}_sorted_md.bam  \
        --METRICS_FILE=${output_folder}/${sample}_md.metrics;

    input_folder=32_mark_duplicates
    output_folder=33_align_indels
    mkdir -p ${output_folder}

    gatk LeftAlignIndels \
        --reference=${reference_fasta} \
        --input=${input_folder}/${sample}_sorted_md.bam  \
        --OUTPUT=${output_folder}/${sample}_sorted_md_lii.bam

    gatk BuildBamIndex \
        --INPUT=${output_folder}/${sample}_sorted_md_lii.bam
done < samples.txt
```

## Variant Calling

Finally we get to the steps where variants are identified.
It is a 2 step process - first each sample is 'genotyped' - an assessment is made of the likelihood of a variation at each possible location.

Secondly the variants are called from the joint information. This avoid the issue with some pipelines where a variant in  two samples falls just above and just below and arbitary cut-off.

```
ploidy=1
input_folder=33_align_indels
output_folder=34_call_haplotypes
mkdir -p ${output_folder}
sample=wildtype

gatk HaplotypeCaller \
    --input=${input_folder}/${sample}_sorted_md_lii.bam \
    --reference=${reference_fasta} \
    --output=${output_folder}/${sample}_g.vcf \
    --emit-ref-confidence GVCF \
    --sample-ploidy ${ploidy} \
    --bam-output=${output_folder}/${sample}_hc.bam \
    --native-pair-hmm-threads 8
```

## Loop

```
while read sample; do 
echo;echo;echo
echo ${sample}
sleep 5
gatk HaplotypeCaller \
    --input=${input_folder}/${sample}_sorted_md_lii.bam \
    --reference=${reference_fasta} \
    --output=${output_folder}/${sample}_g.vcf \
    --emit-ref-confidence GVCF \
    --sample-ploidy ${ploidy} \
    --bam-output=${output_folder}/${sample}_hc.bam \
    --native-pair-hmm-threads 8
done < samples.txt
```

## Joint Variant Calling

```
input_folder=34_call_haplotypes
ouput_folder=34_call_haplotypes
reference_fasta=ncbi_dataset/data/GCF_000011545.1/GCF_000011545.1_ASM1154v1_genomic.fna
sample=wildtype
ploidy=1

unset samples_command
while read sample; do
  samples_command="${samples_command} --variant ${input_folder}/${sample}_g.vcf"
done < samples.txt
echo ${samples_command}

gatk CombineGVCFs \
  --reference=${reference_fasta} \
  --output=${output_folder}/combined_samples.vcf \
  ${samples_command};

gatk GenotypeGVCFs \
  --sample-ploidy ${ploidy} \
  --reference=${reference_fasta} \
  --variant=${output_folder}/combined_samples.vcf \
  --output=${output_folder}/genotyped_samples.vcf ;

gatk VariantsToTable \
  --variant=${output_folder}/genotyped_samples.vcf \
  --output=${output_folder}/genotyped_samples.tsv \
  --fields CHROM \
  --fields POS \
  --fields ID \
  --fields REF \
  --fields ALT \
  --fields QUAL \
  --fields FILTER \
  --genotype-fields GT \
  --genotype-fields AD \
  --genotype-fields DP \
  --genotype-fields GQ \
  --genotype-fields PL;
```

## Add some annotation

```
cd /lustre/projects/Research_Project-BioTraining/ecr2023/${USER}/ncbi_dataset/data/GCF_000011545.1
conda activate /lustre/projects/Research_Project-BioTraining/ecr2023/bioconda-envs/utils/
gt gff3 -sortlines -tidy -retainids genomic.gff > genomic.sorted.gff
bgzip genomic.sorted.gff
tabix genomic.sorted.gff.gz
ls -latr
```

### Coverage by Gene

```
mkdir -p 37_coverage_by_gene
reference_gff=ncbi_dataset/data/GCF_000011545.1/genomic.gff
reference_fasta=ncbi_dataset/data/GCF_000011545.1/GCF_000011545.1_ASM1154v1_genomic.fna
input_folder=33_align_indels

while read sample; do
   echo ${sample}
   bedtools coverage \
    -a ${reference_gff} \
    -b ${input_folder}/${sample}_sorted_md_lii.bam \
    | awk '$3=="CDS"' \
    | sort -t$'\t' -k13 \
    > 37_coverage_by_gene/${sample}_cds_coverage.bed
done < samples.txt
```

## ' fish out incomplete genes'

```
for x in *.bed; do echo $x; perl -F"\t" -alne 'print $_ if $F[12].to_f < 1.0' $x; done
```

### Annotate SNPs

```
input_folder=34_call_haplotypes
l bedtools intersect -loj  -a ${input_folder}/genotyped_samples.vcf -b <( grep -P "\tCDS\t"  ${reference_gff}) > ${input_folder}/genotyped_samples.annotated.vcf
```
