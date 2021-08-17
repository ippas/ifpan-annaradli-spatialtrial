# ifpan-annaradli-spatialtrial

<h2>Description of the project</h2>
<p>The purpose of the project was to determine gene expression changes in different forebrain tissue domains following acute administration of risperidone (<em>i.p.</em>). </p> <p>The details of the methodology are listed below:</p> * 8 adult males (C57BL/6) were housed in 2 cages. (4 mice per cage) <br /> 
* Mice 1-2, 1-4, 2-2 and 2-4 were given risperidone. <br />* Mice 1-1, 1-3, 2-1 and 2-3 were given saline. <br /> * 2 hours after the injection mouse were killed, their brains were extracted, embedded in OCT and frozen in a isopentane bath (isopentane was placed in liquid nitrogen). <br /> * The brains were sectioned (10 um) in a cryostat (CT: -20C, OT: -10C).<br /> * 4 sections of the forebrains (2 per each group) were placed on a Visium Gene Expression slide (10x Genomics). <br /> * The sections were stained with H&E and then used for cDNA synthesis and subsequent library construction. <br /> * The libraries were sequenced in CeGaT, Germany (NovaSeq 6000, PE, Read1:28, Indexi5:10, Indexi7:10, Read2:90). Converting to fastq format and demultiplexing was performed by CeGaT.<br /> 

### samples summary

|Sample ID| Visium Slide Area ID | Mouse ID | Group | Double TT Index Well ID|
|---------| -------------------- | -------- | ------| -----------------------|
|A| A1                   | 1-2      | risperidone | A1 |
|B|B1| 2-2 | risperidone | B1|
|C|C1| 2-3 | saline | C1|
|D|D1| 1-1 | saline | D1 |

### RNA-seq samples summary
|Sample ID|Filename|Sequence type|FastQC failed parameters|
|---|---|---|---|
|A|S3647Nr1.1|Read 1|Sequence Duplication Levels|
|A|S3647Nr1.2|Read 2|Sequence Duplication Levels, Per sequence GC content|
|B|S3647Nr2.1|Read 1|Sequence Duplication Levels|
|B|S3647Nr2.2|Read 2|Sequence Duplication Levels, Per sequence GC content, Overrepresented sequences|
|C|S3647Nr3.1|Read 1|Sequence Duplication Levels|
|C|S3647Nr3.2|Read 2|Sequence Duplication Levels, Per sequence GC content, Overrepresented sequences|
|D|S3647Nr4.1|Read 1|Sequence Duplication Levels|
|D|S3647Nr4.2|Read 2|Sequence Duplication Levels, Per sequence GC content, Overrepresented sequences|

<h2>Labnotes</h2>
Downloaded RNA-seq data. <br />
Performed fastQC and md5sum checking.<br />

### 2021-06-07
Downloaded spaceranger and reference mouse genome from 10x Genomics webiste.<br />

### 2021-06-08
Uploaded the data, spaceranger and reference genome to server. md5sum checked on spaceranger and reference genome mm10 (compliant). Unpacked (`tar`) both genome and spaceranger to opt/. Prepended spaceranger to $PATH on the server. Ran `spaceranger sitecheck` and `spaceranger testrun` successfully. - <br />

### 2021-06-09
md5sum checked on data files (.fastq.gz) uploaded to the server. Renamed the fastq files in the following (Illumina) [convention](https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/NamingConvention_FASTQ-files-swBS.htm): <br /> 
S3647Nr1.1.fastq.gz -> S3647Nr1_S1_L001_R1_001.fastq.gz <br />
S3647Nr1.2.fastq.gz -> S3647Nr1_S1_L001_R2_001.fastq.gz <br />
S3647Nr2.1.fastq.gz -> S3647Nr2_S2_L001_R1_001.fastq.gz <br />
S3647Nr2.2.fastq.gz -> S3647Nr2_S2_L001_R2_001.fastq.gz <br />
S3647Nr3.1.fastq.gz -> S3647Nr3_S3_L001_R1_001.fastq.gz <br />
S3647Nr3.2.fastq.gz -> S3647Nr3_S3_L001_R2_001.fastq.gz <br />
S3647Nr4.1.fastq.gz -> S3647Nr4_S4_L001_R1_001.fastq.gz <br />
S3647Nr4.2.fastq.gz -> S3647Nr4_S4_L001_R2_001.fastq.gz <br />

Prepended spaceranger again. `spaceranger sitecheck` ran successfully. 
Ran first analysis on S3647Nr1 sample (A) - start 13:08, finish 14:19.
Ran analyses on rest of the samples.

### 2021-06-10
Prepended spaceranger again.`spaceranger sitecheck` ran successfully. 
Ran `spaceranger aggr`on 4 samples - start 14:26 , finish 14:41.

### Software 
* spaceranger v1.2.2
* mouse genome mm10 Reference - 2020-A 
* fastQC v0.11.8
* md5sum v8.28
* Loupe Browser v5.0/v5.1.0
* mrc: v4.0.2
* mrp: v4.0.2 
* Anaconda: numpy: 1.15.4, scipy: 1.1.0, pysam: 0.16.0.1, h5py: 2.8.0, pandas: 0.24.2
* STAR: 2.7.2a
* samtools 1.10
* htslib 1.10.2
* Martian Runtime v4.0.2


### Analysis Mateusz
Merging bam together:
```
samtools merge /tmp/merged_four_samples.bam results/runs/S3647Nr1/SPATIAL_RNA_COUNTER_CS/SPATIAL_RNA_COUNTER/_BASIC_SPATIAL_RNA_COUNTER/WRITE_POS_BAM/fork0/join-u5a05c0aed8/files/pos_sorted_single_sample.bam results/runs/S3647Nr2/SPATIAL_RNA_COUNTER_CS/SPATIAL_RNA_COUNTER/_BASIC_SPATIAL_RNA_COUNTER/WRITE_POS_BAM/fork0/join-u6bdec0ddab/files/pos_sorted_single_sample.bam results/runs/S3647Nr3/SPATIAL_RNA_COUNTER_CS/SPATIAL_RNA_COUNTER/_BASIC_SPATIAL_RNA_COUNTER/WRITE_POS_BAM/fork0/join-u6e48c0deab/files/pos_sorted_single_sample.bam results/runs/S3647Nr4/SPATIAL_RNA_COUNTER_CS/SPATIAL_RNA_COUNTER/_BASIC_SPATIAL_RNA_COUNTER/WRITE_POS_BAM/fork0/join-u70d8c0def4/files/pos_sorted_single_sample.bam
```

```
samtools sort /data/merged_four_samples.bam -o /DATA/merget_four_samples_sorted.bam
samtools index /data/merged_four_samples_sorted.bam

```

Performing an analysis in MACS3:
```
docker run -u 1003:1002 -v $PWD:/data/ ubuntu:macs3 macs3 callpeak -t /data/data/merged_four_samples_sorted.bam
-n merged_samples_sorted --outdir /data/data/macs3_results/sorted_bam
```


Annotate peaks:
```
./annotate_peaks.sh
```

Creating a gtf file to create a library
```

```


Creating a corrected library for [spaceranger count](https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/using/count)
```
spaceranger mkref --genome=oprm1_ref_without_strand --fasta=opt/refdata-gex-mm10-2020-A/fasta/genome.fa --genes=opt/refdata-gex-mm10-2020-A/genes/genes.gtf
```

```
peaks2gtf_gene.bed | tail +142 | awk '{sub("$", "+"$18)}; 1' | awk 'BEGIN{FS=OFS"\t"} {gsub(/-1\+/, "-" $18)} 1 {gsub(/1\+/, "+" $18)} 1' | awk '{print $0, $6==$17}' | sed 's/ /\t/' | cut -f18 | awk '{s+=$1} END {print s}'
```

spaceranger_anlysis_analysis.sh contain comeds to run spaceranger for four samples
```
./spaceranger_anlysis_analysis.sh
```



