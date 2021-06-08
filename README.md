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
