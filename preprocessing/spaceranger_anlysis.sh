#!/bin/bash
#REF_PATH=$1
PREFIX=$2
REF_PATH=$1
OUTPUT_PATH=data/spaceranger_results
#REF_PATH=/home/ifpan/projects/ifpan-annaradli-spatialtrial/corrected_reference/
FASTQ_PATH=raw/data/00_untrimmed_fastq/
IMAGE_PATH=raw/data
NUMBER_CORES=6

spaceranger count --id=$PREFIX\S3647Nr1 \
                  --transcriptome=$REF_PATH \
                  --fastqs=$FASTQ_PATH \
                  --image=$IMAGE_PATH/2021-03-25\ VISIUM\ EXPRESSION-A1_A1-2\ Merged.tif \
                  --slide=V10T31-004 \
                  --area=A1 \
                  --loupe-alignment=$IMAGE_PATH/V10T31-004-A1.json \
                  --localcores=$NUMBER_CORES \
                  --sample=S3647Nr1

mv /$PREFIX/S3647Nr1 $OUTPUT_PATH
exit 0
spaceranger count --id=corrected_S3647Nr2 \
                  --transcriptome=$REF_PATH \
                  --fastqs=$FASTQ_PATH \
                  --image=$IMAGE_PATH/2021-03-25\ VISIUM\ EXPRESSION\ B1_B1\ Merged_RAW_ch00.tif  \
                  --slide=V10T31-004 \
                  --area=B1 \
                  --loupe-alignment=$IMAGE_PATH/V10T31-004-B1.json \
                  --localcores=$NUMBER_CORES \
                  --sample=S3647Nr2

mv /$PREFIX/S3647Nr2

spaceranger count --id=corrected_S3647Nr3 \
                  --transcriptome=$REF_PATH \
                  --fastqs=$FASTQ_PATH \
                  --image=$IMAGE_PATH/2021-03-25\ VISIUM\ EXPRESSION\ C1_C1\ Merged_RAW_ch00.tif \
                  --slide=V10T31-004 --area=C1 \
                  --loupe-alignment=$IMAGE_PATH/V10T31-004-C1.json \
                  --localcores=$NUMBER_CORES \
                  --sample=S3647Nr3

mv /$PREFIX/S3647Nr3

spaceranger count --id=corrected_S3647Nr4 \
                  --transcriptome=$REF_PATH \
                  --fastqs=$FASTQ_PATH \
                  --image=$IMAGE_PATH/2021-03-25\ VISIUM\ EXPRESSION\ D1_D1\ Merged_RAW_ch00.tif \
                  --slide=V10T31-004 \
                  --area=D1 \
                  --loupe-alignment=$IMAGE_PATH/V10T31-004-D1.json \
                  --localcores=$NUMBER_CORES \
                  --sample=S3647Nr4

mv /$PREFIX/S3647Nr4



