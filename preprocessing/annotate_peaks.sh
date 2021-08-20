#!/bin/bash
ACCESS_PATH=../data
ACCESS_PATH2=../data/macs3_results

echo "1. prepare file bed with ltr and gene peaks"
bedtools intersect -u -a $ACCESS_PATH2/merged_samples_peaks.narrowPeak -b $ACCESS_PATH/variation_and_repeats.bed  2>/dev/null | awk '{print $0"\tltr"}' > $ACCESS_PATH/tmp_peaks_ltr.bed
bedtools intersect -v -a $ACCESS_PATH2/merged_samples_peaks.narrowPeak -b $ACCESS_PATH/variation_and_repeats.bed  2>/dev/null | awk '{print $0"\tgene"}' > $ACCESS_PATH/tmp_peaks_gene.bed

echo "2. prepare file bamwith minus and plus strand"
samtools view -b -f 16 $ACCESS_PATH/merged_four_samples_sorted.bam > $ACCESS_PATH/tmp_merged_samples_minus.bam
samtools view -b -F 16 $ACCESS_PATH/merged_four_samples_sorted.bam > $ACCESS_PATH/tmp_merged_samples_plus.bam

echo "3. indexing bam files"
samtools index $ACCESS_PATH/tmp_merged_samples_minus.bam
samtools index $ACCESS_PATH/tmp_merged_samples_plus.bam

echo "4. creating file for ltr with minus and plus coverage"
samtools bedcov $ACCESS_PATH/tmp_peaks_ltr.bed $ACCESS_PATH/tmp_merged_samples_minus.bam > $ACCESS_PATH/tmp_peaks_ltr_coverage.minus.bed
samtools bedcov $ACCESS_PATH/tmp_peaks_ltr.bed $ACCESS_PATH/tmp_merged_samples_plus.bam > $ACCESS_PATH/tmp_peaks_ltr_coverage.plus.bed

echo "5. creating file for gene with minus and plus coverage"
samtools bedcov $ACCESS_PATH/tmp_peaks_gene.bed $ACCESS_PATH/tmp_merged_samples_minus.bam > $ACCESS_PATH/tmp_peaks_gene_coverage.minus.bed
samtools bedcov $ACCESS_PATH/tmp_peaks_gene.bed $ACCESS_PATH/tmp_merged_samples_plus.bam > $ACCESS_PATH/tmp_peaks_gene_coverage.plus.bed

echo "6. assessment of strand for ltr"
paste $ACCESS_PATH/tmp_peaks_ltr_coverage.plus.bed $ACCESS_PATH/tmp_peaks_ltr_coverage.minus.bed |
    awk '{print $0"\t+"$12-$24}' |
    awk 'BEGIN{FS=OFS"\t"} {gsub(/+-[0-9]*/, "-" $3)} 1 {gsub(/+[0-9]*/, "+" $3)} 1' |
    awk -F"\t" '{OFS=FS}{ $6=$25 ; print   }' |
    cut -f1-11 > $ACCESS_PATH/tmp_peaks_ltr_strand.bed

echo "7. assessment of strand for gene"
paste $ACCESS_PATH/tmp_peaks_gene_coverage.plus.bed $ACCESS_PATH/tmp_peaks_gene_coverage.minus.bed |
    awk '{print $0"\t+"$12-$24}' |
    awk 'BEGIN{FS=OFS"\t"} {gsub(/+-[0-9]*/, "-" $3)} 1 {gsub(/+[0-9]*/, "+" $3)} 1' |
    awk -F"\t" '{OFS=FS}{ $6=$25 ; print   }' |
    cut -f1-11  > $ACCESS_PATH/tmp_peaks_gene_strand.bed

echo "8. prepare file bed to convert the format gtf for ltr and gene"
bedtools sort -i $ACCESS_PATH/mart_export_v102_mm10.bed  > $ACCESS_PATH/tmp_mart_export_v102_mm10_sorted.bed

bedtools closest -a $ACCESS_PATH/tmp_peaks_gene_strand.bed -b $ACCESS_PATH/tmp_mart_export_v102_mm10_sorted.bed -t first 2>/dev/null > $ACCESS_PATH/tmp_peaks2gtf_gene.bed
bedtools closest -a $ACCESS_PATH/tmp_peaks_ltr_strand.bed -b $ACCESS_PATH/tmp_mart_export_v102_mm10_sorted.bed -t first 2>/dev/null > $ACCESS_PATH/tmp_peaks2gtf_ltr.bed

echo "9. connect file bed for ltr and gene"
cat $ACCESS_PATH/tmp_peaks2gtf_gene.bed $ACCESS_PATH/tmp_peaks2gtf_ltr.bed | grep -Pv "GL4|JH5" > $ACCESS_PATH/tmp_peaks_annotate.bed


echo "10. sort file bed"
bedtools sort -i $ACCESS_PATH/tmp_peaks_annotate.bed > $ACCESS_PATH/peaks_annotate_sorted.bed.test

rm $ACCESS_PATH/tmp*
