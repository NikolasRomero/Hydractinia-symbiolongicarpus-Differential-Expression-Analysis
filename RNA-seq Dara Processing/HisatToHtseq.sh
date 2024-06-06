# BAM file list (Treatment as example)
archivos_bam=("./BAM/Tratamiento/T1.bam" "./BAM/Tratamiento/T2.bam" "./BAM/Tratamiento/T3.bam" "./BAM/Tratamiento/T4.bam" "./BAM/Tratamiento/T5.bam" "./BAM/Tratamiento/T6.bam" "./BAM/Tratamiento/T7.bam" ".$

# Bucle for each BAM file
for bam in ${archivos_bam[@]}
do
    # Extract BAM file name to use it as output name
    base=$(basename $bam .bam)

    # RUN HT-seq-count for each BAM file & save the results as .txt
    htseq-count -f bam -r pos -s no -t exon -i transcript_id $bam ./Features/HSymV2.1.gtf > $base.txt

    # Move the output file to a new directory
    mv $base.txt transcript_counts/
done
