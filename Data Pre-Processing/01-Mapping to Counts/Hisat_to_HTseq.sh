## BAM file list reeplace with your actual BAM files list
archivos_bam=(LIST_OF_BAM)
## Bucle function for each BAM file in list
for bam in ${archivos_bam[@]}
do
    ## Extract the BAM file name to use it as the output name
    base=$(basename $bam .bam)

    ## Run HT-seq-count for each BAM file and save the results as .txt files
    htseq-count -f bam -r pos -s no -t exon -i transcript_id $bam ./Features/HSymV2.1.gtf > $base.txt

    ## Move the output file to a new directory
    mv $base.txt transcript_counts/
done
