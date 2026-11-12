while read -r sampleid; do
    echo "${sampleid} start"
    
    kraken2 \
        --threads 32 \
        --db k2db \
        --gzip-compressed \
        --output "k2report/${sampleid}_k2.kraken" \
        --report "k2report/${sampleid}_k2.txt" \
        "${sampleid}-contigs.fa.gz"
    
    echo "${sampleid} done"
done < sample_id.txt
