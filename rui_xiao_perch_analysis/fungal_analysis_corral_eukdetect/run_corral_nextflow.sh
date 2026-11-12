#!/bin/bash
# run_corral_nextflow.sh
# Save as: run_corral_nextflow.sh && chmod +x run_corral_nextflow.sh && ./run_corral_nextflow.sh


nextflow run wbazant/CORRAL -r main \
    --inputPath in.tsv \
    --resultDir /data/ruicatxiao/perch/corral_results \
    --downloadMethod local --libraryLayout paired \
    --refdb /data/ruicatxiao/perch/eukdb/eukdb \
    --markerToTaxonPath /data/ruicatxiao/perch/eukdb/busco_taxid_link.txt \
    -c ~/nfcore.config
