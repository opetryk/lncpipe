#!/usr/bin/env bash -C -e -u -o pipefail
printf "%s %s\n" sample1_R1.fastq.gz SAMPLE1_PE_1.gz sample1_R2.fastq.gz SAMPLE1_PE_2.gz | while read old_name new_name; do
    [ -f "${new_name}" ] || ln -s $old_name $new_name
done

fastqc \
    --quiet \
    --threads 4 \
    --memory 3840.0 \
    SAMPLE1_PE_1.gz SAMPLE1_PE_2.gz

cat <<-END_VERSIONS > versions.yml
"NFCORE_LNCPIPE:LNCPIPE:FASTQC":
    fastqc: $( fastqc --version | sed '/FastQC v/!d; s/.*v//' )
END_VERSIONS
