#!/usr/local/bin/nextflow

params.file_path
params.file_name
params.publish_path
params.publish_name
params.trim_adapter
params.homo_sapiens_ref

filepairs_ch = Channel.fromFilePairs("${params.file_path}${params.file_name}", flat: true).view()

process PREPROCESS {

    container 'mwielsch/pre-process:v1.1'

    input:
    tuple val(reads_id), path(fastq_R1), path(fastq_R2)
    val trim_adapter
    val homo_sapiens_ref

    output:
    tuple val(reads_id),
          path("${reads_id}_host_removed_R1.fastq.gz"),
          path("${reads_id}_host_removed_R2.fastq.gz"), emit: processed_files
    tuple val(reads_id),
          path("${reads_id}_all_R1_readsin.txt"),
          path("${reads_id}_forward_paired.txt"), emit: count_qc

    script:
    """
    zgrep -E "^@(SAMN[0-9]+|[ES]RR[0-9]+)(\\.[0-9]+)?" $fastq_R1 | wc -l > ${reads_id}_all_R1_readsin.txt

    java -jar /usr/local/bin/trimmomatic.jar PE $fastq_R1 $fastq_R2 \
        ${reads_id}_forward_paired.fq.gz ${reads_id}_forward_unpaired.fq.gz \
        ${reads_id}_reverse_paired.fq.gz ${reads_id}_reverse_unpaired.fq.gz \
        ILLUMINACLIP:$trim_adapter:2:30:10 SLIDINGWINDOW:4:20 MINLEN:50 -threads 12

    zgrep -E "^@(SAMN[0-9]+|[ES]RR[0-9]+)(\\.[0-9]+)?" ${reads_id}_forward_paired.fq.gz | wc -l > ${reads_id}_forward_paired.txt
    zgrep -E "^@(SAMN[0-9]+|[ES]RR[0-9]+)(\\.[0-9]+)?" ${reads_id}_reverse_unpaired.fq.gz | wc -l > ${reads_id}_reverse_not_paired.txt
    zgrep -E "^@(SAMN[0-9]+|[ES]RR[0-9]+)(\\.[0-9]+)?" ${reads_id}_forward_unpaired.fq.gz | wc -l > ${reads_id}_forward_not_paired.txt

    rm ${reads_id}_forward_unpaired.fq.gz
    rm ${reads_id}_reverse_unpaired.fq.gz

    bowtie2 -p 12 -x $homo_sapiens_ref -1 ${reads_id}_forward_paired.fq.gz -2 ${reads_id}_reverse_paired.fq.gz \
    -S ${reads_id}mapped_plus_unmapped.sam --very-sensitive --dovetail

    rm ${reads_id}_forward_paired.fq.gz
    rm ${reads_id}_reverse_paired.fq.gz

    samtools view -bS ${reads_id}mapped_plus_unmapped.sam > ${reads_id}mapped_plus_unmapped.bam
    rm ${reads_id}mapped_plus_unmapped.sam

    samtools view -b -f 12 -F 256 ${reads_id}mapped_plus_unmapped.bam > ${reads_id}_bothReadsUnmapped.bam
    rm ${reads_id}mapped_plus_unmapped.bam

    samtools sort -n  ${reads_id}_bothReadsUnmapped.bam -o ${reads_id}_bothReadsUnmapped_sorted.bam
    rm ${reads_id}_bothReadsUnmapped.bam

    samtools fastq -@ 12 ${reads_id}_bothReadsUnmapped_sorted.bam \
    -1 ${reads_id}_host_removed_R1.fastq.gz \
    -2 ${reads_id}_host_removed_R2.fastq.gz \
    -0 /dev/null -s /dev/null -n

    rm ${reads_id}_bothReadsUnmapped_sorted.bam
    """
}

process SPADES {

    container 'mwielsch/spades-fast:v2.1'

    input:
    tuple val(reads_id),
          path("${reads_id}_host_removed_R1.fastq.gz"),
          path("${reads_id}_host_removed_R2.fastq.gz")

    output:
    tuple val(reads_id),
          path("${reads_id}_Spades_contigs.fasta"), emit: contigs_raw
    tuple val(reads_id),
          path("${reads_id}_final_contigs.txt"), emit: count_qc

    script:
    """
    spades.py --pe1-1 ${reads_id}_host_removed_R1.fastq.gz \
            --pe1-2 ${reads_id}_host_removed_R2.fastq.gz \
            -o \$PWD

    cat scaffolds.fasta | /opt/conda/envs/FAST/bin/faswc > ${reads_id}_all_contigs.txt
    cat scaffolds.fasta | /opt/conda/envs/FAST/bin/faslen | /opt/conda/envs/FAST/bin/fasfilter -t length 2000.. > ${reads_id}_Spades_contigs.fasta
    cat ${reads_id}_Spades_contigs.fasta | /opt/conda/envs/FAST/bin/faswc > ${reads_id}_final_contigs.txt
    """
}

process PREPPILON {

    container 'mwielsch/pre-process:v1.1'

    input:
    tuple val(reads_id),
          path("${reads_id}_Spades_contigs.fasta"),
          path("${reads_id}_host_removed_R1.fastq.gz"),
          path("${reads_id}_host_removed_R2.fastq.gz")

    output:
    tuple val(reads_id),
          path("${reads_id}_sorted.bam"), emit: sorted_bam

    script:
    """
    bowtie2-build ${reads_id}_Spades_contigs.fasta ${reads_id}_contigs_ref

    bowtie2 -x ${reads_id}_contigs_ref -1 ${reads_id}_host_removed_R1.fastq.gz\
    -2 ${reads_id}_host_removed_R2.fastq.gz |\
    samtools view -bS -o ${reads_id}.bam

    samtools sort ${reads_id}.bam -o ${reads_id}_sorted.bam
    samtools index ${reads_id}_sorted.bam
    """
}

process PILON {

    container 'staphb/pilon:1.24'

    publishDir "${params.publish_path}${params.publish_name}", pattern: "${reads_id}.fasta", mode: 'copy'

    input:
    tuple val(reads_id),
          path("${reads_id}_Spades_contigs.fasta"),
          path("${reads_id}_sorted.bam"),
          path("${reads_id}_final_contigs.txt"),
          path("${reads_id}_all_R1_readsin.txt"),
          path("${reads_id}_forward_paired.txt")

    output:
    tuple val(reads_id),
          path("${reads_id}.fasta"), emit: fasta
    path ("${reads_id}_qc_results.tsv"), emit: qc_results

    script:
    """
    java -jar /pilon/pilon.jar --genome ${reads_id}_Spades_contigs.fasta --frags ${reads_id}_sorted.bam --duplicates --output ${reads_id}

    r1_reads=\$(cat "${reads_id}_all_R1_readsin.txt")
    paired_reads=\$(cat "${reads_id}_forward_paired.txt")

    read spades_contigs spades_total_length _ < <(awk '{print \$1, \$2}' "${reads_id}_final_contigs.txt")

    echo -e "${reads_id}\t\${r1_reads}\t\${paired_reads}\t\${spades_contigs}\t\${spades_total_length}" >> "${reads_id}_qc_results.tsv"
    """
}

workflow {
    in_preprocess_ch = filepairs_ch
    preprocess_ch = PREPROCESS(filepairs_ch, params.trim_adapter, params.homo_sapiens_ref)

    in_spades_ch = preprocess_ch.processed_files
    spades_ch = SPADES(in_spades_ch)

    in_preppilon_ch = spades_ch.contigs_raw.join(preprocess_ch.processed_files)
    preppilon_ch = PREPPILON(in_preppilon_ch)

    in_pilon_ch = spades_ch.contigs_raw.join(preppilon_ch.sorted_bam).join(spades_ch.count_qc).join(preprocess_ch.count_qc)
    pilon_ch = PILON(in_pilon_ch)
}
