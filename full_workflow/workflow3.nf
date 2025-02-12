#!/usr/local/bin/nextflow

params.file_path
params.file_name
params.trim_adapter
params.homo_sapiens_ref
params.quast_ref

params.checkm_db
params.saureus_ref
params.genomad_db
params.checkv_db
params.pharokka_db

filepairs_ch = Channel.fromFilePairs("${params.file_path}${params.file_name}", flat: true).view()

include { PREPROCESS; SPADES; PREPPILON; PILON; QUASTCOLLECT; QCRESULTS } from './subworkflow2.nf'

process CHECKM {

    container 'staphb/checkm:1.2.3'

    input:
    path fasta_files
    val checkm_db

    output:
    path "bin_stats_ext.tsv", emit: bin_stats

    script:
    """
    mkdir fasta_dir
    cp *.fasta fasta_dir

    export CHECKM_DATA_PATH=${checkm_db}
    checkm lineage_wf -t $task.cpus -x fasta --reduced_tree fasta_dir \$PWD

    mv \$PWD/storage/bin_stats_ext.tsv \$PWD
    """
}

process RAGOUTASSEMBLY {

    container 'mwielsch/ragout:v1.0'

    input:
    tuple val(reads_id),
          path("${reads_id}.fasta"),
          path("${reads_id}_final_contigs.txt"),
          path("${reads_id}_all_R1_readsin.txt"),
          path("${reads_id}_forward_paired.txt")
    val saureusref
    path bin_stats_path

    output:
    path ("NCBI_SA_2025_${reads_id}.tsv"), emit: qc_results
    tuple val(reads_id),
          path("${reads_id}_scaffolds.fasta"), emit: scaffolds

    script:
    """
    echo ".references = S_aureus" > ragout_${reads_id}.rcp
    echo ".target = ${reads_id}" >> ragout_${reads_id}.rcp
    echo "S_aureus.fasta = ${saureusref}" >> ragout_${reads_id}.rcp
    echo "${reads_id}.fasta = ${reads_id}.fasta" >> ragout_${reads_id}.rcp
    ragout ragout_${reads_id}.rcp --outdir \$PWD --refine

    r1_reads=\$(cat "${reads_id}_all_R1_readsin.txt")
    paired_reads=\$(cat "${reads_id}_forward_paired.txt")
    read spades_contigs spades_total_length _ < <(awk '{print \$1, \$2}' "${reads_id}_final_contigs.txt")

    stats=\$(grep "^${reads_id}" "${bin_stats_path}")
    marker_lineage=\$(echo \${stats} | sed -n "s/.*'marker lineage': '\\([^']*\\)'.*/\\1/p")
    completeness=\$(echo \${stats} | sed -n "s/.*'Completeness': \\([^,]*\\).*/\\1/p")
    contamination=\$(echo \${stats} | sed -n "s/.*'Contamination': \\([^,]*\\).*/\\1/p")
    gc=\$(echo \${stats} | sed -n "s/.*'GC': \\([^,]*\\).*/\\1/p")
    genome_size=\$(echo \${stats} | sed -n "s/.*'Genome size': \\([^,]*\\).*/\\1/p")
    contigs=\$(echo \${stats} | sed -n "s/.*'# contigs': \\([^,]*\\).*/\\1/p")
    longest_contig=\$(echo \${stats} | sed -n "s/.*'Longest contig': \\([^,]*\\).*/\\1/p")
    n50_contigs=\$(echo \${stats} | sed -n "s/.*'N50 (contigs)': \\([^,]*\\).*/\\1/p")
    mean_contig_length=\$(echo \${stats} | sed -n "s/.*'Mean contig length': \\([^,]*\\).*/\\1/p")
    coding_density=\$(echo \${stats} | sed -n "s/.*'Coding density': \\([^,]*\\).*/\\1/p")
    predicted_genes=\$(echo \${stats} | sed -n "s/.*'# predicted genes': \\([^,]*\\).*/\\1/p")

    scaffolds=\$(grep "Scaffolds:" "ragout.log" | awk '{print \$2}')
    draft_genome_length=\$(grep "Scaffolds length:" "ragout.log" | awk '{print \$3}')
    unplaced_fragments=\$(grep "Unplaced fragments:" "ragout.log" | awk '{print \$3}')
    unplaced_length=\$(grep "Unplaced length:" "ragout.log" | awk '{print \$3}')
    unplaced_length_in_PERC=\$(grep "Unplaced length:" "ragout.log" | awk '{print \$4}' | tr -d '()%')
    introduced_Ns_length=\$(grep "Introduced Ns length:" "ragout.log" | awk '{print \$4}')
    introduced_Ns_length_in_PERC=\$(grep "Introduced Ns length:" "ragout.log" | awk '{print \$5}' | tr -d '()%')
    n50_all=\$(grep "Fragments N50:" "ragout.log" | awk '{print \$3}')
    n50_draftGenome=\$(grep "Assembly N50:" "ragout.log" | awk '{print \$3}')

    echo -e "${reads_id}\t\${r1_reads}\t\${paired_reads}\t\${spades_contigs}\t\${spades_total_length}\t\${marker_lineage}\t\${completeness}\t\${contamination}\t\${gc}\t\${genome_size}\t\${contigs}\t\${longest_contig}\t\${n50_contigs}\t\${mean_contig_length}\t\${coding_density}\t\${predicted_genes}\t\${scaffolds}\t\${draft_genome_length}\t\${unplaced_fragments}\t\${unplaced_length}\t\${unplaced_length_in_PERC}\t\${introduced_Ns_length}\t\${introduced_Ns_length_in_PERC}\t\${n50_all}\t\${n50_draftGenome}" > "NCBI_SA_2025_${reads_id}.tsv"
    """
}

process MERGE {

    def timestamp = new java.util.Date().format( 'yyyy-MM-dd_HH-mm')

    publishDir "publish_workflow", pattern: "NCBI_SA_2025_${timestamp}.tsv", mode: 'copy'

    input:
    path result_files

    output:
    path "NCBI_SA_2025_${timestamp}.tsv", emit: final_result

    script:
    """
    echo -e "ID\tR1_reads\tpaired_reads\tSpades_contigs\tspades_total_length\tmarker_lineage\tCompleteness\tContamination\tGC\tcheckM_genomeSize\tcheckM_contigs\tLongest_contig\tN50_contigs\tMean_contig_length\tCoding_density\tPredicted_genes\tS_A_draftGenomes\tdraft_genome_length\tUnplaced_fragments\tUnplaced_length\tUnplaced_length_in_PERC\tIntroduced_Ns_length\tIntroduced_Ns_length_in_PERC\tN50_all\tN50_draftGenome" > NCBI_SA_2025_${timestamp}.tsv

    cat $result_files >> NCBI_SA_2025_${timestamp}.tsv
    """
}

process MLST {

    container 'staphb/mlst:2.23.0-2025-01-01'

    publishDir "publish_workflow", pattern: "${reads_id}_MLST.txt", mode: 'copy'

    input:
    tuple val(reads_id),
          path("${reads_id}_scaffolds.fasta")

    output:
    tuple val(reads_id),
          path("${reads_id}_MLST.txt"), emit: typing_result

    script:
    """
    mlst ${reads_id}_scaffolds.fasta --scheme saureus > ${reads_id}_MLST.txt

    """
}

process AMRFINDER {

    container 'staphb/ncbi-amrfinderplus:4.0.15-2024-12-18.1'

    input:
    tuple val(reads_id),
          path("${reads_id}_scaffolds.fasta")

    output:
    tuple val(reads_id),
          path("${reads_id}_amr.txt"), emit: amr_txt

    script:
    """
    amrfinder --plus -n ${reads_id}_scaffolds.fasta -O Staphylococcus_aureus --print_node > ${reads_id}_AMRFinderPlus.txt

    awk -v id="${reads_id}" 'BEGIN {
    FS="\\t"; OFS="|";
    print "ID", "Start", "Stop", "Strand", "Gene", "Element_type", "Element_subtype", "Target_length", "Reference_sequence_length", "perc_Identity_to_reference_sequence", "Scope", "Product"
    }
    NR > 1 {
    start = \$3
    stop = \$4
    strand = \$5
    gene = \$6
    element_type = \$9
    element_subtype = \$10
    target_length = \$14
    reference_length = \$15
    perc_identity = \$17
    scope = \$7
    product = \$8

    print id, start, stop, strand, gene, element_type, element_subtype, target_length, reference_length, perc_identity, product, scope
    }' "${reads_id}_AMRFinderPlus.txt" > "${reads_id}_amr.txt"

    sed -i 's/ /_/g' "${reads_id}_amr.txt"
    sed -i 's/\\//_/g' "${reads_id}_amr.txt"
    """
}

process PROKKA {

    container 'mwielsch/prokka:v1.2'

    input:
    tuple val(reads_id),
          path("${reads_id}_scaffolds.fasta"),
          path("${reads_id}_MLST.txt")

    output:
    tuple val(reads_id),
          path("${reads_id}_prokka.txt"), emit: prokka_txt

    script:
    """
    prokka --outdir \$PWD/${reads_id} --prefix ${reads_id} ${reads_id}_scaffolds.fasta

    mv \$PWD/${reads_id}/${reads_id}.gff \$PWD

    declare mlst_data
    while IFS=\$'\\t' read -r id organism mlst rest; do
      if [[ "\$mlst" == "-" ]]; then
        mlst_data="unknown"
      else
        mlst_data=\$mlst
      fi
    done < ${reads_id}_MLST.txt

    mlst=\${mlst_data}
    if [[ -z "\$mlst" ]]; then
      mlst="NA"
    fi

    temp_file=\$(mktemp)

    awk '/^##FASTA/ {exit} {print}' "${reads_id}.gff" > "\$temp_file"

    awk -v id="${reads_id}" -v mlst="\$mlst" 'BEGIN {
      FS="\\t"; OFS="|";
      print "ID", "MLST", "Start", "End", "Strand", "Gene", "UniProtKB", "COG", "Product"
    }
    /^#/ { next }
    {
      start = \$4
      end = \$5
      strand = \$7
      gene = "NA"
      uniprot = "NA"
      cog = "NA"
      product = "NA"

      split(\$9, attributes, ";")
      for (i in attributes) {
        split(attributes[i], kv, "=")
        if (kv[1] == "gene") gene = kv[2]
        if (kv[1] == "inference" && kv[2] ~ /UniProtKB:/) {
        match(kv[2], /UniProtKB:([^,;]+)/, arr)
        if (arr[1] != "") {
          uniprot = arr[1]
        }
     }
     if (kv[1] == "db_xref" && kv[2] ~ /COG:/) {
      split(kv[2], cog_info, ":")
      cog = cog_info[2]
     }
     if (kv[1] == "product") product = kv[2]
     }

     print id, "ST" mlst, start, end, strand, gene, uniprot, cog, product
    }' "\$temp_file" > "${reads_id}_prokka.txt"

  rm "\$temp_file"
  sed -i 's/ /_/g' "${reads_id}_prokka.txt"
  sed -i 's/\\//_/g' "${reads_id}_prokka.txt"
    """
}

process GENOMAD {

    container 'staphb/genomad:1.8.1'

    input:
    tuple val(reads_id),
          path("${reads_id}_scaffolds.fasta")
    val genomad_db

    output:
    tuple val(reads_id),
          path("${reads_id}_provirus_*.fna"), emit: provirus

    script:
    """
    genomad end-to-end --cleanup --splits 10 ${reads_id}_scaffolds.fasta \$PWD ${genomad_db}

    mv \$PWD/${reads_id}_scaffolds_summary/${reads_id}_scaffolds_virus.fna \$PWD

    awk -v mw_number="${reads_id}" -v output_dir="\$PWD" '
    BEGIN { RS = ">" ; FS = "\\n" }
    NR > 1 {
      split(\$1, header_parts, "|")
      id = mw_number "_" header_parts[2]
      output_file = output_dir "/" id ".fna"
      print ">" id > output_file
      for (i = 2; i <= NF; i++) {
        printf "%s", \$i > output_file
      }
      printf "\\n" > output_file
    }' "${reads_id}_scaffolds_virus.fna"
    """
}

process PADLOC {

    container 'mwielsch/padloc:v1.4'

    input:
    tuple val(reads_id),
          path("${reads_id}_scaffolds.fasta")

    output:
    tuple val(reads_id),
          path("${reads_id}_padloc.txt"), emit: padloc_txt

    script:
    """
    padloc --fna ${reads_id}_scaffolds.fasta --outdir \$PWD --cpu $task.cpus

    awk -v id="${reads_id}" 'BEGIN {
    FS="\\t"; OFS="|";
    print "ID", "Start", "End", "Strand", "Gene", "Target_coverage", "P_value", "PDLC_id"
    }
    /^#/ { next }
     {
       start = \$4
       end = \$5
       strand = \$7
       gene = "NA"
       target_coverage = "NA"
       p_value = \$6
       p_value = sprintf("%.3e", p_value)
       pdlc_id = "NA"

       split(\$9, attributes, ";")
       for (i in attributes) {
       split(attributes[i], kv, "=")
       if (kv[1] == "Name") gene = kv[2]
       if (kv[1] == "Target.coverage") target_coverage = kv[2]
       if (kv[1] == "HMM.accession") pdlc_id = kv[2]
      }

    print id, start, end, strand, gene, target_coverage, p_value, pdlc_id
    }' "${reads_id}_scaffolds.fasta_padloc.gff" > "${reads_id}_padloc.txt"
    sed -i 's/ /_/g' "${reads_id}_padloc.txt"
    sed -i 's/\\//_/g' "${reads_id}_padloc.txt"
    """
}

process PHAROKKA {

    container 'jinlongru/pharokka:v1.7.3.3'

    input:
    tuple val(reads_id),
          path(provirus)
    val pharokka_db

    output:
    tuple val(reads_id),
          path("${reads_id}*.pharokka_out"), emit: pharokka_out

    script:
    """
    for file in ${reads_id}_provirus_*.fna; do

    id=\$(basename "\${file}" .fna)

    pharokka.py -i \${file} -o \${id}.pharokka_out -d ${pharokka_db} -t 4 --force -p \${id}

    done
    """
}

process PHAROTTER {

    container 'jinlongru/pharokka:v1.7.3.3'

    publishDir "publish_workflow", pattern: "${reads_id}_merged.txt", mode: 'copy'

    input:
    tuple val(reads_id),
          path(provirus),
          path(pharokka_out),
          path("${reads_id}_amr.txt"),
          path("${reads_id}_prokka.txt"),
          path("${reads_id}_padloc.txt")

    output:
    tuple val(reads_id),
          path("${reads_id}_merged.txt"), emit: merged_txt

    script:
    """
    for file in ${reads_id}_provirus_*.fna; do

    id=\$(basename "\${file}" .fna)

    pharokka_plotter.py -i \${file} -n \${id} -p \${id} -o \${id}.pharokka_out -f

    done

    mkdir -p temp

    header_printed=""

    intermediate_file="temp/${reads_id}_combined.gff"

    > "\$intermediate_file"

    for file in ${reads_id}_provirus_*.pharokka_out; do
      sub_id=\$(basename "\$file" .pharokka_out)
      gff_file="\${sub_id}.pharokka_out/\${sub_id}.gff"
      if [[ -f "\$gff_file" ]]; then
        awk '/^##FASTA/ {exit} {print}' "\$gff_file" >> "\$intermediate_file"
      fi
    done

    if [[ -z "\${header_printed}" ]]; then
      echo -e "ID|Start|End|Strand|Phage_start|Phage_stop|Gene_function|Product|vfdb_description|vfdb_species" > "${reads_id}_prophages.txt"
      header_printed=1
    fi

    awk -v id="${reads_id}" 'BEGIN {
      FS="\\t"; OFS="|";
    }
    /^##/ { next }
    {
      split(\$1, arr, "_");
      base_start = arr[3];
      phage_start = \$4;
      phage_stop = \$5;
      start = base_start + phage_start;
      end = base_start + phage_stop;
      strand = \$7;
      gene_function = "NA";
      product = "NA";
      vfdb_description = "NA";
      vfdb_species = "NA";

      split(\$9, attributes, ";");
      for (i in attributes) {
        split(attributes[i], kv, "=");
        gsub(/"/, "", kv[2]);
        if (kv[1] == "function") gene_function = kv[2];
        if (kv[1] == "product") product = kv[2];
        if (kv[1] == "vfdb_description") vfdb_description = kv[2];
        if (kv[1] == "vfdb_species") vfdb_species = kv[2];
      }

    print id, start, end, strand, phage_start, phage_stop, gene_function, product, vfdb_description, vfdb_species;
    }' "\$intermediate_file" >> "${reads_id}_prophages.txt"

    sed -i 's/ /_/g' "${reads_id}_prophages.txt"
    sed -i 's/\\//_/g' "${reads_id}_prophages.txt"

    rm "\$intermediate_file"

    prokka_file="${reads_id}_prokka.txt"
    amr_file="${reads_id}_amr.txt"
    padloc_file="${reads_id}_padloc.txt"
    prophage_file="${reads_id}_prophages.txt"
    output_file="${reads_id}_merged.txt"

    declare -A amr_data
    amr_count=0
    while IFS='|' read -r -a amr_fields; do
      if [[ "\${amr_fields[0]}" == "ID" ]]; then
        continue
      fi
      amr_key="\${amr_fields[1]}_\${amr_fields[3]}"
      amr_data["\$amr_key"]="\${amr_fields[*]}"
      amr_count=\$((amr_count + 1))
    done < "\$amr_file"

    echo "Loaded \$amr_count AMR entries for ID: ${reads_id}"

    declare -A padloc_data
    padloc_count=0
    while IFS='|' read -r -a padloc_fields; do
      if [[ "\${padloc_fields[0]}" == "ID" ]]; then
        continue
      fi
      padloc_key="\${padloc_fields[1]}_\${padloc_fields[3]}"
      padloc_data["\$padloc_key"]="\${padloc_fields[*]}"
      padloc_count=\$((padloc_count + 1))
    done < "\$padloc_file"

    echo "Loaded \$padloc_count Padloc entries for ID: ${reads_id}"

    declare -A prophage_data
    prophage_count=0
    while IFS='|' read -r -a prophage_fields; do
      if [[ "\${prophage_fields[0]}" == "ID" ]]; then
        continue
      fi
      prophage_key="\${prophage_fields[1]}_\${prophage_fields[3]}"
      prophage_data["\$prophage_key"]="\${prophage_fields[*]}"
      prophage_count=\$((prophage_count + 1))
    done < "\$prophage_file"

    echo "Loaded \$prophage_count Prophage entries for ID: ${reads_id}"

    header=\$(head -n 1 "\$prokka_file")
    echo -e "\${header}|amr_Gene|amr_Element_type|amr_Element_subtype|amr_Target_length|amr_Reference_sequence_length|amr_perc_Identity_to_reference_sequence|amr_Product|amr_Scope|padloc_Gene|padloc_Target_coverage|padloc_P_value|padloc_PDLC_id|prophage_Phage_start|prophage_Phage_stop|prophage_Gene_function|prophage_Product|prophage_vfdb_description|prophage_vfdb_species" > "\$output_file"

    tail -n +2 "\$prokka_file" | while IFS=\$'|' read -r -a prokka_fields; do
      prokka_key="\${prokka_fields[2]}_\${prokka_fields[4]}"
      amr_data_to_add="NA|NA|NA|NA|NA|NA|NA|NA"
      padloc_data_to_add="NA|NA|NA|NA"
      prophage_data_to_add="NA|NA|NA|NA|NA|NA"

      for key in "\${!amr_data[@]}"; do
        split_key=(\${key//_/ })
        amr_start=\${split_key[0]}
        amr_strand=\${split_key[1]}

        if [[ "\${prokka_fields[4]}" == "\$amr_strand" && (\${prokka_fields[2]} -ge \$((amr_start - 15)) && \${prokka_fields[2]} -le \$((amr_start + 15))) ]]; then
          amr_entry="\${amr_data[\$key]}"
          amr_fields=(\${amr_entry//|/ })
          amr_data_to_add="\${amr_fields[4]}|\${amr_fields[5]}|\${amr_fields[6]}|\${amr_fields[7]}|\${amr_fields[8]}|\${amr_fields[9]}|\${amr_fields[10]}|\${amr_fields[11]}"
          echo "Match found for AMR entry: \${amr_data[\$key]}"  # Print matched AMR entries
          break
        fi
      done

      for key in "\${!padloc_data[@]}"; do
        split_key=(\${key//_/ })
        padloc_start=\${split_key[0]}
        padloc_strand=\${split_key[1]}

        if [[ "\${prokka_fields[4]}" == "\$padloc_strand" && (\${prokka_fields[2]} -ge \$((padloc_start - 15)) && \${prokka_fields[2]} -le \$((padloc_start + 15))) ]]; then
          padloc_entry="\${padloc_data[\$key]}"
          padloc_fields=(\${padloc_entry//|/ })
          padloc_data_to_add="\${padloc_fields[4]}|\${padloc_fields[5]}|\${padloc_fields[6]}|\${padloc_fields[7]}"
          echo "Match found for Padloc entry: \${padloc_data[\$key]}"
          break
        fi
      done

      for key in "\${!prophage_data[@]}"; do
        split_key=(\${key//_/ })
        prophage_start=\${split_key[0]}
        prophage_strand=\${split_key[1]}

        if [[ "\${prokka_fields[4]}" == "\$prophage_strand" && (\${prokka_fields[2]} -ge \$((prophage_start - 15)) && \${prokka_fields[2]} -le \$((prophage_start + 15))) ]]; then
          prophage_entry="\${prophage_data[\$key]}"
          prophage_fields=(\${prophage_entry//|/ })
          prophage_data_to_add="\${prophage_fields[4]}|\${prophage_fields[5]}|\${prophage_fields[6]}|\${prophage_fields[7]}|\${prophage_fields[8]}|\${prophage_fields[9]}"
          echo "Match found for Prophage entry: \${prophage_data[\$key]}"
          break
        fi
      done

      echo -e "\${prokka_fields[*]}|\$amr_data_to_add|\$padloc_data_to_add|\$prophage_data_to_add" >> "\$output_file"
    done

    sed -i 's/ /|/g' "\$output_file"
    sed -i 's/|/\t/g' "\$output_file"
    """
  }

workflow {
    in_preprocess_ch = filepairs_ch
    preprocess_ch = PREPROCESS(filepairs_ch, params.trim_adapter, params.homo_sapiens_ref)

    in_spades_ch = preprocess_ch.processed_files
    spades_ch = SPADES(in_spades_ch)

    in_preppilon_ch = spades_ch.contigs_raw.join(preprocess_ch.processed_files)
    preppilon_ch = PREPPILON(in_preppilon_ch)

    in_pilon_ch = spades_ch.contigs_raw.join(preppilon_ch.sorted_bam)
    pilon_ch = PILON(in_pilon_ch)

    in_quastcollect_ch = pilon_ch.fasta.join(spades_ch.count_qc).join(preprocess_ch.count_qc)
    quastcollect_ch = QUASTCOLLECT(in_quastcollect_ch, params.quast_ref)

    in_qcresults_ch = quastcollect_ch.qc_results.collect()
    qcresults_ch = QCRESULTS(in_qcresults_ch)

    in_checkm_ch = pilon_ch.fasta.map { it[1] }.collect()
    checkm_ch = CHECKM(in_checkm_ch, params.checkm_db)

    in_ragoutcollect_ch = pilon_ch.fasta.join(spades_ch.count_qc).join(preprocess_ch.count_qc)
    ragoutcollect_ch = RAGOUTASSEMBLY(in_ragoutcollect_ch, params.saureus_ref, checkm_ch.bin_stats.collect())

    in_merge_ch = ragoutcollect_ch.qc_results.collect()
    merge_ch = MERGE(in_merge_ch)

    in_mlst_ch = ragoutcollect_ch.scaffolds
    mlst_ch = MLST(in_mlst_ch)

    in_amrfinder_ch = ragoutcollect_ch.scaffolds
    amrfinder_ch = AMRFINDER(in_amrfinder_ch)

    in_prokka_ch = ragoutcollect_ch.scaffolds.join(mlst_ch.typing_result)
    prokka_ch = PROKKA(in_prokka_ch)

    in_genomad_ch = ragoutcollect_ch.scaffolds
    genomad_ch = GENOMAD(in_genomad_ch,params.genomad_db)

    in_padloc_ch = ragoutcollect_ch.scaffolds
    padloc_ch = PADLOC(in_padloc_ch)

    in_pharokka_ch = genomad_ch.provirus
    pharokka_ch = PHAROKKA(in_pharokka_ch, params.pharokka_db)

    in_pharotter_ch = genomad_ch.provirus.join(pharokka_ch.pharokka_out).join(amrfinder_ch.amr_txt).join(prokka_ch.prokka_txt).join(padloc_ch.padloc_txt)
    pharotter_ch = PHAROTTER(in_pharotter_ch)
}
