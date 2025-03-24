/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_lncpipe_pipeline'
include { checkSamplesAfterGrouping  } from '../subworkflows/local/utils_nfcore_lncpipe_pipeline'
include { samplesheetToList                } from 'plugin/nf-schema'
include { FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS  } from '../subworkflows/nf-core/fastq_qc_trim_filter_setstrandedness'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow LNCPIPE {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    ch_fasta             // channel: path(genome.fasta)
    ch_gtf               // channel: path(genome.gtf)
    ch_fai               // channel: path(genome.fai)
    ch_chrom_sizes       // channel: path(genome.sizes)
    ch_gene_bed          // channel: path(gene.bed)
    ch_transcript_fasta  // channel: path(transcript.fasta)
    ch_star_index        // channel: path(star/index/)
    ch_rsem_index        // channel: path(rsem/index/)
    ch_hisat2_index      // channel: path(hisat2/index/)
    ch_salmon_index      // channel: path(salmon/index/)
    ch_kallisto_index    // channel: [ meta, path(kallisto/index/) ]
    ch_bbsplit_index     // channel: path(bbsplit/index/)
    ch_ribo_db           // channel: path(sortmerna_fasta_list)
    ch_sortmerna_index   // channel: path(sortmerna/index/)
    ch_splicesites       // channel: path(genome.splicesites.txt)
    make_sortmerna_index // boolean: Whether to create an index before running sortmerna
    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // Create channel from input file provided through params.input
    //
    Channel
        .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
        .map {
            meta, fastq_1, fastq_2 ->
                if (!fastq_2) {
                    return [ meta.id, meta + [ single_end:true ], [ fastq_1 ] ]
                } else {
                    return [ meta.id, meta + [ single_end:false ], [ fastq_1, fastq_2 ] ]
                }
        }
        .groupTuple()
        .map { samplesheet ->
            checkSamplesAfterGrouping(samplesheet)
        }
        .set { ch_fastq }

// Checking parameters
// ...





/*
* Step 3: QC (FastQC/AfterQC/Fastp) of raw reads
*/

    //
    // MODULE: FASTP
    //
    // ch_adapters = params.adapters ? params.adapters : []

    // FASTP (
    //     ch_samplesheet,
    //     ch_adapters,
    //     params.discard_trimmed_pass,
    //     params.save_trimmed_fail,
    //     params.save_merged
    // )
    // ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.collect{it[1]}.ifEmpty([]))
    // ch_versions      = ch_versions.mix(FASTP.out.versions.first())


    //
    // Run RNA-seq FASTQ preprocessing subworkflow
    //

    // The subworkflow only has to do Salmon indexing if it discovers 'auto'
    // samples, and if we haven't already made one elsewhere
    salmon_index_available = params.salmon_index || (!params.skip_pseudo_alignment && params.pseudo_aligner == 'salmon')

    FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS (
        ch_fastq,
        ch_fasta,
        ch_transcript_fasta,
        ch_gtf,
        ch_salmon_index,
        ch_sortmerna_index,
        ch_bbsplit_index,
        ch_ribo_db,
        params.skip_bbsplit,
        params.skip_fastqc || params.skip_qc,
        params.skip_trimming,
        params.skip_umi_extract,
        !salmon_index_available,
        !params.sortmerna_index && params.remove_ribo_rna,
        params.trimmer,
        params.min_trimmed_reads,
        params.save_trimmed,
        params.remove_ribo_rna,
        params.with_umi,
        params.umi_discard_read,
        params.stranded_threshold,
        params.unstranded_threshold,
        params.skip_linting,
        false
    )

    ch_multiqc_files                  = ch_multiqc_files.mix(FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS.out.multiqc_files)
    ch_versions                       = ch_versions.mix(FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS.out.versions)
    ch_strand_inferred_filtered_fastq = FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS.out.reads
    ch_trim_read_count                = FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS.out.trim_read_count

    ch_trim_status = ch_trim_read_count
        .map {
            meta, num_reads ->
                return [ meta.id, num_reads > params.min_trimmed_reads.toFloat() ]
        }

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        ch_samplesheet
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'lncpipe_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    // AFTERQC(read_pairs_ch)
    // qc_ch = AFTERQC.out.html
    // read_pairs_ch = AFTERQC.out.reads



/*
* Step 1: Prepare Annotations
*/
/*
* Step 2: Build read aligner (STAR/tophat/HISAT2) index, if not provided
*/
/*
* Step 4: Initialize read alignment (STAR/HISAT2/tophat) <-- no tophat this time
*/
    // if (params.aligner == 'star') {
    //     STAR_ALIGN(read_pairs_ch, params.fasta, params.star_index)
    //     aligned_reads_ch = STAR_ALIGN.out.bam
    // } else if (params.aligner == 'hisat2') {
    //     HISAT2_ALIGN(read_pairs_ch, params.fasta, params.hisat2_index)
    //     aligned_reads_ch = HISAT2_ALIGN.out.bam
    // }


/*
* Step 5: Transcript assembly using Stringtie
*/
/*
* Step 6: Merged GTFs into one
*/
    // if (params.aligner == 'hisat2') {
    // STRINGTIE_ASSEMBLY(aligned_reads_ch, params.gtf)
    // STRINGTIE_MERGE(STRINGTIE_ASSEMBLY.out.gtf.collect(), params.fasta)
    // merged_gtf_ch = STRINGTIE_MERGE.out.merged_gtf
    // } else {
    // CUFFLINKS_ASSEMBLY(aligned_reads_ch, params.fasta, params.gtf)
    // CUFFMERGE(CUFFLINKS_ASSEMBLY.out.gtf.collect(), params.fasta)
    // merged_gtf_ch = CUFFMERGE.out.merged_gtf
    // }

/*
* Step 7: Compare assembled gtf with known annotations (GENCODE)
*/
    //GFFCOMPARE(merged_gtf_ch, params.gtf)


/*
* Step 8: Filter GTFs to distinguish novel lncRNAs
*/
    //IDENTIFY_NOVEL_LNCRNA(GFFCOMPARE.out.tmap, params.fasta, merged_gtf_ch)

/*
* Step 9: Predict coding potential abilities using CPAT and PLEK (CNCI functionality coming soon!)
*/
    // PREDICT_CODING_POTENTIAL_PLEK(IDENTIFY_NOVEL_LNCRNA.out.fasta)
    // PREDICT_CODING_POTENTIAL_CPAT(IDENTIFY_NOVEL_LNCRNA.out.fasta)

/*
* Step 9: Merged and filter lncRNAs based on coding potential (CPAT/PLEK)
*/
    // FILTER_LNCRNA(
    //     PREDICT_CODING_POTENTIAL_PLEK.out.results,
    //     PREDICT_CODING_POTENTIAL_CPAT.out.results,
    //     IDENTIFY_NOVEL_LNCRNA.out.exon_count,
    //     merged_gtf_ch,
    //     params.gtf
    // )

/*
* Step 10: Further filtered lncRNAs with known criterion
*/
    // SUMMARY_AND_CLASSIFICATION(
    //     FILTER_LNCRNA.out.novel_lncrna_gtf,
    //     params.gtf,
    //     params.fasta
    // )

/*
* Step 11: Rerun CPAT to evaluate the results
*/
    // RERUN_CPAT_LNCRNA(SUMMARY_AND_CLASSIFICATION.out.lncrna_fasta)
    // RERUN_CPAT_CODING(SUMMARY_AND_CLASSIFICATION.out.coding_fasta)
    // SECONDARY_STATISTICS(
    //     SUMMARY_AND_CLASSIFICATION.out.gtf,
    //     SUMMARY_AND_CLASSIFICATION.out.lncrna_gtf,
    //     RERUN_CPAT_LNCRNA.out.results,
    //     RERUN_CPAT_CODING.out.results,
    //     SUMMARY_AND_CLASSIFICATION.out.classification
    // )

/*
* Step 11: Quantification step (Kallisto/Htseq)
*/
/*
* Step 12: Generate count matrix for differential expression analysis
*/
    // if (params.quant == "htseq") {
    //     HTSEQ_COUNT(aligned_reads_ch, SUMMARY_AND_CLASSIFICATION.out.final_gtf)
    //     GET_HTSEQ_MATRIX(HTSEQ_COUNT.out.counts.collect(), SUMMARY_AND_CLASSIFICATION.out.final_gtf)
    //     expression_matrix_ch = GET_HTSEQ_MATRIX.out.matrix
    // } else {
    //     KALLISTO_INDEX(SUMMARY_AND_CLASSIFICATION.out.final_fasta)
    //     KALLISTO_QUANT(read_pairs_ch, KALLISTO_INDEX.out.index)
    //     GET_KALLISTO_MATRIX(KALLISTO_QUANT.out.abundance.collect(), SUMMARY_AND_CLASSIFICATION.out.final_gtf)
    //     expression_matrix_ch = GET_KALLISTO_MATRIX.out.matrix
    // }


/*
* Step 13: Perform Differential Expression analysis and generate report
*/
    // LNCPIPEREPORTER(
    //     params.design,
    //     STAR_ALIGN.out.log.collect(),
    //     SECONDARY_STATISTICS.out.stats,
    //     expression_matrix_ch
    // )


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]



}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
