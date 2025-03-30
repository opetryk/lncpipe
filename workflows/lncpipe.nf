/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { GFFCOMPARE                            } from '../modules/nf-core/gffcompare/main'
include { MULTIQC                               } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap                      } from 'plugin/nf-schema'
include { paramsSummaryMultiqc                  } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                } from '../subworkflows/local/utils_nfcore_lncpipe_pipeline'
include { samplesheetToList                     } from 'plugin/nf-schema'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS  } from '../subworkflows/nf-core/fastq_qc_trim_filter_setstrandedness'
include { FASTQ_ALIGN_STAR                      } from '../subworkflows/nf-core/fastq_align_star/main'
include { FASTQ_ALIGN_HISAT2                    } from '../subworkflows/nf-core/fastq_align_hisat2'
include { BAM_DEDUP_UMI as BAM_DEDUP_UMI_STAR   } from '../subworkflows/nf-core/bam_dedup_umi' // I would like to not import these as 2 and just have 1 workflow able to work with both star and hisat2 data.
include { BAM_DEDUP_UMI as BAM_DEDUP_UMI_HISAT2 } from '../subworkflows/nf-core/bam_dedup_umi'
include { STRINGTIE_WORKFLOW                    } from '../subworkflows/local/stringtie'
include { SUBREAD_FEATURECOUNTS                 } from '../../modules/nf-core/subread/featurecounts/main'
include { HTSEQ_COUNT                           } from '../modules/nf-core/htseq/count/main'  


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow LNCPIPE {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    ch_versions          // channel: [ path(versions.yml) ]
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
    ch_kallisto_index    // channel: [ meta, path(kallisto/index/) ] // this should not be ready yet..
    ch_bbsplit_index     // channel: path(bbsplit/index/)
    ch_ribo_db           // channel: path(sortmerna_fasta_list)
    ch_sortmerna_index   // channel: path(sortmerna/index/)
    ch_splicesites       // channel: path(genome.splicesites.txt)
    make_sortmerna_index // boolean: Whether to create an index before running sortmerna
    main:

    ch_multiqc_files = Channel.empty()

    //
    // Run RNA-seq FASTQ preprocessing subworkflow
    //

    // The subworkflow only has to do Salmon indexing if it discovers 'auto'
    // samples, and if we haven't already made one elsewhere
    salmon_index_available = params.salmon_index || (!params.skip_pseudo_alignment && params.pseudo_aligner == 'salmon')

    FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS (
        ch_samplesheet,
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

    /*
    * Step 4: Initialize read alignment (STAR/HISAT2/tophat) <-- no tophat this time
    */

    //
    // SUBWORKFLOW: Alignment with STAR and gene/transcript quantification with Salmon
    //

    ch_genome_bam          = Channel.empty()
    ch_genome_bam_index    = Channel.empty()
    ch_star_log            = Channel.empty()
    ch_unaligned_sequences = Channel.empty()
    ch_transcriptome_bam   = Channel.empty()

    if (!params.skip_alignment && params.aligner == 'star') {
        // Check if an AWS iGenome has been provided to use the appropriate version of STAR
        def is_aws_igenome = false
        if (params.fasta && params.gtf) {
            if ((file(params.fasta).getName() - '.gz' == 'genome.fa') && (file(params.gtf).getName() - '.gz' == 'genes.gtf')) {
                is_aws_igenome = true
            }
        }

        FASTQ_ALIGN_STAR(
            FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS.out.reads,
            ch_star_index.map { [ [:], it ] },
            ch_gtf.map { [ [:], it ] },
            params.star_ignore_sjdbgtf,
            '',
            params.seq_center ?: '',
            ch_fasta.map { [ [:], it ] },
            ch_transcript_fasta.map { [ [:], it ] }
        )

        ch_genome_bam              = FASTQ_ALIGN_STAR.out.bam
        ch_genome_bam_index        = FASTQ_ALIGN_STAR.out.bai
        ch_transcriptome_bam       = FASTQ_ALIGN_STAR.out.orig_bam_transcript
        ch_transcriptome_bai       = FASTQ_ALIGN_STAR.out.bai_transcript
        ch_versions                = ch_versions.mix(FASTQ_ALIGN_STAR.out.versions)

        ch_multiqc_files = ch_multiqc_files
            .mix(FASTQ_ALIGN_STAR.out.stats.collect{it[1]})
            .mix(FASTQ_ALIGN_STAR.out.flagstat.collect{it[1]})
            .mix(FASTQ_ALIGN_STAR.out.idxstats.collect{it[1]})
            .mix(FASTQ_ALIGN_STAR.out.log_final.collect{it[1]})

        if (params.bam_csi_index) {
            ch_genome_bam_index = FASTQ_ALIGN_STAR.out.csi
        }
        ch_versions = ch_versions.mix(FASTQ_ALIGN_STAR.out.versions)

        //
        // SUBWORKFLOW: Remove duplicate reads from BAM file based on UMIs
        //
        if (params.with_umi) {

            BAM_DEDUP_UMI_STAR(
                ch_genome_bam.join(ch_genome_bam_index, by: [0]),
                ch_fasta.map { [ [:], it ] },
                params.umi_dedup_tool,
                params.umitools_dedup_stats,
                params.bam_csi_index,
                ch_transcriptome_bam,
                ch_transcript_fasta.map { [ [:], it ] }
            )

            ch_genome_bam        = BAM_DEDUP_UMI_STAR.out.bam
            ch_transcriptome_bam = BAM_DEDUP_UMI_STAR.out.transcriptome_bam
            ch_genome_bam_index  = BAM_DEDUP_UMI_STAR.out.bai
            ch_versions          = ch_versions.mix(BAM_DEDUP_UMI_STAR.out.versions)

            ch_multiqc_files = ch_multiqc_files
                .mix(BAM_DEDUP_UMI_STAR.out.multiqc_files)

        } else {
            // The deduplicated stats should take priority for MultiQC, but use
            // them straight out of the aligner otherwise

            ch_multiqc_files = ch_multiqc_files
                .mix(FASTQ_ALIGN_STAR.out.stats.collect{it[1]})
                .mix(FASTQ_ALIGN_STAR.out.flagstat.collect{it[1]})
                .mix(FASTQ_ALIGN_STAR.out.idxstats.collect{it[1]})
        }

    }

    //
    // SUBWORKFLOW: Alignment with HISAT2
    //
    if (!params.skip_alignment && params.aligner == 'hisat2') {
        FASTQ_ALIGN_HISAT2 (
            ch_strand_inferred_filtered_fastq,
            ch_hisat2_index.map { [ [:], it ] },
            ch_splicesites.map { [ [:], it ] },
            ch_fasta.map { [ [:], it ] }
        )
        ch_genome_bam          = FASTQ_ALIGN_HISAT2.out.bam
        ch_genome_bam_index    = FASTQ_ALIGN_HISAT2.out.bai
        ch_unaligned_sequences = FASTQ_ALIGN_HISAT2.out.fastq
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_HISAT2.out.summary.collect{it[1]})

        if (params.bam_csi_index) {
            ch_genome_bam_index = FASTQ_ALIGN_HISAT2.out.csi
        }
        ch_versions = ch_versions.mix(FASTQ_ALIGN_HISAT2.out.versions)

        //
        // SUBWORKFLOW: Remove duplicate reads from BAM file based on UMIs
        //

        if (params.with_umi) {

            BAM_DEDUP_UMI_HISAT2(
                ch_genome_bam.join(ch_genome_bam_index, by: [0]),
                ch_fasta.map { [ [:], it ] },
                params.umi_dedup_tool,
                params.umitools_dedup_stats,
                params.bam_csi_index,
                ch_transcriptome_bam,
                ch_transcript_fasta.map { [ [:], it ] }
            )

            ch_genome_bam        = BAM_DEDUP_UMI_HISAT2.out.bam
            ch_genome_bam_index  = BAM_DEDUP_UMI_HISAT2.out.bai
            ch_versions          = ch_versions.mix(BAM_DEDUP_UMI_HISAT2.out.versions)

            ch_multiqc_files = ch_multiqc_files
                .mix(BAM_DEDUP_UMI_HISAT2.out.multiqc_files)
        } else {

            // The deduplicated stats should take priority for MultiQC, but use
            // them straight out of the aligner otherwise
            ch_multiqc_files = ch_multiqc_files
                .mix(FASTQ_ALIGN_HISAT2.out.stats.collect{it[1]})
                .mix(FASTQ_ALIGN_HISAT2.out.flagstat.collect{it[1]})
                .mix(FASTQ_ALIGN_HISAT2.out.idxstats.collect{it[1]})
        }
    }

    //
    // Filter channels to get samples that passed STAR minimum mapping percentage
    //
    if (!params.skip_alignment && params.aligner.contains('star')) {
        ch_star_log
            .map { meta, align_log -> [ meta ] + getStarPercentMapped(params, align_log) }
            .set { ch_percent_mapped }

        // Save status for workflow summary
        ch_map_status = ch_percent_mapped
            .map {
                meta, mapped, pass ->
                    return [ meta.id, pass ]
            }

        ch_percent_mapped
            .branch { meta, mapped, pass ->
                pass: pass
                    return [ "$meta.id\t$mapped" ]
                fail: !pass
                    return [ "$meta.id\t$mapped" ]
            }
            .set { ch_pass_fail_mapped }

        ch_pass_fail_mapped
            .fail
            .collect()
            .map {
                tsv_data ->
                    def header = ["Sample", "STAR uniquely mapped reads (%)"]
                    sample_status_header_multiqc.text + multiqcTsvFromList(tsv_data, header)
            }
            .set { ch_fail_mapping_multiqc }
        ch_multiqc_files = ch_multiqc_files.mix(ch_fail_mapping_multiqc.collectFile(name: 'fail_mapped_samples_mqc.tsv'))
    }


/*
* Step 5: Transcript assembly using Stringtie and merge gtf into one
*/
    STRINGTIE_WORKFLOW (
        ch_genome_bam,
        ch_gtf
    )
    ch_versions = ch_versions.mix(STRINGTIE_WORKFLOW.out.versions)
    ch_merged_gtf = STRINGTIE_WORKFLOW.out.stringtie_gtf_merged
/*
* Step 6: Compare assembled gtf with known annotations
*/
    ch_fasta_meta_fai = ch_fasta.combine(ch_fai)
        .map { fasta, fai ->
            def meta2 = [ id: fasta.getBaseName(), description: 'Genome FASTA with index' ]
            return [ meta2, fasta, fai ]
        }

    ch_reference_gtf = ch_gtf.map { gtf_file ->
        def meta3 = [ id: gtf_file.getBaseName(), description: 'Reference GTF file' ]
        return [ meta3, gtf_file ]
    }
    GFFCOMPARE(ch_merged_gtf, ch_fasta_meta_fai, ch_reference_gtf)
    ch_versions = ch_versions.mix(GFFCOMPARE.out.versions)

    ch_multiqc_files = ch_multiqc_files
        .mix(GFFCOMPARE.out.stats.collect{it[1]})

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
* Step 11: Quantification step (Featurecounts/Htseq)
*/

 if (params.counts == 'featurecounts') {
    ch_feature_counts = ch_genome_bam.combine(ch_gtf) 
 
        SUBREAD_FEATURECOUNTS(
            ch_feature_counts
            )
        ch_versions         = ch_versions.mix(SUBREAD_FEATURECOUNTS.out.versions)
        ch_counts           = SUBREAD_FEATURECOUNTS.out.counts
    }

    if (params.counts == 'htseq') {
        ch_genome_bam_index = ch_genome_bam.combine(ch_genome_bam_index)
        ch_htseq_counts     = ch_genome_bam.combine(ch_genome_bam_index).combine(ch_gtf)
        HTSEQ_COUNT(
            ch_htseq_counts
        )
        ch_versions         = ch_versions.mix(HTSEQ_COUNT.out.versions)
        ch_counts           = HTSEQ_COUNT.out.txt
    }

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
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'lncpipe_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }



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
