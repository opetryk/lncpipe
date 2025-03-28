include { STRINGTIE_STRINGTIE }    from '../../../modules/nf-core/stringtie/stringtie/main'
include { STRINGTIE_MERGE }        from '../../../modules/nf-core/stringtie/merge/main'


workflow STRINGTIE_WORKFLOW {
    take:
        bam_sorted
        ch_gtf

    main:
        ch_versions = Channel.empty()
        ch_stringtie_gtf = Channel.empty()


        //
        // STRINGTIE: Transcript assembly and quantification for RNA-SeQ
        //
        STRINGTIE_STRINGTIE(bam_sorted, ch_gtf)
        ch_stringtie_gtf = STRINGTIE_STRINGTIE.out.transcript_gtf.map { meta, transcript_gtf -> [ transcript_gtf ] }.collect()
        ch_versions = ch_versions.mix(STRINGTIE_STRINGTIE.out.versions)

        //
        // STRINGTIE: Merge transcript assemblies into a non-redundant annotation
        //
        STRINGTIE_MERGE (ch_stringtie_gtf, ch_gtf)
        ch_stringtie_gtf_merged = STRINGTIE_MERGE.out.gtf
            .map { gtf_file ->
                def meta = [ id: 'stringtie_merged' ]
                return [ meta, gtf_file ]
            }
        ch_versions = ch_versions.mix(STRINGTIE_MERGE.out.versions)

    emit:
        stringtie_gtf_merged = ch_stringtie_gtf_merged
        versions      = ch_versions

    }
