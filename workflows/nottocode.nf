#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

//include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { SAMPLESHEET_CHECK      } from '../modules/local/samplesheet_check/main'
include { GTF_FILTER_TPM         } from '../modules/local/gtf_filter_tpm/main'
include { STRINGTIE_MERGE        } from '../modules/nf-core/stringtie/merge/main'
include { GFFCOMPARE             } from '../modules/nf-core/gffcompare/main'
include { COMPARE_TRANSCRIPTOMES } from '../modules/local/compare_transcriptomes/main'
include { GFFREAD                } from '../modules/nf-core/gffread/main'
include { MSTRG_PREP             } from '../modules/local/mstrg/main'
include { CPC2                   } from '../modules/local/cpc2/main'
include { PLEK                   } from '../modules/local/plek/main'
include { TRANSDECODER_LONGORFS  } from '../modules/local/transdecoder/main'
include { HMMER_HMMPRESS         } from '../modules/nf-core/hmmer/hmmpress/main'
include { HMMER_HMMSEARCH        } from '../modules/nf-core/hmmer/hmmsearch/main'
include { SELECT_LNCRNAS         } from '../modules/local/selectlnc/main'
include { XZ_DECOMPRESS          } from '../modules/nf-core/xz/decompress/main'  
include { UNTAR                  } from '../modules/nf-core/untar/main'  
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-validation'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_nottocode_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow NOTTOCODE {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // Verificar parâmetros obrigatórios
    //
    if (!params.input) {
        error "ERROR: missing --input samplesheet.csv"
    }

    if (!params.reference_genome) {
        error "ERROR: missing --reference_genome fasta"
    }

    if (!params.reference_gtf) {
        error "ERROR: missing --reference_gtf annotation file"
    }

    if (!params.pfam_db) {
        error "ERROR: missing --pfam_db pfam file. you can download it with: wget https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz then gunzip Pfam-A.hmm.gz"
    }

    ch_plek_archive = Channel.fromPath(params.plek, checkIfExists: true)
    .map { file -> 
        def meta = [id: 'plek_db']
        [meta, file]
    }

    XZ_DECOMPRESS (
        ch_plek_archive
    )

    UNTAR (
        XZ_DECOMPRESS.out.file
    )
    ch_plek_db = UNTAR.out.untar.map { meta, path -> path }

    //
    // Validate samplesheet
    //
    SAMPLESHEET_CHECK(ch_samplesheet)
    ch_versions = ch_versions.mix(SAMPLESHEET_CHECK.out.versions)

    //
    // Criar canal de entrada a partir da samplesheet
    //
    ch_input = SAMPLESHEET_CHECK.out.csv
        .splitCsv(header:true, sep:',')
        .map { row ->
            def meta = [:]
            meta.id = row.sample
            meta.single_end = true
            
            def gtf_file = file(row.gtf)
            if (!gtf_file.exists()) {
                error("ERROR: Arquivo GTF não encontrado: ${row.gtf}")
            }
            
            return [meta, gtf_file]
        }

    //
    // Refence files
    //
    ch_genome_fasta = file(params.reference_genome)
    ch_reference_gtf = file(params.reference_gtf)

    //
    // Filter GTF by TPM
    //
    if (params.filter_by_tpm) {
        GTF_FILTER_TPM(
            ch_input,
            params.tpm_threshold
        )
        ch_filtered_gtf = GTF_FILTER_TPM.out.filtered_gtf
        ch_versions = ch_versions.mix(GTF_FILTER_TPM.out.versions)
    } else {
        ch_filtered_gtf = ch_input
    }

    ch_filtered_gtf.view { meta, gtf_file ->
        def status = params.filter_by_tpm ? "TPM-filtered" : "original"
        log.info "Using ${status} GTF file: ${gtf_file} for sample: ${meta.id}"
    }

    //
    // Preparar GTFs para StringTie Merge
    //
    ch_gtf_files = ch_filtered_gtf
        .map { meta, gtf -> gtf }  // Extrair apenas os arquivos GTF
        .collect()                 // Coletar todos em uma lista

    // Debug opcional
    ch_gtf_files.view { gtf_list ->
        log.info "StringTie Merge will process ${gtf_list.size()} GTF files: ${gtf_list}"
    }

    

    //
    // StringTie Merge all GTFs
    //
    
    STRINGTIE_MERGE (
        ch_gtf_files,    
        ch_reference_gtf    
    )
    ch_versions = ch_versions.mix(STRINGTIE_MERGE.out.versions)


    //
    // Preapare gffcompare entry
    //

    ch_merged_gtf = STRINGTIE_MERGE.out.gtf.map { gtf ->
        def meta = [id: 'merged_all_samples']
        [meta, gtf]
    }

    //
    // Compare, merge, annotate and estimate accuracy of generated gtf files0
    //

    GFFCOMPARE (
        ch_merged_gtf,
        ch_reference_gtf
    )
    ch_versions = ch_versions.mix(GFFCOMPARE.out.versions)

    //
    // Debug: Ver saídas do GFFCOMPARE
    //

    GFFCOMPARE.out.tmap.view { "TMAP: ${it}" }
    GFFCOMPARE.out.refmap.view { "REFMAP: ${it}" }
    GFFCOMPARE.out.annotated_gtf.view { "ANNOTATED_GTF: ${it}" }

    //
    // Combinar saídas do GFFCOMPARE
    //

    ch_gffcompare_combined = GFFCOMPARE.out.tmap
        .join(GFFCOMPARE.out.refmap, by: 0)      // Join por meta (índice 0)
        .join(GFFCOMPARE.out.annotated_gtf, by: 0)
        .map { meta, tmap, refmap, annotated_gtf ->
            // Verificar se todos os arquivos existem
            log.info "Combining for ${meta.id}: tmap=${tmap}, refmap=${refmap}, annotated=${annotated_gtf}"
            return [meta, tmap, refmap, annotated_gtf]
        }

    //
    // Debug: Ver canal combinado
    //
    ch_gffcompare_combined.view { "Combined channel: ${it}" }


    //
    // Transcipts compare and coding proteins remove
    //

    COMPARE_TRANSCRIPTOMES (
        ch_gffcompare_combined
    )
    ch_versions = ch_versions.mix(COMPARE_TRANSCRIPTOMES.out.versions)

    //
    // Validate, filter, convert and perform various other operations on GFF
    //

    GFFREAD (
        COMPARE_TRANSCRIPTOMES.out.filtered_gtf,
        ch_genome_fasta
    )
    ch_versions = ch_versions.mix(GFFREAD.out.versions)

    //
    // Process transcipts names
    //

    MSTRG_PREP (
        COMPARE_TRANSCRIPTOMES.out.filtered_gtf,
        params.mstrg_prep_script
    )
    ch_versions = ch_versions.mix(MSTRG_PREP.out.versions)

    //
    // Protein coding potential
    //

    CPC2 (
        GFFREAD.out.gffread_fasta,
        file(params.cpc2)
    )
    ch_versions = ch_versions.mix(CPC2.out.versions)


    //
    // PLEK - Coding potential analysis
    //

    PLEK (
        GFFREAD.out.gffread_fasta,
        ch_plek_db,
    )
    ch_versions = ch_versions.mix(PLEK.out.versions)

    //
    // TransDecoder LongOrfs
    //

    TRANSDECODER_LONGORFS (
        GFFREAD.out.gffread_fasta
    )
    ch_versions = ch_versions.mix(TRANSDECODER_LONGORFS.out.versions)   

    //
    // HMMER HMMSEARCH
    //

    ch_pfam_db = Channel.fromPath(params.pfam_db)

    HMMER_HMMPRESS (
        ch_pfam_db
    )
    ch_versions = ch_versions.mix(HMMER_HMMPRESS.out.versions)

    // Definir os parâmetros booleanos
    def write_align = false
    def write_target = false  
    def write_domain = true

    // Combinar os canais corretamente
    // HMMER_HMMPRESS não retorna meta, então precisamos criar o canal corretamente
    // Precisamos do arquivo HMM original, não apenas dos índices
    
    ch_hmmsearch_input = TRANSDECODER_LONGORFS.out.longest_orfs_pep
        .map { meta, seqdb ->
            [meta, params.pfam_db, seqdb, write_align, write_target, write_domain]
        }

    HMMER_HMMSEARCH (
        ch_hmmsearch_input
    )
    ch_versions = ch_versions.mix(HMMER_HMMSEARCH.out.versions)

    SELECT_LNCRNAS (
        MSTRG_PREP.out.gtf_prep,
        CPC2.out.cpc2_results,
        PLEK.out.plek_results,
        HMMER_HMMSEARCH.out.domain_summary,
        ch_reference_gtf
    )
     ch_versions = ch_versions.mix(SELECT_LNCRNAS.out.versions)


    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_pipeline_software_mqc_versions.yml',
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

    // Add the missing channels
    ch_multiqc_replace_names = params.multiqc_replace_names ?
        Channel.fromPath(params.multiqc_replace_names, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_sample_names  = params.multiqc_sample_names ?
        Channel.fromPath(params.multiqc_sample_names, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))

    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
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
        ch_multiqc_replace_names.toList(),  // Added
        ch_multiqc_sample_names.toList()    // Added
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
