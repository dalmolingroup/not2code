#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

//include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { SAMPLESHEET_CHECK      } from '../modules/local/samplesheet_check/main'
include { GTF_FILTER_TPM         } from '../modules/local/gtf_filter_tpm/main'
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
    genome_fasta = file(params.reference_genome)
    reference_gtf = file(params.reference_gtf)

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
        ch_multiqc_logo.toList()
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
