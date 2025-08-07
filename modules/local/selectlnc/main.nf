#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process SELECT_LNCRNAS {
    tag "$meta.id"
    label 'process_medium'
    
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'community.wave.seqera.io/library/bioconductor-rtracklayer_r-tidyverse:4e377b50350447f9' :
        'community.wave.seqera.io/library/bioconductor-rtracklayer_r-tidyverse:4e377b50350447f9' }"
    
    publishDir "${params.outdir}/lncrnas", mode: 'copy'
    
    input:
    tuple val(meta), path(annotated_gtf)
    tuple val(meta2), path(cpc2_results)
    tuple val(meta3), path(plek_results)
    tuple val(meta4), path(pfam_domtblout)
    path reference_gtf
    path lncselect_script  // Adicionar o script como input
    
    output:
    tuple val(meta), path("gtf_lncRNA.gtf"), emit: lncrna_gtf
    tuple val(meta), path("gtf_complete.gtf"), emit: complete_gtf
    tuple val(meta), path("combined_stringtie_ncbi.gtf"), emit: combined_gtf
    tuple val(meta), path("combined_stringtie_ncbi_new_names.gtf"), emit: final_gtf
    //tuple val(meta), path("*.log"), emit: log
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def min_length = task.ext.min_length ?: '200'
    def pfam_evalue = task.ext.pfam_evalue ?: '1e-10'
    
    """
    # Executar o script R com as variáveis de ambiente
    export ANNOTATED_GTF="${annotated_gtf}"
    export CPC2_RESULTS="${cpc2_results}"
    export PLEK_RESULTS="${plek_results}"
    export PFAM_DOMTBLOUT="${pfam_domtblout}"
    export REFERENCE_GTF="${reference_gtf}"
    export PREFIX="${prefix}"
    export MIN_LENGTH="${min_length}"
    export PFAM_EVALUE="${pfam_evalue}"
    
    Rscript ${lncselect_script}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch gtf_lncRNA.gtf
    touch gtf_complete.gtf
    touch combined_stringtie_ncbi.gtf
    touch combined_stringtie_ncbi_new_names.gtf
    touch ${prefix}.log
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: "4.3.0"
        tidyverse: "2.0.0"
        rtracklayer: "1.60.0"
    END_VERSIONS
    """
}