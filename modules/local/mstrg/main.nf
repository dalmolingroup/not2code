#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process MSTRG_PREP {
    tag "$meta.id"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/perl:5.26--2' :
    'quay.io/biocontainers/perl:5.26--2' }"
    
    input:
    tuple val(meta), path(gtf_file)
    path perl_script
    
    output:
    tuple val(meta), path("*_prep.gtf"), emit: gtf_prep
    tuple val(meta), path("*.log"),      emit: log
    path "versions.yml",                 emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output_name = "${prefix}_prep.gtf"
    
    """
    # Log início do processo
    echo "Iniciando processamento de nomes de transcritos com script Perl" > ${prefix}.log
    echo "Arquivo GTF de entrada: ${gtf_file}" >> ${prefix}.log
    echo "Script Perl: ${perl_script}" >> ${prefix}.log
    
    # Executar script Perl para preparar nomes dos transcritos
    perl ${perl_script} ${gtf_file} > ${output_name} 2>&1
    
    # Verificar se o arquivo foi criado com sucesso
    if [ -f "${output_name}" ]; then
        echo "Processamento concluído com sucesso" >> ${prefix}.log
        echo "Arquivo de saída: ${output_name}" >> ${prefix}.log
        
        # Contar número de linhas processadas
        input_lines=\$(wc -l < ${gtf_file})
        output_lines=\$(wc -l < ${output_name})
        echo "Linhas no arquivo de entrada: \$input_lines" >> ${prefix}.log
        echo "Linhas no arquivo de saída: \$output_lines" >> ${prefix}.log
    else
        echo "ERRO: Falha na criação do arquivo de saída" >> ${prefix}.log
        exit 1
    fi
    
    # Criar arquivo de versões
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(perl --version | grep -oP 'v\\K[0-9.]+')
    END_VERSIONS
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_prep.gtf
    touch ${prefix}.log
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: "5.32.1"
    END_VERSIONS
    """
}