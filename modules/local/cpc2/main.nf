#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process CPC2 {
    tag "$meta.id"
    label 'process_medium'
    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/cpc2:1.0.1--hdfd78af_0' :
    'quay.io/biocontainers/cpc2:1.0.1--hdfd78af_0' }"
    
    publishDir "${params.outdir}/cpc2", mode: 'copy'
    
    input:
    tuple val(meta), path(fasta_file)
    
    output:
    tuple val(meta), path("*.txt"), emit: cpc2_results
    tuple val(meta), path("*.log"), emit: log
    path "versions.yml",            emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output_file = "${prefix}_CPC2_output.txt"
    
    """
    # Log início do processo
    echo "Iniciando análise de potencial codificante com CPC2" > ${prefix}.log
    echo "Arquivo FASTA de entrada: ${fasta_file}" >> ${prefix}.log
    echo "Arquivo de saída: ${output_file}" >> ${prefix}.log
    echo "Data/hora de início: \$(date)" >> ${prefix}.log
    
    # Contar número de sequências no arquivo FASTA
    seq_count=\$(grep -c "^>" ${fasta_file})
    echo "Número de sequências a serem analisadas: \$seq_count" >> ${prefix}.log
    
    # Executar CPC2
    echo "Executando CPC2..." >> ${prefix}.log
    CPC2.py \\
        -i ${fasta_file} \\
        -o ${output_file} \\
        ${args} \\
        2>&1 | tee -a ${prefix}.log
    
    # Verificar se o arquivo de saída foi criado
    if [ -f "${output_file}" ]; then
        echo "Análise CPC2 concluída com sucesso" >> ${prefix}.log
        
        # Contar resultados
        total_results=\$(tail -n +2 ${output_file} | wc -l)
        coding_count=\$(tail -n +2 ${output_file} | awk '\$8=="coding"' | wc -l)
        noncoding_count=\$(tail -n +2 ${output_file} | awk '\$8=="noncoding"' | wc -l)
        
        echo "Total de transcritos analisados: \$total_results" >> ${prefix}.log
        echo "Transcritos codificantes: \$coding_count" >> ${prefix}.log
        echo "Transcritos não-codificantes: \$noncoding_count" >> ${prefix}.log
    else
        echo "ERRO: Arquivo de saída não foi criado" >> ${prefix}.log
        exit 1
    fi
    
    echo "Data/hora de término: \$(date)" >> ${prefix}.log
    
    # Criar arquivo de versões
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cpc2: \$(CPC2.py --version 2>&1 | grep -oP 'CPC2-\\K[0-9.]+' || echo "1.0.1")
        python: \$(python --version 2>&1 | grep -oP 'Python \\K[0-9.]+')
    END_VERSIONS
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_CPC2_output.txt
    touch ${prefix}.log
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cpc2: "1.0.1"
        python: "2.7.18"
    END_VERSIONS
    """
}