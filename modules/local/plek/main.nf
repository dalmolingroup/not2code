#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process PLEK {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'community.wave.seqera.io/library/cxx-compiler_python:25b9eddb2d84e6d6' :
    'community.wave.seqera.io/library/cxx-compiler_python:25b9eddb2d84e6d6' }"

    publishDir "${params.outdir}/plek", mode: 'copy'
    
    input:
    tuple val(meta), path(fasta_file)
    path plek_script
    
    output:
    tuple val(meta), path("*.txt"), emit: plek_results
    tuple val(meta), path("*.log"), emit: log
    path "versions.yml",            emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output_file = "${prefix}_PLEK_output.txt"
    def threads = task.cpus ?: 6
    
    """
    # Inicializar log
    echo "=== Análise PLEK ===" > ${prefix}.log
    echo "Data/hora: \$(date)" >> ${prefix}.log
    echo "Arquivo FASTA: ${fasta_file}" >> ${prefix}.log
    echo "Arquivo de saída: ${output_file}" >> ${prefix}.log
    echo "Threads: ${threads}" >> ${prefix}.log
    
    # Verificar arquivo FASTA
    seq_count=\$(grep -c "^>" ${fasta_file} 2>/dev/null || echo "0")
    echo "Número de sequências: \$seq_count" >> ${prefix}.log
    
    if [ "\$seq_count" -eq 0 ]; then
        echo "ERRO: Nenhuma sequência no FASTA" >> ${prefix}.log
        exit 1
    fi
    
    # Executar PLEK - ignorar exit code e verificar sucesso pelo arquivo de saída
    echo "Executando PLEK..." >> ${prefix}.log
    
    python ${plek_script}/PLEK.py \\
        -fasta ${fasta_file} \\
        -out ${output_file} \\
        -thread ${threads} \\
        ${args} >> ${prefix}.log 2>&1 || true
    
    # Verificar se o PLEK executou com sucesso baseado no arquivo de saída
    if [ ! -f "${output_file}" ]; then
        echo "ERRO: PLEK falhou - arquivo de saída não foi criado" >> ${prefix}.log
        exit 1
    fi
    
    if [ ! -s "${output_file}" ]; then
        echo "ERRO: PLEK falhou - arquivo de saída está vazio" >> ${prefix}.log
        exit 1
    fi
    
    # Verificar se o arquivo tem conteúdo válido (pelo menos 2 linhas)
    line_count=\$(wc -l < ${output_file} 2>/dev/null || echo "0")
    if [ "\$line_count" -lt 2 ]; then
        echo "ERRO: PLEK falhou - arquivo de saída inválido (\$line_count linhas)" >> ${prefix}.log
        exit 1
    fi
    
    echo "PLEK executado com sucesso" >> ${prefix}.log
    echo "Arquivo de saída: ${output_file} (\$line_count linhas)" >> ${prefix}.log
    
    # Análise básica dos resultados
    total_results=\$(tail -n +2 ${output_file} | wc -l)
    coding_count=\$(tail -n +2 ${output_file} | grep -c "^Coding" || echo "0")
    noncoding_count=\$(tail -n +2 ${output_file} | grep -c "^Non-coding" || echo "0")
    
    echo "=== Resultados ===" >> ${prefix}.log
    echo "Total de transcritos analisados: \$total_results" >> ${prefix}.log
    echo "Transcritos codificantes: \$coding_count" >> ${prefix}.log
    echo "Transcritos não-codificantes: \$noncoding_count" >> ${prefix}.log
    
    # Calcular porcentagem usando uma abordagem mais simples
    if [ "\$total_results" -gt 0 ]; then
        # Usar Python para calcular a porcentagem (mais confiável)
        coding_percent=\$(python3 -c "print(f'{(\$coding_count * 100.0 / \$total_results):.2f}')")
        echo "Percentual codificante: \$coding_percent%" >> ${prefix}.log
    fi
    
    echo "Análise concluída com sucesso" >> ${prefix}.log
    
    # Criar arquivo de versões
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plek: "1.2"
        python: \$(python --version 2>&1 | grep -oP 'Python \\K[0-9.]+' || echo "unknown")
    END_VERSIONS
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "STUB mode" > ${prefix}.log
    echo -e "Transcript_ID\\tLabel\\tScore" > ${prefix}_PLEK_output.txt
    echo -e "test_transcript_1\\tCoding\\t0.8" >> ${prefix}_PLEK_output.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plek: "1.2"
        python: "3.10.5"
    END_VERSIONS
    """
}