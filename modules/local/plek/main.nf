#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process PLEK {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'community.wave.seqera.io/library/plek_numpy_scikit-learn_scipy:ebea7cd00ea7ee71' :
    'community.wave.seqera.io/library/plek_numpy_scikit-learn_scipy:ebea7cd00ea7ee71' }"

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
    def threads = task.cpus ?: 16
    
    """
    # Log início do processo
    echo "Iniciando análise de potencial codificante com PLEK" > ${prefix}.log
    echo "Arquivo FASTA de entrada: ${fasta_file}" >> ${prefix}.log
    echo "Script PLEK: ${plek_script}" >> ${prefix}.log
    echo "Arquivo de saída: ${output_file}" >> ${prefix}.log
    echo "Número de threads: ${threads}" >> ${prefix}.log
    echo "Data/hora de início: \$(date)" >> ${prefix}.log
    
    # Contar número de sequências no arquivo FASTA
    seq_count=\$(grep -c "^>" ${fasta_file})
    echo "Número de sequências a serem analisadas: \$seq_count" >> ${prefix}.log
    
    # Verificar se o script PLEK existe
    if [ ! -f "${plek_script}" ]; then
        echo "ERRO: Script PLEK não encontrado: ${plek_script}" >> ${prefix}.log
        exit 1
    fi
    
    # Executar PLEK
    echo "Executando PLEK..." >> ${prefix}.log
    python ${plek_script} \
        -fasta ${fasta_file} \
        -out ${output_file} \
        -thread ${threads} \
        ${args} \
        2>&1 | tee -a ${prefix}.log
    
    # Verificar se o arquivo de saída foi criado
    if [ -f "${output_file}" ]; then
        echo "Análise PLEK concluída com sucesso" >> ${prefix}.log
        
        # Contar resultados
        total_results=\$(tail -n +2 ${output_file} | wc -l)
        coding_count=\$(tail -n +2 ${output_file} | awk '\$2=="Coding"' | wc -l)
        noncoding_count=\$(tail -n +2 ${output_file} | awk '\$2=="Non-coding"' | wc -l)
        
        echo "Total de transcritos analisados: \$total_results" >> ${prefix}.log
        echo "Transcritos codificantes: \$coding_count" >> ${prefix}.log
        echo "Transcritos não-codificantes: \$noncoding_count" >> ${prefix}.log
        
        # Mostrar estatísticas dos scores
        echo "Estatísticas dos scores PLEK:" >> ${prefix}.log
        tail -n +2 ${output_file} | awk '{print \$3}' | sort -n | awk '
        BEGIN { sum = 0; count = 0; }
        { 
            values[count] = \$1; 
            sum += \$1; 
            count++; 
        }
        END {
            if (count > 0) {
                mean = sum / count;
                if (count % 2 == 1) {
                    median = values[int(count/2)];
                } else {
                    median = (values[count/2-1] + values[count/2]) / 2;
                }
                print "  Score médio: " mean;
                print "  Score mediano: " median;
                print "  Score mínimo: " values[0];
                print "  Score máximo: " values[count-1];
            }
        }' >> ${prefix}.log
        
    else
        echo "ERRO: Arquivo de saída não foi criado" >> ${prefix}.log
        exit 1
    fi
    
    echo "Data/hora de término: \$(date)" >> ${prefix}.log
    
    # Criar arquivo de versões
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plek: "1.2"
        python: \$(python --version 2>&1 | grep -oP 'Python \\K[0-9.]+')
    END_VERSIONS
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_PLEK_output.txt
    touch ${prefix}.log
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plek: "1.2"
        python: "3.10.5"
    END_VERSIONS
    """
}