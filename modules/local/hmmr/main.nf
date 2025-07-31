#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process HMMER_HMMSEARCH {
    tag "$meta.id"
    label 'process_high'
    
    container 'quay.io/biocontainers/hmmer:3.4--hb6cb901_4'
    
    publishDir "${params.outdir}/hmmer", mode: 'copy'
    
    input:
    tuple val(meta), path(longest_orfs_pep)
    path pfam_db
    
    output:
    tuple val(meta), path("*.domtblout"), emit: domtblout
    tuple val(meta), path("*.out"), emit: hmmsearch_out
    tuple val(meta), path("*.log"), emit: log
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def domtblout_file = "${prefix}_pfam_hits_longest_orfs.domtblout"
    def output_file = "${prefix}_hmmsearch.out"
    def evalue = task.ext.evalue ?: '1e-10'
    def threads = task.cpus ?: 1
    
    """
    set -euo pipefail
    
    # Log início do processo
    echo "=== Iniciando busca de domínios Pfam com HMMER ===" > ${prefix}.log
    echo "Data/hora: \$(date)" >> ${prefix}.log
    echo "Arquivo de peptídeos: ${longest_orfs_pep}" >> ${prefix}.log
    echo "Banco Pfam: ${pfam_db}" >> ${prefix}.log
    echo "E-value threshold: ${evalue}" >> ${prefix}.log
    echo "Threads: ${threads}" >> ${prefix}.log
    
    # Verificar recursos disponíveis
    echo "=== Recursos disponíveis ===" >> ${prefix}.log
    echo "Memória: \$(free -h | grep '^Mem:' | awk '{print \$2 " total, " \$7 " disponível"}')" >> ${prefix}.log
    echo "Espaço em disco: \$(df -h . | tail -1 | awk '{print \$4 " disponível"}')" >> ${prefix}.log
    
    # Verificar arquivo de entrada (peptídeos)
    if [ ! -f "${longest_orfs_pep}" ]; then
        echo "ERRO: Arquivo de peptídeos não encontrado: ${longest_orfs_pep}" >> ${prefix}.log
        exit 1
    fi
    
    pep_count=\$(grep -c "^>" ${longest_orfs_pep} 2>/dev/null || echo "0")
    echo "Número de sequências peptídicas: \$pep_count" >> ${prefix}.log
    
    if [ "\$pep_count" -eq 0 ]; then
        echo "ERRO: Nenhuma sequência peptídica encontrada no arquivo" >> ${prefix}.log
        exit 1
    fi
    
    file_size=\$(du -h ${longest_orfs_pep} | cut -f1)
    echo "Tamanho do arquivo de peptídeos: \$file_size" >> ${prefix}.log
    
    # Verificar banco Pfam
    if [ ! -f "${pfam_db}" ]; then
        echo "ERRO: Banco Pfam não encontrado: ${pfam_db}" >> ${prefix}.log
        exit 1
    fi
    
    pfam_size=\$(du -h ${pfam_db} | cut -f1)
    echo "Tamanho do banco Pfam: \$pfam_size" >> ${prefix}.log
    
    # Verificar se o banco está indexado (hmmpress)
    if [ ! -f "${pfam_db}.h3f" ]; then
        echo "Indexando banco Pfam com hmmpress..." >> ${prefix}.log
        hmmpress ${pfam_db} >> ${prefix}.log 2>&1
        
        if [ \$? -ne 0 ]; then
            echo "ERRO: Falha ao indexar banco Pfam" >> ${prefix}.log
            exit 1
        fi
        
        echo "Banco Pfam indexado com sucesso" >> ${prefix}.log
        echo "Arquivos de índice criados:" >> ${prefix}.log
        ls -la ${pfam_db}.* >> ${prefix}.log
    else
        echo "Banco Pfam já está indexado" >> ${prefix}.log
    fi
    
    # Verificar informações do banco Pfam
    echo "=== Informações do banco Pfam ===" >> ${prefix}.log
    hmmstat ${pfam_db} | head -10 >> ${prefix}.log 2>&1
    
    # Executar hmmsearch
    echo "=== Executando hmmsearch ===" >> ${prefix}.log
    echo "Comando: hmmsearch --cpu ${threads} -E ${evalue} --domtblout ${domtblout_file} -o ${output_file} ${args} ${pfam_db} ${longest_orfs_pep}" >> ${prefix}.log
    echo "Início da execução: \$(date)" >> ${prefix}.log
    
    # Executar com timeout de 4 horas
    timeout 14400 hmmsearch \\
        --cpu ${threads} \\
        -E ${evalue} \\
        --domtblout ${domtblout_file} \\
        -o ${output_file} \\
        ${args} \\
        ${pfam_db} \\
        ${longest_orfs_pep} >> ${prefix}.log 2>&1
    
    hmmsearch_exit=\$?
    echo "Fim da execução: \$(date)" >> ${prefix}.log
    echo "Exit code do hmmsearch: \$hmmsearch_exit" >> ${prefix}.log
    
    # Verificar resultado da execução
    if [ \$hmmsearch_exit -eq 124 ]; then
        echo "ERRO: hmmsearch foi interrompido por timeout (4 horas)" >> ${prefix}.log
        exit 1
    elif [ \$hmmsearch_exit -ne 0 ]; then
        echo "ERRO: hmmsearch falhou (exit code: \$hmmsearch_exit)" >> ${prefix}.log
        exit 1
    fi
    
    # Verificar se os arquivos de saída foram criados
    if [ ! -f "${domtblout_file}" ] || [ ! -f "${output_file}" ]; then
        echo "ERRO: Arquivos de saída não foram criados" >> ${prefix}.log
        echo "Arquivos presentes no diretório:" >> ${prefix}.log
        ls -la >> ${prefix}.log
        exit 1
    fi
    
    echo "=== hmmsearch concluído com sucesso ===" >> ${prefix}.log
    
    # Análise dos resultados
    if [ -s "${domtblout_file}" ]; then
        significant_hits=\$(grep -v "^#" ${domtblout_file} | wc -l)
        echo "Hits significativos encontrados: \$significant_hits" >> ${prefix}.log
        
        if [ "\$significant_hits" -gt 0 ]; then
            unique_domains=\$(grep -v "^#" ${domtblout_file} | awk '{print \$4}' | sort | uniq | wc -l)
            sequences_with_hits=\$(grep -v "^#" ${domtblout_file} | awk '{print \$1}' | sort | uniq | wc -l)
            
            echo "Domínios Pfam únicos identificados: \$unique_domains" >> ${prefix}.log
            echo "Sequências com hits Pfam: \$sequences_with_hits" >> ${prefix}.log
            
            # Top 10 domínios mais frequentes
            echo "Top 10 domínios Pfam mais frequentes:" >> ${prefix}.log
            grep -v "^#" ${domtblout_file} | awk '{print \$4}' | sort | uniq -c | sort -nr | head -10 >> ${prefix}.log
            
            # Estatísticas dos E-values
            echo "Estatísticas dos E-values:" >> ${prefix}.log
            grep -v "^#" ${domtblout_file} | awk '{print \$7}' | sort -g | awk '
            BEGIN { count = 0; }
            { 
                evalues[count] = \$1; 
                count++; 
            }
            END {
                if (count > 0) {
                    if (count % 2 == 1) {
                        median = evalues[int(count/2)];
                    } else {
                        median = (evalues[count/2-1] + evalues[count/2]) / 2;
                    }
                    print "  E-value mediano: " median;
                    print "  E-value mínimo: " evalues[0];
                    print "  E-value máximo: " evalues[count-1];
                }
            }' >> ${prefix}.log
        fi
    else
        echo "Nenhum hit significativo encontrado" >> ${prefix}.log
    fi
    
    echo "Processo finalizado com sucesso: \$(date)" >> ${prefix}.log
    
    # Criar arquivo de versões
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hmmer: \$(hmmsearch -h | grep -oP 'HMMER \\K[0-9.]+' | head -1)
    END_VERSIONS
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_pfam_hits_longest_orfs.domtblout
    touch ${prefix}_hmmsearch.out
    touch ${prefix}.log
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hmmer: "3.3.2"
    END_VERSIONS
    """
}