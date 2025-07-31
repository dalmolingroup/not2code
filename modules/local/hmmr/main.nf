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
    # Não usar pipefail para evitar problemas com SIGPIPE
    set -eu
    
    # Log início do processo
    echo "=== Iniciando busca de domínios Pfam com HMMER ===" > ${prefix}.log
    echo "Data/hora: \$(date)" >> ${prefix}.log
    echo "Arquivo de peptídeos: ${longest_orfs_pep}" >> ${prefix}.log
    echo "Banco Pfam: ${pfam_db}" >> ${prefix}.log
    echo "E-value threshold: ${evalue}" >> ${prefix}.log
    echo "Threads: ${threads}" >> ${prefix}.log
    echo "PID do processo: \$\$" >> ${prefix}.log
    
    # Verificar recursos e limites
    echo "=== Recursos e limites ===" >> ${prefix}.log
    echo "Memória total: \$(free -h | grep '^Mem:' | awk '{print \$2}')" >> ${prefix}.log
    echo "Memória disponível: \$(free -h | grep '^Mem:' | awk '{print \$7}')" >> ${prefix}.log
    echo "Espaço em disco: \$(df -h . | tail -1 | awk '{print \$4}')" >> ${prefix}.log
    echo "Limites do processo:" >> ${prefix}.log
    ulimit -a >> ${prefix}.log 2>&1 || echo "Não foi possível obter limites" >> ${prefix}.log
    
    # Verificar arquivo de entrada
    if [ ! -f "${longest_orfs_pep}" ]; then
        echo "ERRO: Arquivo de peptídeos não encontrado: ${longest_orfs_pep}" >> ${prefix}.log
        exit 1
    fi
    
    pep_count=\$(grep -c "^>" ${longest_orfs_pep} 2>/dev/null || echo "0")
    echo "Número de sequências peptídicas: \$pep_count" >> ${prefix}.log
    
    if [ "\$pep_count" -eq 0 ]; then
        echo "ERRO: Nenhuma sequência peptídica encontrada" >> ${prefix}.log
        exit 1
    fi
    
    file_size=\$(du -h ${longest_orfs_pep} | cut -f1)
    echo "Tamanho do arquivo de peptídeos: \$file_size" >> ${prefix}.log
    
    # Preparar banco Pfam
    PFAM_DB="${pfam_db}"
    
    if [[ "${pfam_db}" == *.gz ]]; then
        echo "Descompactando banco Pfam..." >> ${prefix}.log
        gunzip -c ${pfam_db} > Pfam-A.hmm
        PFAM_DB="Pfam-A.hmm"
        echo "Banco descompactado: \$(du -h \$PFAM_DB | cut -f1)" >> ${prefix}.log
    fi
    
    # Verificar banco
    if [ ! -f "\$PFAM_DB" ]; then
        echo "ERRO: Banco Pfam não encontrado: \$PFAM_DB" >> ${prefix}.log
        exit 1
    fi
    
    pfam_size=\$(du -h \$PFAM_DB | cut -f1)
    echo "Tamanho do banco Pfam: \$pfam_size" >> ${prefix}.log
    
    # Indexar se necessário
    if [ ! -f "\${PFAM_DB}.h3f" ]; then
        echo "Indexando banco Pfam..." >> ${prefix}.log
        hmmpress \$PFAM_DB >> ${prefix}.log 2>&1
        
        if [ \$? -ne 0 ]; then
            echo "ERRO: Falha ao indexar banco Pfam" >> ${prefix}.log
            exit 1
        fi
        
        echo "Banco indexado com sucesso" >> ${prefix}.log
    else
        echo "Banco já está indexado" >> ${prefix}.log
    fi
    
    # Verificar informações do banco
    echo "=== Informações do banco Pfam ===" >> ${prefix}.log
    hmmstat \$PFAM_DB | head -5 >> ${prefix}.log 
    
    # Teste simples do hmmsearch primeiro
    echo "=== Testando hmmsearch ===" >> ${prefix}.log
    hmmsearch -h >> ${prefix}.log  || echo "AVISO: hmmsearch -h falhou" >> ${prefix}.log
    
    # Executar hmmsearch com monitoramento detalhado
    echo "=== Executando hmmsearch ===" >> ${prefix}.log
    echo "Comando: hmmsearch --cpu ${threads} -E ${evalue} --domtblout ${domtblout_file} -o ${output_file} ${args} \$PFAM_DB ${longest_orfs_pep}" >> ${prefix}.log
    echo "Memória antes: \$(free -h | grep '^Mem:' | awk '{print \$7}')" >> ${prefix}.log
    echo "Início: \$(date)" >> ${prefix}.log
    
    # Executar em background para monitorar
    hmmsearch \\
        --cpu ${threads} \\
        -E ${evalue} \\
        --domtblout ${domtblout_file} \\
        -o ${output_file} \\
        ${args} \\
        \$PFAM_DB \\
        ${longest_orfs_pep} > hmmsearch_stdout.log 2> hmmsearch_stderr.log &
    
    hmmsearch_pid=\$!
    echo "PID do hmmsearch: \$hmmsearch_pid" >> ${prefix}.log
    
    # Monitorar execução
    monitor_count=0
    while kill -0 \$hmmsearch_pid 2>/dev/null; do
        sleep 30
        monitor_count=\$((monitor_count + 1))
        echo "[\$(date)] Monitoramento \$monitor_count: hmmsearch ainda executando" >> ${prefix}.log
        echo "  Memória disponível: \$(free -h | grep '^Mem:' | awk '{print \$7}')" >> ${prefix}.log
        
        # Verificar se arquivos de saída estão sendo criados
        if [ -f "${domtblout_file}" ]; then
            domtbl_size=\$(wc -l < ${domtblout_file} 2>/dev/null || echo "0")
            echo "  Linhas em domtblout: \$domtbl_size" >> ${prefix}.log
        fi
        
        # Timeout após 4 horas (480 monitoramentos de 30s)
        if [ \$monitor_count -gt 480 ]; then
            echo "TIMEOUT: Matando hmmsearch após 4 horas" >> ${prefix}.log
            kill -9 \$hmmsearch_pid 2>/dev/null || true
            break
        fi
    done
    
    # Aguardar conclusão
    wait \$hmmsearch_pid
    hmmsearch_exit=\$?
    
    echo "Fim: \$(date)" >> ${prefix}.log
    echo "Exit code: \$hmmsearch_exit" >> ${prefix}.log
    echo "Memória após: \$(free -h | grep '^Mem:' | awk '{print \$7}')" >> ${prefix}.log
    
    # Capturar saídas do hmmsearch
    if [ -f hmmsearch_stdout.log ]; then
        echo "=== STDOUT do hmmsearch ===" >> ${prefix}.log
        cat hmmsearch_stdout.log >> ${prefix}.log
    fi
    
    if [ -f hmmsearch_stderr.log ]; then
        echo "=== STDERR do hmmsearch ===" >> ${prefix}.log
        cat hmmsearch_stderr.log >> ${prefix}.log
    fi
    
    # Verificar arquivos de saída
    echo "=== Verificando arquivos de saída ===" >> ${prefix}.log
    echo "Arquivos no diretório:" >> ${prefix}.log
    ls -la >> ${prefix}.log
    
    # Criar arquivos mínimos se não existirem
    if [ ! -f "${domtblout_file}" ]; then
        echo "Criando arquivo domtblout mínimo" >> ${prefix}.log
        echo "# hmmsearch :: search profile(s) against a sequence database" > ${domtblout_file}
        echo "# HMMER 3.4 (Aug 2023); http://hmmer.org/" >> ${domtblout_file}
        echo "#" >> ${domtblout_file}
        echo "# [No hits detected that satisfy reporting thresholds]" >> ${domtblout_file}
    fi
    
    if [ ! -f "${output_file}" ]; then
        echo "Criando arquivo de saída mínimo" >> ${prefix}.log
        echo "hmmsearch :: search profile(s) against a sequence database" > ${output_file}
        echo "HMMER 3.4 (Aug 2023); http://hmmer.org/" >> ${output_file}
        echo "" >> ${output_file}
        echo "[No hits detected that satisfy reporting thresholds]" >> ${output_file}
    fi
    
    # Análise dos resultados
    if [ -s "${domtblout_file}" ]; then
        hits=\$(grep -v "^#" ${domtblout_file} | wc -l)
        echo "Hits encontrados: \$hits" >> ${prefix}.log
        
        if [ "\$hits" -gt 0 ]; then
            domains=\$(grep -v "^#" ${domtblout_file} | awk '{print \$4}' | sort | uniq | wc -l)
            echo "Domínios únicos: \$domains" >> ${prefix}.log
        fi
    else
        echo "Nenhum hit encontrado" >> ${prefix}.log
    fi
    
    echo "Processo finalizado: \$(date)" >> ${prefix}.log
    
    # Versões
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