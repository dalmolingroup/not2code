#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process HMMER_HMMSEARCH {
    tag "$meta.id"
    label 'process_high'
    
    container 'quay.io/biocontainers/hmmer:3.4--hb6cb901_4'
    
    publishDir "${params.outdir}/hmmer", mode: 'copy'
    
    input:
    tuple val(meta), path(longest_orfs_pep)
    val pfam_db_path
    
    output:
    tuple val(meta), path("*.domtblout"), emit: domtblout
    tuple val(meta), path("*.out"), emit: hmmsearch_out
    tuple val(meta), path("*.log"), emit: log
    path "Pfam-A.hmm*", emit: pfam_db, optional: true
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
    def pfam_url = "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"
    
    """
    # Log início do processo
    echo "Iniciando busca de domínios Pfam com HMMER" > ${prefix}.log
    echo "Arquivo de peptídeos: ${longest_orfs_pep}" >> ${prefix}.log
    echo "Caminho do banco Pfam especificado: ${pfam_db_path}" >> ${prefix}.log
    echo "E-value threshold: ${evalue}" >> ${prefix}.log
    echo "Número de threads: ${threads}" >> ${prefix}.log
    echo "Data/hora de início: \$(date)" >> ${prefix}.log
    
    # Verificar se o banco Pfam existe no caminho especificado
    PFAM_DB=""
    if [ -f "${pfam_db_path}" ]; then
        echo "Banco Pfam encontrado no caminho especificado: ${pfam_db_path}" >> ${prefix}.log
        PFAM_DB="${pfam_db_path}"
    elif [ -f "Pfam-A.hmm" ]; then
        echo "Banco Pfam encontrado no diretório de trabalho: Pfam-A.hmm" >> ${prefix}.log
        PFAM_DB="Pfam-A.hmm"
    else
        echo "Banco Pfam não encontrado. Iniciando download..." >> ${prefix}.log
        echo "URL de download: ${pfam_url}" >> ${prefix}.log
        
        # Download do banco Pfam
        echo "Baixando Pfam-A.hmm.gz..." >> ${prefix}.log
        wget -O Pfam-A.hmm.gz "${pfam_url}" 2>&1 | tee -a ${prefix}.log
        
        if [ \$? -eq 0 ] && [ -f "Pfam-A.hmm.gz" ]; then
            echo "Download concluído com sucesso" >> ${prefix}.log
            
            # Verificar tamanho do arquivo baixado
            file_size=\$(du -h Pfam-A.hmm.gz | cut -f1)
            echo "Tamanho do arquivo baixado: \$file_size" >> ${prefix}.log
            
            # Descompactar o arquivo
            echo "Descompactando Pfam-A.hmm.gz..." >> ${prefix}.log
            gunzip Pfam-A.hmm.gz 2>&1 | tee -a ${prefix}.log
            
            if [ \$? -eq 0 ] && [ -f "Pfam-A.hmm" ]; then
                echo "Descompactação concluída com sucesso" >> ${prefix}.log
                PFAM_DB="Pfam-A.hmm"
                
                # Verificar tamanho do arquivo descompactado
                file_size_unzip=\$(du -h Pfam-A.hmm | cut -f1)
                echo "Tamanho do arquivo descompactado: \$file_size_unzip" >> ${prefix}.log
            else
                echo "ERRO: Falha na descompactação do arquivo Pfam" >> ${prefix}.log
                exit 1
            fi
        else
            echo "ERRO: Falha no download do banco Pfam" >> ${prefix}.log
            exit 1
        fi
    fi
    
    echo "Usando banco Pfam: \$PFAM_DB" >> ${prefix}.log
    
    # Contar número de sequências peptídicas
    pep_count=\$(grep -c "^>" ${longest_orfs_pep})
    echo "Número de sequências peptídicas a serem analisadas: \$pep_count" >> ${prefix}.log
    
    # Verificar se o banco Pfam está pressionado
    if [ ! -f "\${PFAM_DB}.h3f" ]; then
        echo "Pressionando banco de dados Pfam..." >> ${prefix}.log
        hmmpress \$PFAM_DB 2>&1 | tee -a ${prefix}.log
        
        if [ \$? -eq 0 ]; then
            echo "Banco Pfam pressionado com sucesso" >> ${prefix}.log
            
            # Listar arquivos de índice criados
            echo "Arquivos de índice criados:" >> ${prefix}.log
            ls -la \${PFAM_DB}.* >> ${prefix}.log
        else
            echo "ERRO: Falha ao pressionar o banco Pfam" >> ${prefix}.log
            exit 1
        fi
    else
        echo "Banco de dados Pfam já está pressionado" >> ${prefix}.log
    fi
    
    # Verificar informações do banco Pfam
    echo "Informações do banco Pfam:" >> ${prefix}.log
    hmmstat \$PFAM_DB | head -10 >> ${prefix}.log
    
    # Executar hmmsearch
    echo "Executando hmmsearch..." >> ${prefix}.log
    hmmsearch \\
        --cpu ${threads} \\
        -E ${evalue} \\
        --domtblout ${domtblout_file} \\
        -o ${output_file} \\
        ${args} \\
        \$PFAM_DB \\
        ${longest_orfs_pep} \\
        2>&1 | tee -a ${prefix}.log
    
    # Verificar se os arquivos de saída foram criados
    if [ -f "${domtblout_file}" ] && [ -f "${output_file}" ]; then
        echo "HMMER hmmsearch concluído com sucesso" >> ${prefix}.log
        
        # Contar hits significativos
        significant_hits=\$(grep -v "^#" ${domtblout_file} | wc -l)
        echo "Número de hits significativos encontrados: \$significant_hits" >> ${prefix}.log
        
        # Contar domínios únicos
        unique_domains=\$(grep -v "^#" ${domtblout_file} | awk '{print \$4}' | sort | uniq | wc -l)
        echo "Número de domínios Pfam únicos identificados: \$unique_domains" >> ${prefix}.log
        
        # Contar sequências com hits
        sequences_with_hits=\$(grep -v "^#" ${domtblout_file} | awk '{print \$1}' | sort | uniq | wc -l)
        echo "Número de sequências com hits Pfam: \$sequences_with_hits" >> ${prefix}.log
        
        # Mostrar top 10 domínios mais frequentes
        if [ \$significant_hits -gt 0 ]; then
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
        echo "ERRO: Arquivos de saída não foram criados" >> ${prefix}.log
        exit 1
    fi
    
    echo "Data/hora de término: \$(date)" >> ${prefix}.log
    
    # Criar arquivo de versões
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hmmer: \$(hmmsearch -h | grep -oP 'HMMER \\K[0-9.]+' | head -1)
        wget: \$(wget --version | head -1 | grep -oP 'Wget \\K[0-9.]+')
    END_VERSIONS
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_pfam_hits_longest_orfs.domtblout
    touch ${prefix}_hmmsearch.out
    touch ${prefix}.log
    touch Pfam-A.hmm
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hmmer: "3.3.2"
        wget: "1.21.3"
    END_VERSIONS
    """
}