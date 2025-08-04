#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process TRANSDECODER_LONGORFS {
    tag "$meta.id"
    label 'process_medium'
    
    container 'quay.io/biocontainers/transdecoder:5.5.0--pl5321hdfd78af_5'
    
    publishDir "${params.outdir}/transdecoder", mode: 'copy'
    
    input:
    tuple val(meta), path(fasta_file)
    
    output:
    tuple val(meta), path("${prefix}.transdecoder_dir/"), emit: transdecoder_dir
    tuple val(meta), path("${prefix}.transdecoder_dir/longest_orfs.pep"), emit: longest_orfs_pep
    tuple val(meta), path("${prefix}.transdecoder_dir/longest_orfs.gff3"), emit: longest_orfs_gff3
    tuple val(meta), path("${prefix}.transdecoder_dir/longest_orfs.cds"), emit: longest_orfs_cds
    tuple val(meta), path("*.log"), emit: log
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def output_dir = "${prefix}.transdecoder_dir"
    
    """
    # Log início do processo
    echo "Iniciando identificação de ORFs longas com TransDecoder.LongOrfs" > ${prefix}.log
    echo "Arquivo FASTA de entrada: ${fasta_file}" >> ${prefix}.log
    echo "Diretório de saída: ${output_dir}" >> ${prefix}.log
    echo "Data/hora de início: \$(date)" >> ${prefix}.log
    
    # Contar número de sequências no arquivo FASTA
    seq_count=\$(grep -c "^>" ${fasta_file})
    echo "Número de sequências a serem analisadas: \$seq_count" >> ${prefix}.log
    
    # Executar TransDecoder.LongOrfs
    echo "Executando TransDecoder.LongOrfs..." >> ${prefix}.log
    TransDecoder.LongOrfs \\
        -t ${fasta_file} \\
        --output_dir ${output_dir} \\
        ${args} \\
        2>&1 | tee -a ${prefix}.log
    
    # Verificar se os arquivos de saída foram criados
    if [ -d "${output_dir}" ]; then
        echo "TransDecoder.LongOrfs concluído com sucesso" >> ${prefix}.log
        
        # Contar ORFs identificadas
        if [ -f "${output_dir}/longest_orfs.pep" ]; then
            orf_count=\$(grep -c "^>" ${output_dir}/longest_orfs.pep)
            echo "Número de ORFs longas identificadas: \$orf_count" >> ${prefix}.log
        fi
        
        # Listar arquivos criados
        echo "Arquivos criados no diretório TransDecoder:" >> ${prefix}.log
        ls -la ${output_dir}/ >> ${prefix}.log
        
        # Estatísticas dos comprimentos das ORFs
        if [ -f "${output_dir}/longest_orfs.pep" ]; then
            echo "Estatísticas dos comprimentos das ORFs:" >> ${prefix}.log
            grep -v "^>" ${output_dir}/longest_orfs.pep | awk '{print length(\$0)}' | sort -n | awk '
            BEGIN { sum = 0; count = 0; }
            { 
                lengths[count] = \$1; 
                sum += \$1; 
                count++; 
            }
            END {
                if (count > 0) {
                    mean = sum / count;
                    if (count % 2 == 1) {
                        median = lengths[int(count/2)];
                    } else {
                        median = (lengths[count/2-1] + lengths[count/2]) / 2;
                    }
                    print "  Comprimento médio: " mean " aa";
                    print "  Comprimento mediano: " median " aa";
                    print "  Comprimento mínimo: " lengths[0] " aa";
                    print "  Comprimento máximo: " lengths[count-1] " aa";
                }
            }' >> ${prefix}.log
        fi
        
    else
        echo "ERRO: Diretório de saída não foi criado" >> ${prefix}.log
        exit 1
    fi
    
    echo "Data/hora de término: \$(date)" >> ${prefix}.log
    
    # Criar arquivo de versões
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        transdecoder: \$(TransDecoder.LongOrfs --version 2>&1 | grep -oP 'TransDecoder_v\\K[0-9.]+' || echo "5.5.0")
    END_VERSIONS
    """
    
    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}.transdecoder_dir
    touch ${prefix}.transdecoder_dir/longest_orfs.pep
    touch ${prefix}.transdecoder_dir/longest_orfs.gff3
    touch ${prefix}.transdecoder_dir/longest_orfs.cds
    touch ${prefix}.log
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        transdecoder: "5.5.0"
    END_VERSIONS
    """
}