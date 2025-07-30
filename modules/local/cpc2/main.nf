#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process CPC2 {
    tag "$meta.id"
    label 'process_medium'
    
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'community.wave.seqera.io/library/cpc2:1.0.1--f139e3bfc23f65fe' :
        'community.wave.seqera.io/library/cpc2:1.0.1--f139e3bfc23f65fe' }"
    
    publishDir "${params.outdir}/cpc2", mode: 'copy'
    
    input:
    tuple val(meta), path(fasta_file)
    path cpc2_script
    
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
    echo "Script CPC2: ${cpc2_script}" >> ${prefix}.log
    
    # Contar número de sequências no arquivo FASTA
    seq_count=\$(grep -c "^>" ${fasta_file})
    echo "Número de sequências a serem analisadas: \$seq_count" >> ${prefix}.log
    
    # Copiar executáveis SVM para o diretório atual
    echo "Copiando executáveis SVM..." >> ${prefix}.log
    cp /opt/conda/bin/svm-scale ./svm-scale 2>&1 | tee -a ${prefix}.log
    cp /opt/conda/bin/svm-predict ./svm-predict 2>&1 | tee -a ${prefix}.log
    chmod +x ./svm-scale ./svm-predict
    
    # Verificar se os executáveis foram copiados
    echo "Verificando executáveis SVM:" >> ${prefix}.log
    ls -la ./svm-* >> ${prefix}.log
    
    # Adicionar diretório atual ao PATH
    export PATH="./:\$PATH"
    echo "PATH atualizado: \$PATH" >> ${prefix}.log
    
    # Criar script Python para modificar o CPC2.py corretamente
    cat > fix_cpc2_paths.py << 'EOF'
#!/usr/bin/env python3
import sys
import re

def fix_cpc2_script(input_file, output_file):
    """Corrigir caminhos dos executáveis SVM no script CPC2"""
    
    with open(input_file, 'r') as f:
        content = f.read()
    
    # Substituir referências aos executáveis SVM
    # Procurar por padrões como 'svm-scale' e 'svm-predict' e substituir por './svm-scale' e './svm-predict'
    
    # Padrão 1: chamadas diretas aos executáveis
    content = re.sub(r'(["\']?)svm-scale(["\']?)', r'\\1./svm-scale\\2', content)
    content = re.sub(r'(["\']?)svm-predict(["\']?)', r'\\1./svm-predict\\2', content)
    
    # Padrão 2: em comandos subprocess ou os.system
    content = re.sub(r'(subprocess\\.call\\(["\']?)svm-scale', r'\\1./svm-scale', content)
    content = re.sub(r'(subprocess\\.call\\(["\']?)svm-predict', r'\\1./svm-predict', content)
    content = re.sub(r'(os\\.system\\(["\']?)svm-scale', r'\\1./svm-scale', content)
    content = re.sub(r'(os\\.system\\(["\']?)svm-predict', r'\\1./svm-predict', content)
    
    # Padrão 3: em variáveis que armazenam o nome do executável
    content = re.sub(r'(=\\s*["\']?)svm-scale(["\']?)', r'\\1./svm-scale\\2', content)
    content = re.sub(r'(=\\s*["\']?)svm-predict(["\']?)', r'\\1./svm-predict\\2', content)
    
    with open(output_file, 'w') as f:
        f.write(content)
    
    print(f"Script CPC2 modificado: {input_file} -> {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Uso: python fix_cpc2_paths.py <input_file> <output_file>")
        sys.exit(1)
    
    fix_cpc2_script(sys.argv[1], sys.argv[2])
EOF
    
    # Executar o script de correção
    echo "Modificando script CPC2 para usar executáveis locais..." >> ${prefix}.log
    python fix_cpc2_paths.py ${cpc2_script} ./CPC2_fixed.py 2>&1 | tee -a ${prefix}.log
    
    # Verificar se o script foi criado
    if [ -f "./CPC2_fixed.py" ]; then
        echo "Script CPC2 modificado com sucesso" >> ${prefix}.log
        
        # Mostrar algumas linhas modificadas para verificação
        echo "Verificando modificações (primeiras ocorrências de svm):" >> ${prefix}.log
        grep -n "svm" ./CPC2_fixed.py | head -5 >> ${prefix}.log || echo "Nenhuma referência a svm encontrada" >> ${prefix}.log
    else
        echo "ERRO: Falha ao criar script CPC2 modificado" >> ${prefix}.log
        exit 1
    fi

    # Executar CPC2
    echo "Executando CPC2..." >> ${prefix}.log
    python ./CPC2_fixed.py \\
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
        cpc2: \$(python ./CPC2_fixed.py --version 2>&1 | grep -oP 'CPC2-\\K[0-9.]+' || echo "1.0.1")
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