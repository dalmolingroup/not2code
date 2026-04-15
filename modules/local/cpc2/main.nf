process CPC2 {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::cpc2=1.0.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cpc2:1.0.1--hdfd78af_0' :
        'quay.io/biocontainers/cpc2:1.0.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.txt"), emit: cpc2_results
    tuple val(meta), path("*.log"), emit: log
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Log início do processo
    echo "Iniciando análise de potencial codificante com CPC2" > ${prefix}.log
    echo "Arquivo FASTA de entrada: ${fasta}" >> ${prefix}.log
    echo "Arquivo de saída: ${prefix}_CPC2_output.txt" >> ${prefix}.log
    echo "Data/hora de início: \$(date)" >> ${prefix}.log

    # Contar número de sequências no arquivo FASTA
    seq_count=\$(grep -c "^>" ${fasta})
    echo "Número de sequências a serem analisadas: \$seq_count" >> ${prefix}.log

    # CPC2 instalado via bioconda - disponível diretamente no PATH como CPC2.py
    CPC2.py \\
        -i ${fasta} \\
        -o ${prefix} \\
        ${args} \\
        2>> ${prefix}.log

    echo "Data/hora de término: \$(date)" >> ${prefix}.log

    # Criar arquivo de versões
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cpc2: \$(CPC2.py --version 2>&1 | grep -oP 'CPC2-\\K[0-9.]+' || echo "1.0.1")
        python: \$(python --version 2>&1 | grep -oP 'Python \\K[0-9.]+')
        libsvm: \$(svm-predict 2>&1 | head -1 | grep -oP 'libsvm version \\K[0-9.]+' || echo "3.25")
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
        python: "3.9"
        libsvm: "3.25"
    END_VERSIONS
    """
}