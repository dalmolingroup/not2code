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

    # CPC2.py hardcodes svm-predict/svm-scale at ../libs/libsvm/libsvm-3.18/ relative
    # to its own location (/usr/local/bin/). This path is non-writable when Docker runs
    # as the host user (-u uid:gid, set in nextflow.config docker.runOptions).
    # It also imports seqio as a sibling Python module from /usr/local/bin/.
    #
    # Fix: copy CPC2.py to the work dir, patch lib_dir to a writable local cpc2_libs/
    # subdir (with symlinks to svm-predict/svm-scale), patch data_dir to the real
    # /usr/local/data/ (where cpc2.model and cpc2.range live), and export PYTHONPATH
    # so seqio.py can be found.
    mkdir -p cpc2_libs/libsvm/libsvm-3.18/
    ln -sf \$(which svm-predict) cpc2_libs/libsvm/libsvm-3.18/svm-predict
    ln -sf \$(which svm-scale)   cpc2_libs/libsvm/libsvm-3.18/svm-scale

    cp \$(which CPC2.py) ./CPC2_patched.py
    chmod +x ./CPC2_patched.py
    ABS_LIBS=\$(realpath cpc2_libs)

    sed -i "s|lib_dir = script_dir + \\"/../libs/\\"|lib_dir = '\${ABS_LIBS}/'|"   ./CPC2_patched.py
    sed -i 's|data_dir = script_dir + "/../data/"|data_dir = "/usr/local/data/"|'   ./CPC2_patched.py

    PYTHONPATH=/usr/local/bin python ./CPC2_patched.py \\
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