process GTF_FILTER_TPM {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::gawk=5.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/gawk:5.1.0--h5b5514e_1' :
    'quay.io/biocontainers/gawk:5.1.0--h5b5514e_1' }"

    input:
    tuple val(meta), path(gtf)
    val(tpm_threshold)

    output:
    tuple val(meta), path("*_filter_TPM.gtf"), emit: filtered_gtf
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def threshold = tpm_threshold ?: 0.1
    """
    # Create temporary file for valid transcript IDs
    valid_ids=\$(mktemp)

    # Extract transcript_id from transcripts with TPM >= threshold
    awk '\$3 == "transcript" && match(\$0, /transcript_id "([^"]+)".*TPM "([0-9.]+)"/, a) { 
        if (a[2] >= ${threshold}) print a[1]; 
    }' "${gtf}" > "\$valid_ids"

    # Filter GTF and keep transcripts and their exons
    awk -v valid_ids="\$valid_ids" '
    BEGIN {
        while ((getline < valid_ids) > 0) valid[\$1] = 1;
        close(valid_ids);
    }
    \$3 == "transcript" && match(\$0, /transcript_id "([^"]+)"/, a) { 
        keep = (a[1] in valid); 
    }
    keep || (\$3 == "exon" && match(\$0, /transcript_id "([^"]+)"/, a) && (a[1] in valid))
    ' "${gtf}" > "${prefix}_filter_TPM.gtf"

    # Remove temporary file
    rm -f "\$valid_ids"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(awk --version | head -n1 | sed 's/GNU Awk //; s/,.*//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_filter_TPM.gtf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(awk --version | head -n1 | sed 's/GNU Awk //; s/,.*//' || echo "5.1.0")
    END_VERSIONS
    """
}
