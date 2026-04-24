process REPORT {
    tag "$meta.id"
    label 'process_single'
    
    // We use a container that has R, tidyverse, and quarto. 
    // We'll install rtracklayer and plotly dynamically if they aren't present.
    // rocker/verse has tidyverse and quarto pre-installed.
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'community.wave.seqera.io/library/bioconductor-rtracklayer_quarto_r-base_r-ggplot2_pruned:6894bafbf31d3979' :
        'community.wave.seqera.io/library/bioconductor-rtracklayer_quarto_r-base_r-ggplot2_pruned:6894bafbf31d3979' }"
    
    publishDir "${params.outdir}/report", mode: 'copy'
    
    input:
    tuple val(meta), path(gtf_complete)
    tuple val(meta2), path(combined_gtf)
    path quarto_script
    
    output:
    path "*.html", emit: report
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def min_length = task.ext.min_length ?: '200'
    
    """
    # Install required packages if not present in rocker/verse
    #Rscript -e 'if (!requireNamespace("plotly", quietly = TRUE)) install.packages("plotly", repos="http://cran.us.r-project.org")'
    #Rscript -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos="http://cran.us.r-project.org")'
    #Rscript -e 'if (!requireNamespace("rtracklayer", quietly = TRUE)) BiocManager::install("rtracklayer")'

    # Export environmental variables for Quarto to pick up
    export GTF_COMPLETE="${gtf_complete}"
    export COMBINED_GTF="${combined_gtf}"
    export MIN_LENGTH="${min_length}"
    
    # Render the Quarto document
    quarto render ${quarto_script} --to html --output not2code_report_${meta.id}.html
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quarto: \$(quarto --version)
        r-base: \$(Rscript -e 'cat(as.character(getRversion()))')
    END_VERSIONS
    """
}
