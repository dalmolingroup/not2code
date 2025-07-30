#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process COMPARE_TRANSCRIPTOMES {
    tag "$meta.id"
    label 'process_medium'
    
    container 'docker://rocker/tidyverse:4.3.0'
    
    input:
    tuple val(meta), 
    path(tmap_file), 
    path(refmap_file), 
    path(annotated_gtf)
    
    output:
    tuple val(meta), path("*.remove_prot.gtf"), emit: filtered_gtf
    tuple val(meta), path("*.log"),             emit: log
    path "versions.yml",                        emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    #!/usr/bin/env Rscript
    
    # Carregar bibliotecas necessárias
    library(tidyverse)
    library(rtracklayer)
    
    # Configurar log
    log_file <- "${prefix}.log"
    sink(log_file, append = TRUE, split = TRUE)
    
    cat("Iniciando análise de comparação de transcriptomas\\n")
    cat("Arquivos de entrada:\\n")
    cat("  TMAP: ${tmap_file}\\n")
    cat("  REFMAP: ${refmap_file}\\n")
    cat("  GTF anotado: ${annotated_gtf}\\n")
    
    # Ler arquivos de entrada
    cat("Lendo arquivo TMAP...\\n")
    dup_reference_tmap <- read.table("${tmap_file}", header = TRUE)
    
    cat("Lendo arquivo REFMAP...\\n")
    dup_reference_refmap <- read.table("${refmap_file}", header = TRUE)
    
    cat("Lendo arquivo GTF anotado...\\n")
    dup_reference_annot <- as.data.frame(rtracklayer::import("${annotated_gtf}"))
    
    # Análise dos códigos de classe
    cat("Distribuição dos códigos de classe:\\n")
    class_table <- table(dup_reference_tmap\$class_code)
    print(class_table)
    
    # Verificar se todos os ref_ids com class_code "=" estão no refmap
    equal_refs <- dup_reference_tmap %>% 
        filter(class_code == "=") %>% 
        pull(ref_id)
    
    verification <- all(equal_refs %in% dup_reference_refmap\$ref_id)
    cat("Verificação - todos os ref_ids com class_code '=' estão no refmap:", verification, "\\n")
    
    # Identificar classes diferentes de "="
    classes_different <- dup_reference_tmap %>% 
        filter(class_code != "=") %>% 
        pull(qry_id)
    
    cat("Número de transcritos com class_code diferente de '=':", length(classes_different), "\\n")
    
    # Verificar se o número de classes diferentes corresponde à diferença entre tmap e refmap
    expected_diff <- nrow(dup_reference_tmap) - nrow(dup_reference_refmap)
    actual_diff <- length(classes_different)
    cat("Diferença esperada:", expected_diff, "\\n")
    cat("Diferença atual:", actual_diff, "\\n")
    cat("Verificação da diferença:", actual_diff == expected_diff, "\\n")
    
    # Filtrar GTF removendo proteínas codificantes
    cat("Filtrando GTF para remover proteínas codificantes...\\n")
    gtf_remove_prot <- dup_reference_annot %>% 
        filter(transcript_id %in% classes_different) %>% 
        filter(strand %in% c("+", "-"))
    
    # Contar transcritos únicos
    unique_transcripts <- gtf_remove_prot %>% 
        filter(type == "transcript") %>% 
        pull(transcript_id) %>% 
        unique()
    
    cat("Número de transcritos únicos após filtragem:", length(unique_transcripts), "\\n")
    
    # Exportar GTF filtrado
    output_file <- "${prefix}.remove_prot.gtf"
    cat("Exportando GTF filtrado para:", output_file, "\\n")
    
    rtracklayer::export(gtf_remove_prot, output_file, format = "gtf")
    
    cat("Análise concluída com sucesso!\\n")
    sink()
    
    # Criar arquivo de versões
    cat('"${task.process}":\\n', file="versions.yml")
    cat('    r-base: "', R.version.string, '"\\n', file="versions.yml", append=TRUE)
    cat('    tidyverse: "', packageVersion("tidyverse"), '"\\n', file="versions.yml", append=TRUE)
    cat('    rtracklayer: "', packageVersion("rtracklayer"), '"\\n', file="versions.yml", append=TRUE)
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.remove_prot.gtf
    touch ${prefix}.log
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: "4.3.0"
        tidyverse: "2.0.0"
        rtracklayer: "1.60.0"
    END_VERSIONS
    """
}
