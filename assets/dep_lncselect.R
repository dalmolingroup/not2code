 #!/usr/bin/env Rscript
    
    # Carregar bibliotecas necessárias
    library(tidyverse)
    library(rtracklayer)

    # Obter variáveis de ambiente
    annotated_gtf <- Sys.getenv("ANNOTATED_GTF")
    cpc2_results <- Sys.getenv("CPC2_RESULTS")
    plek_results <- Sys.getenv("PLEK_RESULTS")
    pfam_domtblout <- Sys.getenv("PFAM_DOMTBLOUT")
    reference_gtf <- Sys.getenv("REFERENCE_GTF")
    prefix <- Sys.getenv("PREFIX")
    min_length <- as.numeric(Sys.getenv("MIN_LENGTH"))
    pfam_evalue <- as.numeric(Sys.getenv("PFAM_EVALUE"))

    # Configurar log
    log_file <- paste0(prefix, ".log")
    sink(log_file, append = TRUE, split = TRUE)
    
    cat("=== SELEÇÃO DE lncRNAs ===\n")
    cat("Iniciando seleção de lncRNAs e combinação com anotação de referência\n")
    cat("Data/hora de início:", as.character(Sys.time()), "\n")
    cat("\n")

    # Ler arquivos de entrada
    cat("Lendo arquivos de entrada...\n")
    cat("  GTF anotado:", annotated_gtf, "\n")
    cat("  Resultados CPC2:", cpc2_results, "\n")
    cat("  Resultados PLEK:", plek_results, "\n")
    cat("  Domínios Pfam:", pfam_domtblout, "\n")
    cat("  GTF de referência:", reference_gtf, "\n")
    cat("\n")

    # Carregar GTF anotado
    dup_reference_annot <- as.data.frame(rtracklayer::import(annotated_gtf))
    cat("GTF anotado carregado:", nrow(dup_reference_annot), "linhas\n")

    # Carregar resultados CPC2
    CPC2 <- read.delim(cpc2_results)
    colnames(CPC2) <- c("transcript_id", "transcript_length", "peptide_length", "Fickett_score",
                        "pI", "ORF_integrity", "coding_probability", "CPC2")
    cat("Resultados CPC2 carregados:", nrow(CPC2), "transcritos\n")
    cat("  Codificantes:", sum(CPC2$CPC2 == "coding"), "\n")
    cat("  Não-codificantes:", sum(CPC2$CPC2 == "noncoding"), "\n")

    # Carregar resultados PLEK
    PLEK <- read.delim(plek_results, header = FALSE)
    PLEK <- PLEK %>% mutate(transcript_id = str_extract(V3, "(?<=\\>).*?(?= loc:)"))
    colnames(PLEK) <- c("PLEK", "value", "features", "transcript_id")
    cat("Resultados PLEK carregados:", nrow(PLEK), "transcritos\n")
    cat("  Codificantes:", sum(PLEK$PLEK == "Coding"), "\n")
    cat("  Não-codificantes:", sum(PLEK$PLEK == "Non-coding"), "\n")

    # Carregar domínios Pfam
    pfam <- readr::read_table(pfam_domtblout, comment = "#", col_names = FALSE)
    cat("Domínios Pfam carregados:", nrow(pfam), "hits\n")

    # Identificar transcritos com domínios significativos
    transcripts_with_domain <- as.data.frame(pfam) %>%
    dplyr::filter(X7 <= pfam_evalue) %>%
    mutate(transcript_id = str_remove(X1, "\\.p[0-9]+$")) %>%
    distinct(transcript_id) %>%
    mutate(domain_present = "yes")

    cat("Transcritos com domínios Pfam significativos (E-value <=", pfam_evalue, "):", nrow(transcripts_with_domain), "\n")
    
    # Combinar todos os dados
    cat("\nCombinando dados...\n")
    gtf_complete <- merge(dup_reference_annot, transcripts_with_domain, by = "transcript_id", all.x = TRUE)
    gtf_complete <- merge(gtf_complete, PLEK[,c("transcript_id", "PLEK")], by = "transcript_id", all.x = TRUE)
    gtf_complete <- merge(gtf_complete, CPC2[,c("transcript_id", "transcript_length", "peptide_length", "CPC2")], by = "transcript_id", all.x = TRUE)
    
    cat("GTF completo criado:", nrow(gtf_complete), "linhas\n")
    
    # Aplicar filtros para seleção de lncRNAs
    cat("\nAplicando filtros para seleção de lncRNAs...\n")
    cat("Critérios:\n")
    cat("  - class_code em: x, i, u\n")
    cat("  - Comprimento >= min_length nt\n")
    cat("  - PLEK == 'Non-coding' E CPC2 == 'noncoding'\n")
    cat("  - Sem domínios Pfam significativos\n")
    cat("\n")
    
    # Contar transcritos em cada etapa do filtro
    step1 <- gtf_complete %>% dplyr::filter(class_code %in% c("x", "i", "u"))
    cat("Após filtro class_code:", length(unique(step1$transcript_id)), "transcritos únicos\n")
    
    step2 <- step1 %>% dplyr::filter(transcript_length >= min_length)
    cat("Após filtro comprimento:", length(unique(step2$transcript_id)), "transcritos únicos\n")
    
    step3 <- step2 %>% dplyr::filter(PLEK == "Non-coding" & CPC2 == "noncoding")
    cat("Após filtros CPC2 e PLEK:", length(unique(step3$transcript_id)), "transcritos únicos\n")
    
    step4 <- step3 %>% dplyr::filter(is.na(domain_present))
    cat("Após filtro domínios Pfam:", length(unique(step4$transcript_id)), "transcritos únicos\n")
    
    # Selecionar lncRNAs finais
    lncRNA <- step4 %>% pull(transcript_id) %>% unique()
    cat("\nlncRNAs selecionados:", length(lncRNA), "transcritos\n")
    
    # Criar GTF de lncRNAs
    gtf_lncRNA <- gtf_complete %>% 
      filter(transcript_id %in% lncRNA) %>%
      mutate(transcript_biotype = ifelse(type == "transcript", "lncRNA", ""))
    
    cat("GTF lncRNA criado:", nrow(gtf_lncRNA), "linhas\n")
    
    # Verificar genes com diferentes class_codes
    cat("\nVerificando genes com diferentes class_codes...\n")
    lncRNAs_check <- gtf_lncRNA %>%
      dplyr::filter(type == "transcript") %>%
      dplyr::select(class_code, gene_id) %>%
      distinct()
    
    duplicated_genes <- nrow(lncRNAs_check) - length(unique(gtf_lncRNA$gene_id))
    cat("Genes com diferentes class_codes:", duplicated_genes, "\n")
    
    # Exportar GTFs intermediários
    cat("\nExportando GTFs...\n")
    rtracklayer::export(gtf_lncRNA, "gtf_lncRNA.gtf")
    rtracklayer::export(gtf_complete, "gtf_complete.gtf")
    cat("  gtf_lncRNA.gtf exportado\n")
    cat("  gtf_complete.gtf exportado\n")
    
    # Combinar com anotação de referência
    cat("\nCombinando com anotação de referência...\n")
    reference <- as.data.frame(rtracklayer::import(reference_gtf))
    reference_filter <- reference %>% filter(type %in% c("transcript", "exon"))
    cat("Referência carregada:", nrow(reference_filter), "linhas (transcript + exon)\n")
    
    combined_stringtie_ncbi <- dplyr::bind_rows(reference_filter, gtf_lncRNA)

    combined_stringtie_ncbi$Dbxref <- sapply(combined_stringtie_ncbi$Dbxref, function(x) {
      if (is.null(x)) return(NA_character_)
      paste(x, collapse = ",")
    })

    combined_stringtie_ncbi$end_range <- sapply(combined_stringtie_ncbi$end_range, function(x) {
      if (is.null(x)) return(NA_character_)
      paste(x, collapse = ",")
    })

    combined_stringtie_ncbi$start_range <- sapply(combined_stringtie_ncbi$start_range, function(x) {
      if (is.null(x)) return(NA_character_)
      paste(x, collapse = ",")
    })

    combined_stringtie_ncbi$Parent <- sapply(combined_stringtie_ncbi$Parent, function(x) {
      if (is.null(x)) return(NA_character_)
      paste(x, collapse = ",")
    })

    combined_stringtie_ncbi$Note <- sapply(combined_stringtie_ncbi$Note, function(x) {
      if (is.null(x)) return(NA_character_)
      paste(x, collapse = ",")
    })

    combined_stringtie_ncbi$gene_id[is.na(combined_stringtie_ncbi$gene_id)] <- "unknown_gene"
    combined_stringtie_ncbi$transcript_id[is.na(combined_stringtie_ncbi$transcript_id)] <- "unknown_transcript"
    combined_stringtie_ncbi$exon_number[is.na(combined_stringtie_ncbi$exon_number)] <- "0"

    # Criar a coluna `attribute` no formato GTF
    combined_stringtie_ncbi$attribute <- paste0(
      'gene_id "', combined_stringtie_ncbi$gene_id, '"; ',
      'transcript_id "', combined_stringtie_ncbi$transcript_id, '"; ',
      'exon_number "', combined_stringtie_ncbi$exon_number, '";'
    )

    # Criar o GRanges com os campos obrigatórios
    combined_stringtie_ncbi <- GRanges(
      seqnames = combined_stringtie_ncbi$seqnames,
      ranges = IRanges(start = combined_stringtie_ncbi$start, end = combined_stringtie_ncbi$end),
      strand = combined_stringtie_ncbi$strand,
      
      # metadados (evite nomes reservados)
      source = combined_stringtie_ncbi$source,
      type = combined_stringtie_ncbi$type,
      class_code = combined_stringtie_ncbi$class_code,
      transcript_id = combined_stringtie_ncbi$transcript_id,
      gene_id = combined_stringtie_ncbi$gene_id,
      db_xref = sapply(combined_stringtie_ncbi$Dbxref, paste, collapse = ","),
      gbkey = combined_stringtie_ncbi$gbkey,
      locus_tag = combined_stringtie_ncbi$locus_tag,
      partial = combined_stringtie_ncbi$partial,
      orig_protein_id = combined_stringtie_ncbi$orig_protein_id,
      orig_transcript_id = combined_stringtie_ncbi$orig_transcript_id,
      product = combined_stringtie_ncbi$product,
      transcript_biotype = combined_stringtie_ncbi$transcript_biotype,
      exon_number = combined_stringtie_ncbi$exon_number,
      pseudo = combined_stringtie_ncbi$pseudo,
      PLEK = combined_stringtie_ncbi$PLEK,
      transcript_length = combined_stringtie_ncbi$transcript_length,
      peptide_length = combined_stringtie_ncbi$peptide_length,
      CPC2 = combined_stringtie_ncbi$CPC2,
      xloc = combined_stringtie_ncbi$xloc,
      tss_id = combined_stringtie_ncbi$tss_id,
      cmp_ref = combined_stringtie_ncbi$cmp_ref,
      transcript_id_stringtie = NA_character_,
      gene_id_stringtie = NA_character_,
      score = combined_stringtie_ncbi$score,
      phase = combined_stringtie_ncbi$phase,
      attribute = combined_stringtie_ncbi$attribute
    )

    cat("GTF combinado criado:", nrow(combined_stringtie_ncbi), "linhas\n")

    # Verificar sobreposição
    overlap_check <- sum(gtf_lncRNA$transcript_id %in% reference$transcript_id)
    cat("Sobreposição transcript_id (esperado 0):", overlap_check, "\n")

    rtracklayer::export(combined_stringtie_ncbi, "combined_stringtie_ncbi.gtf")
    cat("  combined_stringtie_ncbi.gtf exportado\n")

    # Renomear IDs
    cat("\nRenomeando IDs...\n")
    combined_stringtie_ncbi_new_names <- as.data.frame(rtracklayer::import("combined_stringtie_ncbi.gtf"))

    new_names <- combined_stringtie_ncbi_new_names %>%
      dplyr::filter(type == "transcript") %>%
      mutate(transcript_id_new = ifelse(source == "StringTie",
                                        paste0("Pb18_lncRNA", class_code, gsub("MSTRG.", "_", transcript_id)),
                                        transcript_id)) %>%
      mutate(gene_id_new = ifelse(source == "StringTie",
                                  paste0("Pb18_lncRNA", gsub("MSTRG.", "_", gene_id)),
                                  gene_id)) %>%
      dplyr::select(transcript_id, transcript_id_new, gene_id, gene_id_new) %>% distinct()

    cat("Novos nomes criados para:", nrow(new_names), "transcritos\n")

    # Aplicar novos nomes
    combined_stringtie_ncbi_new_names <- left_join(combined_stringtie_ncbi_new_names,
                                                  new_names %>% dplyr::select(transcript_id, transcript_id_new) %>% distinct(), 
                                                  by = "transcript_id")

    combined_stringtie_ncbi_new_names <- left_join(combined_stringtie_ncbi_new_names,
                                                  new_names %>% dplyr::select(gene_id, gene_id_new) %>% distinct(), 
                                                  by = "gene_id")

    # Reorganizar colunas
    combined_stringtie_ncbi_new_names <- combined_stringtie_ncbi_new_names %>%
      dplyr::rename(transcript_id_stringtie = transcript_id) %>%
      dplyr::rename(gene_id_stringtie = gene_id) %>%
      dplyr::rename(transcript_id = transcript_id_new) %>%
      dplyr::rename(gene_id = gene_id_new) %>%
      dplyr::select(seqnames, start, end, width, strand,
                    source, type, score, phase, gene_id,
                    transcript_id, db_xref, gbkey, locus_tag, partial,
                    orig_protein_id, orig_transcript_id, product, transcript_biotype, exon_number,
                    pseudo, PLEK, transcript_length, peptide_length, CPC2,
                    xloc, class_code, tss_id, cmp_ref,
                    transcript_id_stringtie, gene_id_stringtie)

    
    rtracklayer::export(combined_stringtie_ncbi_new_names, "combined_stringtie_ncbi_new_names.gtf")
    cat("  combined_stringtie_ncbi_new_names.gtf exportado\n")
    
    # Estatísticas finais
    cat("\n=== ESTATÍSTICAS FINAIS ===\n")
    cat("Total de lncRNAs identificados:", length(lncRNA), "\n")
    
    # Estatísticas por class_code
    class_stats <- gtf_lncRNA %>% 
      filter(type == "transcript") %>% 
      count(class_code) %>% 
      arrange(desc(n))
    
    cat("Distribuição por class_code:\n")
    for(i in 1:nrow(class_stats)) {
      cat("  ", class_stats$class_code[i], ":", class_stats$n[i], "\n")
    }
    
    # Estatísticas de comprimento
    length_stats <- gtf_lncRNA %>% 
      filter(type == "transcript") %>% 
      summarise(
        mean_length = mean(transcript_length, na.rm = TRUE),
        median_length = median(transcript_length, na.rm = TRUE),
        min_length = min(transcript_length, na.rm = TRUE),
        max_length = max(transcript_length, na.rm = TRUE)
      )
    
    cat("Estatísticas de comprimento:\n")
    cat("  Média:", round(length_stats$mean_length, 2), "nt\n")
    cat("  Mediana:", length_stats$median_length, "nt\n")
    cat("  Mínimo:", length_stats$min_length, "nt\n")
    cat("  Máximo:", length_stats$max_length, "nt\n")
    
    cat("\nProcessamento concluído com sucesso!\n")
    cat("Data/hora de término:", as.character(Sys.time()), "\n")
    sink()
    
    cat('"', Sys.getenv("TASK_PROCESS", "SELECT_LNCRNAS"), '":\n', file="versions.yml", sep="")
    cat('    r-base: "', R.version.string, '"\n', file="versions.yml", append=TRUE, sep="")
    cat('    tidyverse: "', as.character(packageVersion("tidyverse")), '"\n', file = "versions.yml", append = TRUE, sep="")
    cat('    rtracklayer: "', as.character(packageVersion("rtracklayer")), '"\n', file = "versions.yml", append = TRUE, sep="")

    sink()