################################################################################
###################### Select lncRNAs (back to shell) ##########################
################################################################################
library(tidyverse)
library(dplyr)
library(stringr)

    annotated_gtf <- Sys.getenv("ANNOTATED_GTF")
    cpc2_results <- Sys.getenv("CPC2_RESULTS")
    plek_results <- Sys.getenv("PLEK_RESULTS")
    pfam_domtblout <- Sys.getenv("PFAM_DOMTBLOUT")
    reference_gtf <- Sys.getenv("REFERENCE_GTF")
    prefix <- Sys.getenv("PREFIX")
    min_length <- as.numeric(Sys.getenv("MIN_LENGTH"))
    pfam_evalue <- as.numeric(Sys.getenv("PFAM_EVALUE"))

dup_reference_annot <- as.data.frame(rtracklayer::import(annotated_gtf))

# dup_reference_annot <- dup_reference_annot %>%
#   mutate(
#     gene_id = ifelse(
#       str_detect(gene_id, "gene-PADG_"),
#       str_extract(gene_id, "gene-PADG_[0-9]+"),
#       gene_id  # mantém MSTRG.xxx se não tiver gene-PADG
#     )
#   )

CPC2 <- read.delim(cpc2_results)
colnames(CPC2) <- c("transcript_id", "transcript_length", "peptide_length", "Fickett_score",
                    "pI", "ORF_integrity", "coding_probability", "CPC2")

PLEK <- read.delim(plek_results, header = F)
PLEK <- PLEK %>% mutate(transcript_id = str_extract(V3, "(?<=\\>).*?(?= loc:)"))
colnames(PLEK) <- c("PLEK", "value", "features", "transcript_id")

pfam <- readr::read_table(pfam_domtblout, comment = "#", col_names = FALSE)



# colnames(pfam) <- c(
#  "target_name", "target_accession", "tlen",
#  "query_name", "query_accession", "qlen",
#  "full_seq_Evalue", "full_seq_score", "full_seq_bias",
#  "domain_num", "domain_total",
#  "domain_c_Evalue", "domain_i_Evalue", "domain_score", "domain_bias",
#  "hmm_from", "hmm_to", "ali_from", "ali_to", "env_from", "env_to",
#  "acc", "description"
# )


transcripts_with_domain <- as.data.frame(pfam) %>%
  dplyr::filter(X7 <= pfam_evalue) %>%
  mutate(transcript_id = str_remove(X1, "\\.p[0-9]+$")) %>% # Remove .p1, .p2 etc
  distinct(transcript_id) %>% 
  mutate(domain_present = "yes")

gtf_complete <- merge(dup_reference_annot, transcripts_with_domain, by = "transcript_id", all.x = T)
gtf_complete <- merge(gtf_complete, PLEK[,c("transcript_id", "PLEK")], by = "transcript_id", all.x = T)
gtf_complete <- merge(gtf_complete, CPC2[,c("transcript_id", "transcript_length", "peptide_length", "CPC2")], by = "transcript_id", all.x = T)

rm(pfam, CPC2, PLEK, transcripts_with_domain)

lncRNA <- gtf_complete %>% 
  dplyr::filter(class_code %in% c("x", "i", "u")) %>% 
  dplyr::filter(transcript_length >= min_length) %>% 
  dplyr::filter(CPC2 == "noncoding" & PLEK == "Non-coding") %>% 
  dplyr::filter(is.na(domain_present)) %>% pull(transcript_id) %>% unique()

gtf_lncRNA <- gtf_complete %>% filter(transcript_id %in% lncRNA) %>% 
  mutate(transcript_biotype = ifelse(type == "transcript", "lncRNA", ""))


rtracklayer::export(gtf_lncRNA, "gtf_lncRNA.gtf", format = "gtf")
rtracklayer::export(gtf_complete, "gtf_complete.gtf", format = "gtf")

rm(dup_reference_annot, gtf_complete, lncRNA)

################################################################################
#################### Combine reference gtf and stringtie #######################
################################################################################

## Check gene equal with different class
lncRNAs_check <- gtf_lncRNA %>% 
  dplyr::filter(type == "transcript") %>% 
  dplyr::select(class_code, gene_id) %>% 
  distinct()

## 15 genes with different class
nrow(lncRNAs_check) - length(unique(gtf_lncRNA$gene_id))

## Select genes
cut_lncRNA <- lncRNAs_check %>% dplyr::filter(duplicated(lncRNAs_check$gene_id)) %>% pull(gene_id)

#gtf_lncRNA %>% dplyr::filter(gene_id %in% cut_lncRNA) %>% 
#  dplyr::filter(type == "transcript") %>% 
#  dplyr::select(start, end, strand, gene_id, transcript_id, class_code, xloc, transcript_length)



rm(lncRNAs_check, cut_lncRNA)

## Filter transcript and exon from reference
reference <- as.data.frame(rtracklayer::import(reference_gtf))
reference_filter <- reference %>% filter(type %in% c("transcript", "exon"))

#gtf_lncRNA <- GenomicRanges::makeGRangesFromDataFrame(gtf_lncRNA, keep.extra.columns = TRUE)

# reference_filter <- makeGRangesFromDataFrame(
#   reference_filter,
#   seqnames.field = "seqnames",
#   start.field = "start",
#   end.field = "end",
#   strand.field = "strand",
#   keep.extra.columns = TRUE
# )

combined_stringtie_ncbi <- dplyr::bind_rows(reference_filter, gtf_lncRNA)
# combined_stringtie_ncbi <- c(reference_filter, gtf_lncRNA)



sum(gtf_lncRNA$transcript_id %in% reference$transcript_id) #0 is expected

# Exportar como novo GTF
rtracklayer::export(combined_stringtie_ncbi, "combined_stringtie_ncbi.gtf", format = "gtf")

#rm(combined_stringtie_ncbi, gtf_lncRNA, reference, reference_filter)

################################################################################
################################# Change name ##################################
################################################################################

combined_stringtie_ncbi_new_names <- as.data.frame(combined_stringtie_ncbi)

new_names <- combined_stringtie_ncbi_new_names %>% 
  dplyr::filter(type == "transcript") %>% 
  mutate(transcript_id_new = ifelse(source == "StringTie",
                                    paste0("Pb18_lncRNA", class_code, gsub("MSTRG.", "_", transcript_id)), 
                                    transcript_id)) %>% 
  mutate(gene_id_new = ifelse(source == "StringTie",
                              paste0("Pb18_lncRNA", gsub("MSTRG.", "_", gene_id)), 
                              gene_id)) %>% 
  dplyr::select(transcript_id, transcript_id_new, gene_id, gene_id_new) %>% distinct()


combined_stringtie_ncbi_new_names <- left_join(combined_stringtie_ncbi_new_names, 
                                               new_names %>% dplyr::select(transcript_id, transcript_id_new) %>% distinct(), by = "transcript_id")

combined_stringtie_ncbi_new_names <- left_join(combined_stringtie_ncbi_new_names, 
                                               new_names %>% dplyr::select(gene_id, gene_id_new) %>% distinct(), by = "gene_id")

combined_stringtie_ncbi_new_names <- combined_stringtie_ncbi_new_names %>% 
  dplyr::rename(transcript_id_stringtie = transcript_id) %>% 
  dplyr::rename(gene_id_stringtie = gene_id) %>% 
  dplyr::rename(transcript_id = transcript_id_new) %>%
  dplyr::rename(gene_id = gene_id_new) %>% 
  dplyr::select(seqnames,start, end,width,strand,
                source, type,score,phase,gene_id,
                transcript_id,db_xref,gbkey,locus_tag,partial,
                orig_protein_id,orig_transcript_id,product, transcript_biotype,exon_number,
                pseudo,PLEK,transcript_length,peptide_length,CPC2,
                xloc,class_code,tss_id,gene_name,cmp_ref,
                transcript_id_stringtie,gene_id_stringtie)

combined_stringtie_ncbi_new_names %>% dplyr::select(gene_id,gene_id_stringtie,
                                                    transcript_id, transcript_id_stringtie, class_code)

rtracklayer::export(combined_stringtie_ncbi_new_names, "combined_stringtie_ncbi_new_names.gtf", format = "gtf")

    cat('"', Sys.getenv("TASK_PROCESS", "SELECT_LNCRNAS"), '":\n', file="versions.yml", sep="")
    cat('    r-base: "', R.version.string, '"\n', file="versions.yml", append=TRUE, sep="")
    cat('    tidyverse: "', as.character(packageVersion("tidyverse")), '"\n', file = "versions.yml", append = TRUE, sep="")
    cat('    rtracklayer: "', as.character(packageVersion("rtracklayer")), '"\n', file = "versions.yml", append = TRUE, sep="")
    