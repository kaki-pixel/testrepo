# Using IsoformSwitchAnalyze
# Note: Vignette did looking for biggest isoform switches im just looking at the ones with biggest log2fcp


## Set working directory
setwd("C:/Users/kaki/OneDrive - King's College London/Year 3/6BBG0303 Research Project/GTEx analysis")

# Load packages
library(data.table)
library(matrixStats)
library(GenomicFeatures)
library(openxlsx)
library(dplyr)
library(IsoformSwitchAnalyzeR)
library(rtracklayer)
library(ggplot2)

#vjbnsjnvweaa

# Load DTU and Data results, significant transcripts and genes
drimDTU_mito <- read.table("DRIMSeq_Results/HLV_L_MITO_DRIM_results.txt", header = TRUE, sep = "\t")
drimData_mito <- read.table("DRIMSeq_Results/HLV_L_MITO_DRIM_means_and_proportions.txt", header = TRUE, sep = "\t")
mito_sig_t <- read.table("TS_Results/HLV_L_MITO_sig_transcripts.txt", header = TRUE, sep = "\t")
mito_all_t_of_sig_genes <- read.table("TS_Results/HLV_L_MITO_all_transcripts_of_sig_genes.txt", header = TRUE, sep = "\t")

# genes that have significant DTU that has transcript driving it, rank by gene p value
sig_genes_ranked <- drimDTU_mito[drimDTU_mito$transcript < 0.05,
                                 c("GeneID", "GeneSymbol", "gene")]
sig_genes_ranked <- unique(sig_genes_ranked)
sig_genes_ranked <- sig_genes_ranked[order(sig_genes_ranked$gene), ]

# # genes that have significant DTU, IGNORE THIS
# sig_genes_ranked <- drimDTU_mito[drimDTU_mito$gene < 0.05, 
#                                  c("GeneID", "GeneSymbol", "gene")] 
# sig_genes_ranked <- unique(sig_genes_ranked)
# sig_genes_ranked <- sig_genes_ranked[order(sig_genes_ranked$gene), ]

# Load objects
counts_for_isa <- read.table("ISA/HLV_L_counts_for_isa.txt", header = TRUE, sep = "\t")
colnames(counts_for_isa) <- gsub("\\.", "-", colnames(counts_for_isa))
sample_info <- read.table("ISA/HLV_L_sample_info.txt", header = TRUE, sep = "\t")

# Set the isoform_id column as row names
rownames(counts_for_isa) <- counts_for_isa$isoform_id

# # MAKE NEW GTF with only transcript that sample have because kept getting error
# isoform_ids <- rownames(counts_for_isa)   # your isoform IDs from the count matrix
# ori_gtf <- import("gencode.v26.annotation.gtf")
# sub_gtf <- ori_gtf[ori_gtf$transcript_id %in% isoform_ids] # make gtf file that only contains transcripts which samples have
# export(sub_gtf, "v26_filt_annotation.gtf")

# Reformat sample info for design matrix
design_matrix <- data.frame(sampleID = sample_info$sample_id,
                            condition = sample_info$tissue)
design_matrix$condition <- gsub(" - ", "_", design_matrix$condition)
design_matrix$condition <- gsub(" ", "_", design_matrix$condition)

comp_to_make <- data.frame(condition_1 = "Heart_Left_Ventricle",
                           condition_2 = "Lung")

# # Create switchAnalyzeRlist
# aSwitchList <- importRdata(
#   isoformCountMatrix = counts_for_isa,
#   isoformRepExpression = NULL,
#   designMatrix = design_matrix,
#   isoformExonAnnoation = "v26_filt_annotation.gtf",
#   isoformNtFasta = "gencode.v26.transcripts.fa.gz",
#   # removeNonConvensionalChr = TRUE,
#   comparisonsToMake = comp_to_make,
#   addAnnotatedORFs = TRUE,
#   showProgress = TRUE)
# 
#saveRDS(aSwitchList, "ISA/HLV_L_switchAnalyzeRlist.rds")

# LOAD IT IN IF STOPPED HERE
aSwitchList <- readRDS("ISA/HLV_L_switchAnalyzeRlist.rds")

# Check the isoformFeatures name
head(aSwitchList$isoformFeatures)

# Strip version numbers to match with the drimseq
aSwitchList$isoformFeatures$isoform_id_clean <- gsub("\\..*", "", aSwitchList$isoformFeatures$isoform_id)

# Add transcript level q-values from DRIMSeq results
aSwitchList$isoformFeatures$isoform_switch_q_value <- drimDTU_mito$transcript[
  match(aSwitchList$isoformFeatures$isoform_id_clean, drimDTU_mito$TranscriptID)]

# Add gene level q-values from DRIMSeq results
aSwitchList$isoformFeatures$gene_switch_q_value <- drimDTU_mito$gene[
  match(aSwitchList$isoformFeatures$isoform_id_clean, drimDTU_mito$TranscriptID)]

# Get list of top genes (in terms of magnitude of switching)
top_genes <- unique(mito_sig_t[order(mito_sig_t$abs_log2fcp, decreasing = TRUE), "GeneSymbol"])

# Subset the switchlist to keep the significant genes
aSwitchList_top <- subsetSwitchAnalyzeRlist(aSwitchList,
                                            aSwitchList$isoformFeatures$gene_name %in% top_genes)

#saveRDS(aSwitchList_top, "ISA/HLV_L_switchAnalyzeRlist_top.rds")

# Extract sequences, ran this
aSwitchList_top <- extractSequence(
  aSwitchList_top,
  removeLongAAseq = FALSE,
  alsoSplitFastaFile = FALSE,
  pathToOutput = "ISA/extracted_seq_complete",
  writeToFile = FALSE)

# # Extract sequences (for web server), did not run this
# aSwitchList_top <- extractSequence(
#   aSwitchList_top,
#   removeLongAAseq = TRUE,
#   alsoSplitFastaFile = TRUE,
#   pathToOutput = "ISA/extracted_seq_complete",
#   writeToFile = TRUE)

# check number
length(aSwitchList_top$ntSequence) # 301
length(aSwitchList_top$aaSequence) # 167 because some are non-coding

# # Subset switchlist to candidate genes
# aSwitchList_final <- subsetSwitchAnalyzeRlist(
#   aSwitchList_top,
#   aSwitchList_top$isoformFeatures$gene_name %in% final_genes)
# 
# # Extract sequences OF FINAL CANDIDATES for PFAM
# aSwitchList_final <- extractSequence(
#   aSwitchList_final,
#   removeLongAAseq = FALSE,
#   alsoSplitFastaFile = FALSE,
#   pathToOutput = "ISA/pfam_extracted_seq_complete/",
#   writeToFile = TRUE)
# 
# # Import PFAM analysis results
# aSwitchList_final <- analyzePFAM(
#   switchAnalyzeRlist = aSwitchList_final,
#   pathToPFAMresultFile = "ISA/pfam_results.txt",
#   showProgress = TRUE,
#   quiet = FALSE)
# 
# # Analyse functional consequences
# aSwitchList_final <- analyzeSwitchConsequences(
#   aSwitchList_final,
#   consequencesToAnalyze = c("isoform_length", "ORF_length", "NMD_status", "domains_identified"),
#   showProgress = TRUE)

# Look for pairs with AAs in both tissues
# get aa ids
aa_t_ids <- gsub("\\..*", "", names(aSwitchList_top$aaSequence))

interest_enst <- mito_sig_t[, c("GeneSymbol", "TranscriptName", "TranscriptID", 
                                "Class", "dominant_in", "gene", "transcript", 
                                "log2fcp", "abs_log2fcp")]

interest_enst$has_aa <- interest_enst$TranscriptID %in% aa_t_ids

# Keep rows where transcript has aa sequence
interest_enst <- interest_enst[interest_enst$has_aa == TRUE, ]

# Gene pairs
paired_genes <- interest_enst %>%
  group_by(GeneSymbol) %>%
  summarise(n_tissues = n_distinct(dominant_in)) %>%
  filter(n_tissues == 2) # so will deffo have at least 2 transcripts

set_of_pairs_waa <- interest_enst[interest_enst$GeneSymbol %in% paired_genes$GeneSymbol, ] %>%
  arrange(gene, GeneSymbol, dominant_in)

# saveRDS(set_of_pairs_waa, "ISA/set_of_pairs_info.rds")

# Get the specific transcript IDs for each gene pair
final_transcript_ids <- set_of_pairs_waa$TranscriptID
#saveRDS(final_transcript_ids, "ISA/final_transcript_ids.rds")

# Subset switchlist to only those transcripts
aSwitchList_final_pairs <- subsetSwitchAnalyzeRlist(
  aSwitchList_top,
  aSwitchList_top$isoformFeatures$isoform_id_clean %in% final_transcript_ids)

# Get AA sequences for TARGETP, but its fine to analyze PFAM with the previous one

aSwitchList_final_pairs <- extractSequence(
  aSwitchList_final_pairs,
  removeLongAAseq = FALSE,
  alsoSplitFastaFile = FALSE,
  pathToOutput = "ISA/pfam_pair_extracted_seq_complete/",
  writeToFile = TRUE)

# Import PFAM analysis results
aSwitchList_final_pairs <- analyzePFAM(
  switchAnalyzeRlist = aSwitchList_final_pairs,
  pathToPFAMresultFile = "ISA/pfam_results.txt",
  showProgress = TRUE,
  quiet = FALSE)

# Analyse functional consequences (for coding pot need cpc2, intron reten need alt splicing)
aSwitchList_final_pairs <- analyzeSwitchConsequ
ences(
  aSwitchList_final_pairs,
  consequencesToAnalyze = c("isoform_length", "ORF_length", "NMD_status", "domains_identified"),
  showProgress = TRUE)

# SAVE THE FINAL
saveRDS(aSwitchList_final_pairs, "ISA/HLV_L_switchAnalyzeRlist_final_pairs.rds")
aSwitchList_final_pairs <- readRDS("ISA/HLV_L_switchAnalyzeRlist_final_pairs.rds")

switchPlot(aSwitchList_final_pairs, gene = "IDH3B")
switchPlot(aSwitchList_final_pairs, gene = "SLC25A3")
switchPlot(aSwitchList_final_pairs, gene = "ATP5G1")
switchPlot(aSwitchList_final_pairs, gene = "IVD")
switchPlot(aSwitchList_final_pairs, gene = "DCAKD")
switchPlot(aSwitchList_final_pairs, gene = "CCDC51")
switchPlot(aSwitchList_final_pairs, gene = "TSFM")
switchPlot(aSwitchList_final_pairs, gene = "MECR")
switchPlot(aSwitchList_final_pairs, gene = "TMEM143")
switchPlot(aSwitchList_final_pairs, gene = "HEMK1")
switchPlot(aSwitchList_final_pairs, gene = "MTPAP")
switchPlot(aSwitchList_final_pairs, gene = "AMT")

png("MTPAP_switchPlot.png", width = 2200, height = 2000, res = 300)
switchPlot(aSwitchList_final_pairs, gene = "MTPAP")
dev.off()


#switchPlotTopSwitches(aSwitchList_final_pairs)
#extractSwitchSummary(aSwitchList_final_pairs, filterForConsequences = TRUE)
# 
# # global consequence analysis
# extractConsequenceSummary(aSwitchList_final_pairs)
# extractConsequenceEnrichment(aSwitchList_final_pairs)
# extractConsequenceGenomeWide(aSwitchList_final_pairs)

#dev.off()
switchPlot(aSwitchList_final_pairs, gene = "IDH3B")

head(aSwitchList_final_pairs$isoformFeatures)

# MAKE GRAPH FOR RESULTS
library(ggplot2)
library(dplyr)
library(tidyr)

# put dom first
plot_data <- isa_info %>%
  mutate(dominant = ifelse(dIF < 0, "Heart_LV", "Lung")) %>%
  select(gene_name, isoform_id_clean, IF1, IF2, dominant, gene_switch_q_value) %>%
  arrange(gene_switch_q_value, gene_name, dominant) %>%
  mutate(isoform_id_clean = factor(isoform_id_clean, levels = rev(unique(isoform_id_clean))))

# Long format for plotting
plot_long <- plot_data %>%
  pivot_longer(cols = c(IF1, IF2), names_to = "tissue", values_to = "IF") %>%
  mutate(tissue = ifelse(tissue == "IF1", "Heart - Left Ventricle", "Lung"))

p_dumbbell <- ggplot(plot_long, aes(x = IF, y = isoform_id_clean)) +
  geom_line(aes(group = isoform_id_clean), color = "grey75", linewidth = 0.8) +
  geom_point(aes(color = tissue), size = 3.5) +
  scale_color_manual(values = c("Heart - Left Ventricle" = "steelblue", 
                                "Lung" = "darkorange")) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  scale_y_discrete(position = "right") +  # moves transcript IDs to right
  facet_grid(gene_name ~ ., scales = "free_y", space = "free",
             switch = "y") +  # moves gene names to left
  labs(x = "Isoform fraction (IF)", y = NULL, color = NULL,
       title = "Isoform Switching between Heart (Left Ventricle) and Lung") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "top", #put legend on top
        strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0, face = "italic", size = 9),
        axis.text.y = element_text(size = 8),
        panel.spacing = unit(0.3, "lines"),
        plot.title = element_text(hjust = 0.5))

print(p_dumbbell)

ggsave("ISA/dumbbell_isoform_switch.png", p_dumbbell, width = 10, height = 7, dpi = 300)

# # Get FINAL LIST OF CANDIDATE GENES THIS ONE INCLUDED ALL THE TRANSCRIPTS IN THE SWITCHLIST
# final_genes <- unique(set_of_pairs_waa$GeneSymbol)
# 
# # Subset switchlist to candidate genes
# aSwitchList_final <- subsetSwitchAnalyzeRlist(
#   aSwitchList_top,
#   aSwitchList_top$isoformFeatures$gene_name %in% final_genes)
# 
# # Extract sequences OF FINAL CANDIDATES for PFAM
# aSwitchList_final <- extractSequence(
#   aSwitchList_final,
#   removeLongAAseq = FALSE,
#   alsoSplitFastaFile = FALSE,
#   pathToOutput = "ISA/pfam_extracted_seq_complete/",
#   writeToFile = TRUE)
# 
# # Import PFAM analysis results
# aSwitchList_final <- analyzePFAM(
#   switchAnalyzeRlist = aSwitchList_final,
#   pathToPFAMresultFile = "ISA/pfam_results.txt",
#   showProgress = TRUE,
#   quiet = FALSE)
# 
# # Analyse functional consequences
# aSwitchList_final <- analyzeSwitchConsequences(
#   aSwitchList_final,
#   consequencesToAnalyze = c("isoform_length", "ORF_length", "NMD_status", "domains_identified"),
#   showProgress = TRUE)
# 
# final_genes
# 
# switchPlot(aSwitchList_final, gene = "IDH3B")
# switchPlot(aSwitchList_final, gene = "SLC25A3")
# switchPlot(aSwitchList_final, gene = "ATP5G1")
# switchPlot(aSwitchList_final, gene = "IVD")
# switchPlot(aSwitchList_final, gene = "DCAKD")
# switchPlot(aSwitchList_final, gene = "CCDC51")
# switchPlot(aSwitchList_final, gene = "TSFM")
# switchPlot(aSwitchList_final, gene = "MECR")
# switchPlot(aSwitchList_final, gene = "TMEM143")
# switchPlot(aSwitchList_final, gene = "HEMK1")
# switchPlot(aSwitchList_final, gene = "MTPAP")
# switchPlot(aSwitchList_final, gene = "AMT")
# switchPlotTopSwitches(aSwitchList_final)

?analyzePFAM
?analyzeSwitchConsequences

# # check the unchanged usage example
# mito_all[mito_all$TranscriptID == "ENST00000395338", 
#          c("GeneSymbol", "TranscriptName", "TranscriptID", "transcript", "log2fcp")]

# # paired_genes <- names(table(interest_enst$GeneSymbol)[table(interest_enst$GeneSymbol) >= 2])
# # set_of_pairs_waa <- interest_enst[interest_enst$GeneSymbol %in% paired_genes, ] %>%
# #arrange(GeneSymbol, dominant_in)
# # get all isoform IDs in switchlist
# all_isoforms <- aSwitchList_top_60$isoformFeatures$isoform_id
# 
# # get isoforms that have AA sequences
# has_aa <- names(aSwitchList_top_60$aaSequence)
# 
# # find ones WITHOUT aa
# missing_aa <- all_isoforms[!all_isoforms %in% has_aa]
# head(missing_aa) # try ENST00000460872.1
# 
# # check if this transcript has ORF info in your switchlist
# aSwitchList_top_60$orfAnalysis[
#   aSwitchList_top_60$orfAnalysis$isoform_id == "ENST00000460872.1",
#   c("isoform_id", "orfTransciptStart", "orfTransciptEnd", "orfTransciptLength")]
# 
# # Check if the transcript exists in orfAnalysis at all
# "ENST00000460872.1" %in% aSwitchList_top_60$orfAnalysis$isoform_id
# 
# # Check ORF analysis for both transcripts
# aSwitchList_top_60$orfAnalysis[
#   aSwitchList_top_60$orfAnalysis$isoform_id %in% c("ENST00000341695.9", "ENST00000348706.9"),
# ]
# 
# # Check the nucleotide sequences slot
# aSwitchList_top_60$ntSequence[
#   names(aSwitchList_top_60$ntSequence) %in% c("ENST00000341695.9", "ENST00000348706.9")
# ]
# 
# # Check if sequences were in the original switchlist before any filtering THESE ARE THE PROTEIN CODING ONES BUT NO SEQ
# aSwitchList_top_60$isoformFeatures[
#   aSwitchList_top_60$isoformFeatures$isoform_id %in% c("ENST00000341695.9", "ENST00000348706.9"),
#   c("isoform_id", "gene_name")
# ]
# aSwitchList_top_60$ntSequence["ENST00000341695.9"]
# "ENST00000341695.9" %in% names(aSwitchList_top_60$ntSequence)
# 
# # Or check the aa sequences slot
# aSwitchList_top_60$aaSequence[
#   names(aSwitchList_top_60$aaSequence) %in% c("ENST00000341695.9", "ENST00000348706.9")
# ]
# length(aSwitchList_top_60$ntSequence)
# head(names(aSwitchList_top_60$ntSequence))
# 
# length(aSwitchList_top_60$ntSequence)
# nrow(aSwitchList_top_60$orfAnalysis)
# 
# aa_ids <- gsub("\\..*", "", names(aSwitchList_top_60$aaSequence))
# interest_enst2 <- mito_sig_t[mito_sig_t$GeneSymbol %in% top_60_genes, 
#                              c("GeneSymbol", "TranscriptName", "TranscriptID", "Class", "dominant_in")]
# interest_enst2$has_aa <- interest_enst2$TranscriptID %in% aa_ids
# interest_enst2
# 
# # keep only transcripts with AA sequences
# complete <- interest_enst2[interest_enst2$has_aa == TRUE, ]
# 
# # find genes that have BOTH heart and lung transcripts
# paired_genes <- names(table(complete$GeneSymbol)[table(complete$GeneSymbol) >= 2])
# 
# # show only those genes
# set_of_pairs_waa <- complete[complete$GeneSymbol %in% paired_genes, ] %>%
#   arrange(GeneSymbol, dominant_in)
# 
# switching_pairs <- c(set_of_pairs_waa$TranscriptID)
# 
# aSwitchList_top_60$isoformFeatures %>%
#   mutate(isoform_clean = gsub("\\..*", "", isoform_id)) %>%
#   filter(isoform_clean %in% switching_pairs) %>%
#   select(gene_name, isoform_id, IF1, IF2, condition_1, condition_2) %>%
#   arrange(gene_name)
# 
# switchPlot(aSwitchList_extended, gene = "SLC25A3")
# switchPlot(aSwitchList_extended, gene = "ATP5G1")
# switchPlotTopSwitches(aSwitchList_top_60)
# 
# 
# aSwitchList_top_60 <- analyzeSwitchConsequences(
#   aSwitchList_top_60,
#   consequencesToAnalyze = c("isoform_length", "ORF_length", "NMD_status", "domains_identified"),
#   showProgress = TRUE)
# 
# # Run CPC2 and PFAM
# # Downloaded results
# # IMPORT IT BACK IN 
# aSwitchList_top_6_analysed <- analyzeCPC2(
#   switchAnalyzeRlist = aSwitchList_top_6,
#   pathToCPC2resultFile = "isa_extracted_seq/result_cpc2.txt",
#   removeNoncodinORFs = TRUE)

