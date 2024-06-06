```{r}

library(apeglm)
library(ape)
library(Biostrings)
library(DESeq2)
library(circlize)
library(clusterProfiler)
library(ComplexHeatmap)
library(CorLevelPlot)
library(corrplot)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggtree)
library(ggVennDiagram)
library(GO.db)
library(gridExtra)
library(openxlsx)
library(phangorn)
library(pheatmap)
library(RColorBrewer)
library(tidyr)
library(tidyverse)
library(treeio)
library(WGCNA)


#Upload data, filtering and set data

matriz_cruda <- read.delim("Route_to_your_file", row.names = 1)
counts <- matriz_cruda[, c("C1", "C2", "C3", "C4", "C5", "T1", "T2", "T3", "T4", "T5")]
mcounts=counts[1e5>apply(counts,1,max),]
mcounts=mcounts[0<apply(mcounts[,1:10],1,var),]

conditions <- factor(c(rep("Control", 5), rep("Treatment", 5)))
metadata <- data.frame(condition = conditions)
names <- c("C1", "C2", "C3", "C4", "C5", "T1", "T2", "T3", "T4", "T5")
rownames(metadata) <- names
all(rownames(metadata) == colnames(counts))


#Differential expression analysis

dds <- DESeqDataSetFromMatrix(countData = mcounts, colData = metadata, design = ~ condition)
dds <- estimateSizeFactors(dds)
dds
sizeFactors(dds)

head(counts(dds,normalized=TRUE))
boxplot(counts(dds,normalized=TRUE),las=2)

dds <- DESeq(dds, test="Wald")
res <- results(dds)
resultsNames(dds)
plotDispEsts(dds)
plotMA(res, ylim=c(-6,10))
abline(h=c(-2,2), col="#FF2400", lwd=2)

resLFC <- lfcShrink(dds, coef="condition_Treatment_vs_Control", type="apeglm")

plotMA(resLFC, ylim=c(-6,10))
abline(h=c(-2,2), col="#FF2400", lwd=2)

Padj=data.frame(p=res$padj[resLFC$padj<0.95])
ggplot(Padj, aes(x=p)) + geom_histogram(breaks = seq(0,0.95,0.05)) +
  labs(y='frecuencia') + scale_x_continuous(breaks = seq(0,0.95,0.05))

#Rlog normalizatión & corrplot

rld <- rlog(dds, blind=FALSE)

boxplot(assay(rld), pch=20, las=2)

corr <- corrplot(cor(assay(rld)), cl.lim = c(0.95,1), is.corr = F, tl.col = "black",
                 col = colorRampPalette(c("#0099FF","ivory","#FF0090"))(200))

resOrdered <- res[order(res$padj),]
resOrdered

resSig <- subset(res, padj < 0.05)


#Quality Control

sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(counts)
colnames(sampleDistMatrix) <- colnames(counts)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors, fontsize = 16)

pcaData <- plotPCA(rld, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pca_plot <- ggplot(pcaData, aes(PC1, PC2, color = condition, label = rownames(pcaData))) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed()
pca_plot <- pca_plot + scale_color_manual(values = c("#0099FF", "#FF0090"))
pca_plot <- pca_plot + 
  stat_ellipse(aes(group = condition, color = condition), type = "norm", level = 0.95,
               linetype = "dashed", fill = NA)
pca_plot <- pca_plot + geom_text_repel(
  aes(label = colnames(counts)),
  nudge_x = 0.1,
  nudge_y = 0.1,
  segment.size = 0.2,
  box.padding = 0.5,
  max.overlaps = Inf) 
pca_plot

SigGenesLvV <- as.data.frame(res) %>%
  rownames_to_column("GeneID")  %>%
  rename(logFC=log2FoldChange, FDR=padj)
filtTab <- SigGenesLvV %>% 
  filter(!is.na(FDR)) %>% 
  mutate(`-log10(FDR)` = -log10(FDR))
custom_colors <- c("Positive" = "#C21807", "Negative" = "#2823bc", "Non-significant" = "grey")
ggplot(filtTab, aes(x = logFC, y = `-log10(FDR)`)) + 
  geom_point(aes(colour = factor(ifelse(FDR > 0.05 | abs(logFC) < 1, "Non-significant", ifelse(logFC > 0, "Positive", "Negative")))), size = 2) +
  scale_x_continuous(limits = c(-4, 6)) +
  scale_y_continuous(limits = c(0, 20)) +
  scale_colour_manual(values = custom_colors) +
  labs(colour = "FDR Category") +
  theme(legend.position = "top") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "darkred", linewidth = 1) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "darkblue", linewidth = 1)


#Export DEG ti anotate and import result

DEGs = "DEGs.tsv"

write.table(resSig, file = DEGs, sep = "\t", quote = FALSE)

degs <- read.delim("Route_to_your_file", row.names = NULL)
degs <- na.omit(degs)

no_filt_degs <-  subset(no_filt_degs, abs(log2FoldChange) > 1)

names(degs)[names(degs) == "row.names"] <- "gene_id"
names(degs)[names(degs) == "X.baseMean"] <- "baseMean"

#Sequence extraction

fastap <- "Route_to_your_file"
protein_sequences <- readAAStringSet(fastap)
pids <- degs$protein_id
sequence_names <- names(protein_sequences)
extracted_pids <- gsub("^([^ ]+) .*", "\\1", sequence_names)
matching_indices <- match(pids, extracted_pids)
selected_proteins <- protein_sequences[matching_indices]

Interpro_degs <- "Interpro_degs.fasta"
writeXStringSet(selected_proteins, file = Interpro_degs)

#Gene Ontology enrichment analysis and GSEA

gaf <- read.delim("Route_to_your_file", header = TRUE, comment.char = "!", stringsAsFactors = FALSE)

del <- "GeneID:"
degs$GeneID <- gsub(del, "", degs$GeneID)
saved_degs <- degs
degs <- subset(degs, abs(log2FoldChange) > 1)

annotations <- gaf[, c("GeneID", "GO_ID", "Aspect", "With.From")]
ids_to_gaf <- degs$GeneID
annotated_degs <- annotations[annotations$GeneID %in% ids_to_gaf, ]

number_of_degs <- "GeneID"
unique_degs <- unique(annotated_degs[[number_of_degs]])
length(unique_degs)

unique_go_ids <- na.omit(unique(annotated_degs$GO_ID))
go_annotation <- Term(GOTERM[unique_go_ids])
terms_to_id <- as.data.frame(go_annotation)
terms_to_id <- terms_to_id %>% 
  rownames_to_column(var = "RowNames")
colnames(terms_to_id) <- c("GO_ID", "GO_TERM")

ids_GO <- unique(gaf$GO_ID)
ids_GO_annotation <- Term(GOTERM[ids_GO])
terms_to_GO <- as.data.frame(ids_GO_annotation)
terms_to_GO <- terms_to_GO %>% 
  rownames_to_column(var = "RowNames")
colnames(terms_to_GO) <- c("GO_ID", "GO_TERM")

term2gene=gaf[, c("GO_ID", "GeneID")]
term2name=terms_to_GO[, c("GO_ID", "GO_TERM")]


cluster_degs = degs[, "log2FoldChange"]
names(cluster_degs) = as.character(degs[, "GeneID"])

genes <- names(cluster_degs)[abs(cluster_degs) > 0.5]
sorted_genes = sort(cluster_degs, decreasing = TRUE)

xx = enricher(genes, TERM2GENE = term2gene, TERM2NAME = term2name)

yy = GSEA(sorted_genes, TERM2GENE = term2gene, TERM2NAME = term2name)

xx_df <- as.data.frame(xx)
yy_df <- as.data.frame(yy)

a <- enrichKEGG(gene = genes, 
                organism = "hsy",
                keyType = "kegg",
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH",
                qvalueCutoff = 0.2,
                use_internal_data = FALSE)
a_df <- as.data.frame(a)

b <- gseKEGG(geneList = sorted_genes,
             organism = 'hsy',
             pvalueCutoff = 0.05,
             verbose  = FALSE,
             use_internal_data = FALSE)
b_df <- as.data.frame(b)

c <- enrichMKEGG(gene = genes,
                 organism = 'hsy',
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.2)
c_df <- as.data.frame(c)

#Dotplots

xx_dot <- dotplot(xx, showCategory = 13) + ggtitle("Enrichment Analysis GO")
xx_dot

yy_dot <- dotplot(yy, showCategory = 17) + ggtitle("GSEA GO")
yy_dot

a_dot <- dotplot(a, showCategory = 20) + ggtitle("Enrichment Analysis KEGG")
a_dot

b_dot <- dotplot(b, showCategory = 10) + ggtitle("GSEA  KEGG")
b_dot

#Interesting Domains

InterPro_IDs <- read.delim("Route_to_your_file", row.names = NULL, header = FALSE)
Interpro_colnames <- c("Protein_id", "MD5_digest", "lenght", "db", "db_id", "description", "start_location", "stop_location", "score", "status", "date", "Interpro_IDs", "domain_name", "Interpro_GO_IDs", "GO_name")
colnames(InterPro_IDs) <- Interpro_colnames

domains <- InterPro_IDs[, c("Protein_id", "db", "db_id", "description", "Interpro_IDs", "domain_name")]

id2seq <- degs[, c("gene_id", "GeneID", "protein_id")]
go2dom <- yy_df
rownames(go2dom) <- NULL
ko2dom <- xx_df
rownames(ko2dom) <- NULL

rows_EA_GO <- c(4, 7, 12, 13)
rows_GSEA_GO <- c(3, 4, 6, 8, 9, 10, 11, 13, 16)

GSEA_GO_rows <- go2dom[rows_GSEA_GO, "core_enrichment"]
combined_values_GSEA <- paste(GSEA_GO_rows, collapse = "/")
split_values_GSEA <- unlist(strsplit(combined_values_GSEA, "/"))
go2domains_GSEA <- unique(split_values_GSEA)

EA_GO_rows <- ko2dom[rows_EA_GO, "geneID"]
combined_values_EA <- paste(EA_GO_rows, collapse = "/")
split_values_EA <- unlist(strsplit(combined_values_EA, "/"))
go2domains_EA <- unique(split_values_EA)

go2domains_1 <- union(go2domains_GSEA, go2domains_EA)

id2prot <- id2seq[id2seq$GeneID %in% go2domains_1, ]

interesting_domains <- domains[domains$Protein_id %in% id2prot$protein_id, ]

annotated_degs <- annotated_degs %>% 
  mutate(GeneID = as.character(GeneID))

unannotated_degs <- anti_join(degs, annotated_degs, by = "GeneID")

unnanotated_domains <- domains[domains$Protein_id %in% unannotated_degs$protein_id, ]

process_terms <- function(terms_names, terms, table, columns) {
  results <- list()
  
  for (j in seq_along(columns)) {
    column <- columns[j]
    
    for (i in seq_along(terms)) {
      term <- terms[i]
      term_name <- terms_names[i]
      
      
      count <- length(grep(term, table[[column]], ignore.case = TRUE))
      
      
      rows <- table[grepl(term, table[[column]], ignore.case = TRUE), ]
      
     
      unique_vals <- unique(rows$Protein_id)
      
      
      results[[paste0(term_name, "_count_", column)]] <- count
      results[[paste0(term_name, "_rows_", column)]] <- rows
      results[[paste0("unique_", term_name, "_", column)]] <- unique_vals
    }
  }
  
  return(results)
}

terms <- c("A2M", "ABD", "Actinoporin", "Anemone", "Antibacterial", "Antibiotic", "Antimicrobial", "Ankyrin", "AP-1", "Arabinose", "ATF", "B2M", "Bacterial", "Bactericidal", "BGRP", "bZIP", "C-type", "CAP ", "CARD", "Caspase", "CCL", "Chemokine", "Clathrin", "CLR", "Cnid", "Coagulation", "Complement", "Coral", "CTL", "CUB", "CXCR", "Cysteine.*Rich", "Cytolysin", "Cytotoxic", "DEAD", "Death", "Defensin", "EGF", "ERK", "Fibrinogen", "Fibronectin", "Flagellin", "FN3", "FNIII", "FOS", "Fucos", "Galactos", "Glucan", "Glycoside.*Hydrolase", "Glycosyl", "GPCR", "Granzyme", "Heat.*Shock", "Histamine", "HSP", "Ig-", "Immune", "Immunity", "Immunoglobulin", "Integrin", "Interferon", "INF", "Interleukin", "IPT ", "IPS", "ITAM", "ITIM", "JAK", "JNK", "JUN", "Kazal", "LBP", "LDL", "Lectin", "Leucine.*Rich", "Lipopolysaccharide", "Lipoprotein", "Litaf", "LRR", "Lyzosyme", "Mannose", "MAP K", "MAPK", "Mitogen", "Muramic", "Myeloid", "NACHT", "NF.*Kappa", "NF-Y", "NKAP", "NLR", "NOD", "Pathogen", "Pattern", "Permeability", "Peptidoglycan", "Perforin", "PGB", "PGRP", "Phospholipase", "Plant", "Pore", "PPG", "Protease", "Protein .*tyrosine", "PRR", "PYD", "RBL ", "Retino", "Rhamnos", "Receptor", "Receptor-like", "Scavenger", "SRC", "SEFIR", "SH2", "SH3", "ShKT", "SPINK", "SRC", "SRCR", "Sushi", "Tetramine", "TGF", "Threalose", "TIR", "TNF", "Tox-", "Toxin", "Transactivator", "Trypsin", "Tumor", "Tyrosin.*ase", "TX", "VWA", "WD", "Willebrand", "Wnt")

terms_names <- c("A2M", "ABD", "Actinoporin", "Anemone", "Antibacterial", "Antibiotic", "Antimicrobial", "Ankyrin", "AP_1", "Arabinose", "ATF", "B2M", "Bacterial", "Bactericidal", "BGRP", "bZIP", "C_type_lectins", "CAP", "CARD", "Caspase", "CCL", "Chemokine", "Clathrin", "CLR", "Cnid", "Coagulation", "Complement", "Coral", "CTL", "CUB", "CXCR", "Cysteine_Rich", "Cytolysin", "Cytotoxic", "DEAD", "Death", "Defensin", "EGF", "ERK", "Fibrinogen", "Fibronectin", "Flagellin", "FN3", "FNIII", "FOS", "Fucos", "Galactose", "Glucan", "Glycoside_Hydrolase", "Glycosyl", "GPCR", "Granzyme", "Heat_Shock", "Histamine", "HSP", "Ig", "Immune", "Immunity", "Immunoglobulin", "Integrin", "Interferon", "INF", "Interleukin", "IPT", "IPS", "ITAM", "ITIM", "JAK", "JNK", "JUN", "Kazal", "LBP", "LDL", "Lectin", "Leucine_Rich", "Lipopolysaccharide", "Lipoprotein", "Litaf", "LRR", "Lyzosyme", "Mannose", "MAP_K", "MAPK", "Mitogen", "Muramic", "Myeloid", "NACHT", "NF_Kappa", "NF_Y", "NKAP", "NLR", "NOD", "Pathogen", "Pattern", "Permeability", "Peptidoglycan", "Perforin", "PGB", "PGRP", "Phospholipase", "Plant", "Pore", "PPG", "Protease", "Protein_Tyrosine", "PRR", "PYD", "RBL", "Retino", "Rhamnose", "Receptor", "Receptor_like", "Scavenger", "SRC", "SEFIR", "SH2", "SH3", "ShKT", "SPINK", "SRC", "SRCR", "Sushi", "Tetramine", "TGF", "Threalose", "TIR", "TNF", "Tox", "Toxin", "Transactivator", "Trypsin", "Tumor", "Tyrosin_ase", "TX", "VWA", "WD", "Willebrand", "Wnt")

columns <- c("description", "domain_name")

results_unnanot <- process_terms(terms_names, terms, unnanotated_domains, columns)
results_anot <- process_terms(terms_names, terms, interesting_domains, columns)

matching_indices_anot <- grep(common_part, names(results_anot))
matching_objects_anot <- results_anot[matching_indices_anot]
all_terms_annot <- unlist(matching_objects_anot)
um_unique_terms_annot <- length(unique(all_terms_annot))

matching_indices_unnanot <- grep(common_part, names(results_unnanot))
matching_objects_unnanot <- results_unnanot[matching_indices_unnanot]
all_terms_unnannot <- unlist(matching_objects_unnanot)
um_unique_terms_unnannot <- length(unique(all_terms_unnannot))

result_df_annot <- data.frame()


for (i in seq_along(matching_objects_anot)) {
  name <- names(matching_objects_anot)[i]
  num_terms <- length(matching_objects_anot[[i]])
  row <- data.frame(name = name, num_terms = num_terms)
  result_df_annot <- rbind(result_df_annot, row)
}

result_df_unnannot <- data.frame()

for (i in seq_along(matching_objects_unnanot)) {
  name <- names(matching_objects_unnanot)[i]
  num_terms <- length(matching_objects_unnanot[[i]])
  row <- data.frame(name = name, num_terms = num_terms)
  result_df_unnannot <- rbind(result_df_unnannot, row)
}

result_df_annot$name <- gsub("unique_", "", result_df_annot$name)
result_df_annot$name <- gsub("_description", "", result_df_annot$name)
result_df_annot$name <- gsub("_domain_name", "", result_df_annot$name)
sum_df_annot <- aggregate(num_terms ~ name, data = result_df_annot, sum)

result_df_unnannot$name <- gsub("unique_", "", result_df_unnannot$name)
result_df_unnannot$name <- gsub("_description", "", result_df_unnannot$name)
result_df_unnannot$name <- gsub("_domain_name", "", result_df_unnannot$name)
sum_df_unnannot <- aggregate(num_terms ~ name, data = result_df_unnannot, sum)

comparisson <- merge(sum_df_annot, sum_df_unnannot, by = "name", suffixes = c("_Annotated Genes", "_Unnanotated Genes"))

names(comparisson)[names(comparisson) == "num_terms_Annotated Genes"] <- "Annotated Genes"
names(comparisson)[names(comparisson) == "num_terms_Unnanotated Genes"] <- "Unnannotated Genes"

plot_domains <- comparisson %>%
  pivot_longer(cols = ends_with("Genes"), names_to = "Annotation", values_to = "Number of IDCP")

plot_domains_filtered <- plot_domains %>%
  filter(`Number of IDCP` != 0)

p <- ggplot(plot_domains_filtered, aes(x = name, y = `Number of IDCP`, fill = Annotation)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Comparison of Gene Annotations", x = "Domain", y = "Number of IDCP") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))
p


only_annot <- subset(sum_df_annot, num_terms != 0)
only_unnannot <- subset(sum_df_unnannot, num_terms != 0)

venn_annot <- only_annot$name
venn_unnannot <- only_unnannot$name

venn_lists <- list(Annotated = venn_annot, Unnannotated = venn_unnannot)
venn <- Venn(venn_lists)
venn_data <- process_data(venn)
venn_plot <- ggplot() +
  geom_polygon(aes(Y, X, fill = name, group = id), 
               data = venn_regionedge(venn_data),
               show.legend = FALSE) +
  geom_path(aes(Y, X, color = id, group = id), 
            data = venn_setedge(venn_data), 
            linewidth = 1,
            show.legend = FALSE) +
  geom_text(aes(Y, X, label = name),
            size = 7,
            data = venn_setlabel(venn_data)) +
  geom_label(aes(Y, X, label = count),
             size = 8,
             data = venn_regionlabel(venn_data),
             alpha = 0) +
  scale_fill_manual(values = c("Annotated" = "#F8766D", "Unnannotated" = "#00BFC4", "Annotated/Unnannotated" = "#7c9a98")) + 
  scale_color_manual(values = c("1"= "black", "2" = "black")) +
  coord_equal() +
  theme_void() 

venn_plot


annot_count <- Filter(function(x) any(x != 0), updated_matching_objects_annot)
unnanot_count <- Filter(function(x) any(x != 0), updated_matching_objects_unnannot)

text_for_list1 <- "_description"
text_for_list2 <- "_domain_name"

annot_description_i <- grep(text_for_list1, names(matching_objects_anot))
annot_domain_name_i <- grep(text_for_list2, names(matching_objects_anot))

unnannot_description_i <- grep(text_for_list1, names(matching_objects_unnanot))
unnannot_domain_name_i <- grep(text_for_list2, names(matching_objects_unnanot))

annot_description <- matching_objects_anot[annot_description_i]
annot_domain_name <- matching_objects_anot[annot_domain_name_i]

unnannot_description <- matching_objects_unnanot[unnannot_description_i]
unnannot_domain_name <- matching_objects_unnanot[unnannot_domain_name_i]

names(annot_description) <- gsub(text_for_list1, "", names(annot_description))
names(annot_domain_name) <- gsub(text_for_list2, "", names(annot_domain_name))

names(unnannot_description) <- gsub(text_for_list1, "", names(unnannot_description))
names(unnannot_domain_name) <- gsub(text_for_list2, "", names(unnannot_domain_name))

updated_matching_objects_annot <- list()

for (name in names(annot_description)) {
  updated_matching_objects_annot[[name]] <- c(annot_description[[name]], annot_domain_name[[name]])
}

updated_matching_objects_unnannot <- list()

for (name in names(unnannot_description)) {
  updated_matching_objects_unnannot[[name]] <- c(unnannot_description[[name]], unnannot_domain_name[[name]])
}


merged_list <- list()

for (name in names(updated_matching_objects_annot)) {
  merged_list[[name]] <- c(updated_matching_objects_annot[[name]], updated_matching_objects_unnannot[[name]])
}

text_for_list3 <- "unique_"

names(merged_list) <- gsub(text_for_list3, "", names(merged_list))

original_list <- merged_list

merged_list <- Filter(function(x) any(x != 0), merged_list)

extract_rows <- function(merged_list, domains) {
  
  extracted_list <- list()
  
  for (name in names(merged_list)) {
    extracted_rows <- domains[domains$Protein_id %in% merged_list[[name]], ]
    extracted_list[[name]] <- extracted_rows
  }
  
  return(extracted_list)
}


extracted_list <- extract_rows(merged_list, domains)



save_as_excel <- function(extracted_list, file_name) {
  
  wb <- createWorkbook()
  
  for (name in names(extracted_list)) {
    addWorksheet(wb, sheetName = name)
    
    writeData(wb, sheet = name, extracted_list[[name]])
  }
  
  saveWorkbook(wb, file_name)
}

excel_file <- "extracted_domains.xlsx"

save_as_excel(extracted_list, excel_file)


combined_vector <- unlist(merged_list)
unique_list_terms <- unique(combined_vector)
num_unique__list_terms <- length(unique_list_terms)

#Sequence extraction to blast
domains2seq <- c("XP_057308713.1", "XP_057306849.1", "XP_057296087.1", "XP_057292180.1", "XP_057309179.1", "XP_057298019.1", "XP_057290494.1", "XP_057289544.1", "XP_057312802.1", "XP_057303477.1", "XP_057315269.1", "XP_057307724.1", "XP_057297863.1", "XP_057293581.1", "XP_057313813.1", "XP_057316842.1", "XP_057298326.1", "XP_057311725.1", "XP_057313020.1", "XP_057290626.1", "XP_057297702.1", "XP_057294781.1", "XP_057308464.1", "XP_057295560.1", "XP_057305314.1", "XP_057304961.1", "XP_057301477.1", "XP_057298407.1", "XP_057303477.1", "XP_057296629.1", "XP_057292664.1")
domains2seq <- unique(domains2seq)
ids2seq <- degs[degs$protein_id %in% domains2seq, ]

fastap <- "~/2023-II/Estadistica Genómica/ncbi_dataset/ncbi_dataset/data/GCF_029227915.1/protein.faa"
protein_sequences <- readAAStringSet(fastap)
pids_in <- ids2seq$protein_id
sequence_names <- names(protein_sequences)
extracted_pids <- gsub("^([^ ]+) .*", "\\1", sequence_names)
matching_indices_in <- match(pids_in, extracted_pids)
selected_proteins_in <- protein_sequences[matching_indices_in]

Interesting_proteins <- "Interesting_final.fasta"
writeXStringSet(selected_proteins_in, file = Interesting_proteins)

#Heatmap
heat_rows <- ids2seq[, c("gene_id", "protein_id")]
row_lab <- data.frame(Column_Name = ids2seq$protein_id)
row_labels <- heat_rows$gene_id
rownames(row_lab) <- row_labels


select_rows <- match(heat_rows$gene_id, rownames(assay(rld)))
select_rows <- select_rows[!is.na(select_rows)]

heat_df <- as.data.frame(colData(dds)["condition"])

heatmap_colors <- colorRampPalette(c("blue", "ivory", "red"))(30)
annotation_colors <- list(condition = c("Control" = "#0099FF", "Treatment" = "#FF0090"))
annotation_row <- data.frame(protein_id = heat_rows$protein_id)


ph <- pheatmap(assay(rld)[select_rows,], cluster_rows = TRUE, show_rownames = TRUE,
               cluster_cols = TRUE, 
               annotation_col = heat_df, 
               labels_row = row_lab$Column_Name, 
               color = heatmap_colors, 
               annotation_colors = annotation_colors)

ph


#Coexpression network

gsg <- goodSamplesGenes(t(mcounts))
summary(gsg)
gsg$allOK
table(gsg$goodGenes)
table(gsg$goodSamples)

Wdds <- DESeqDataSetFromMatrix(countData = mcounts,
                               colData = metadata,
                               design = ~ 1)

Wdds75 <- Wdds[rowSums(counts(Wdds) >= 15) >= 8,]
nrow(Wdds75)

Wdds_norm <- vst(Wdds75)

tcounts <- assay(Wdds_norm) %>% 
  t()
power <- c(c(1:10), seq(from = 12, to = 50, by = 2))

sft <- pickSoftThreshold(tcounts,
                         powerVector = power,
                         networkType = "signed",
                         verbose = 5)

sft_data <- sft$fitIndices

a1 <- ggplot(sft_data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = "#C21807") + 
  labs(x = "Power", y = "Scale free topology model fit, signed R^2") +
  theme_classic()

a2 <- ggplot(sft_data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = "Power", y = "Mean Conectivity") +
  theme_classic()

grid.arrange(a1, a2, ncol = 2)

tcounts[] <- sapply(tcounts, as.numeric)
soft_power <- 20
temp_cor <- cor
cor <- WGCNA::cor

bwnet <- blockwiseModules(tcounts,
                          maxBlockSize = 15000,
                          TOMType = "signed",
                          power = soft_power,
                          mergeCutHeight = 0.25,
                          numericLabels = FALSE,
                          randomSeed = 1234,
                          saveTOMs = TRUE,
                          verbose = 3)
cor <- temp_cor
module_eigengenes <- bwnet$MEs
head(module_eigengenes)
table(bwnet$colors)

plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwnet$unmergedColors, bwnet$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)

colData <- metadata

traits <- colData %>%
  mutate(state_bin = ifelse(grepl("Treatment", condition), 1, 0)) %>%
  dplyr::select(2)

nsamples <- nrow(tcounts)
ngenes <- ncol(tcounts)

module_trair_corr <- cor(module_eigengenes, traits, use = 'p')
module_trair_corr_p <- corPvalueStudent(module_trair_corr, nsamples)

heatmap_data <- merge(module_eigengenes, traits, by = 'row.names')
head(heatmap_data)

heatmap_data <- heatmap_data %>%
  column_to_rownames(var = 'Row.names')
names(heatmap_data)

CLP <- CorLevelPlot(heatmap_data,
                    x = names(heatmap_data)[33],
                    y = names(heatmap_data)[1:32],
                    col = c("blue", "skyblue", "ivory", "pink", "red"))

module_gene_mapping <- as.data.frame(bwnet$colors)

#Positive

turquoise <- module_gene_mapping %>% 
  filter(`bwnet$colors` == "turquoise")
red <- module_gene_mapping %>% 
  filter(`bwnet$colors` == "red")
tan <- module_gene_mapping %>% 
  filter(`bwnet$colors` == "tan")

#Negative

yellow <- module_gene_mapping %>% 
  filter(`bwnet$colors` == "yellow")
blue <- module_gene_mapping %>% 
  filter(`bwnet$colors` == "blue")
darkred <- module_gene_mapping %>% 
  filter(`bwnet$colors` == "darkred")

module_membership_measure <- cor(module_eigengenes, tcounts, use = 'p')
module_membership_measure_pvals <- corPvalueStudent(module_membership_measure, nsamples)

gene_signf_cor <- cor(tcounts, traits$state_bin, use = 'p')
gene_signf_cor_pvals <- corPvalueStudent(gene_signf_cor, nsamples)


#Intersections

#Postives

turquoise2degs <- rownames(turquoise)
turquoise_genes <- degs[degs$gene_id %in% turquoise2degs, ]
turquoise_gene2domain <- turquoise_genes$protein_id
turquoise_domains <- domains[domains$Protein_id %in% turquoise_gene2domain, ]

red2degs <- rownames(red)
red_genes <- degs[degs$gene_id %in% red2degs, ]
red_gene2domain <- red_genes$protein_id
red_domains <- domains[domains$Protein_id %in% red_gene2domain, ]

tan2degs <- rownames(tan)
tan_genes <- degs[degs$gene_id %in% tan2degs, ]
tan_gene2domain <- tan_genes$protein_id
tan_domains <- domains[domains$Protein_id %in% tan_gene2domain, ]

#Negatives

yellow2degs <- rownames(yellow)
yellow_genes <- degs[degs$gene_id %in% yellow2degs, ]
yellow_gene2domain <- yellow_genes$protein_id
yellow_domains <- domains[domains$Protein_id %in% yellow_gene2domain, ]

blue2degs <- rownames(blue)
blue_genes <- degs[degs$gene_id %in% blue2degs, ]
blue_gene2domain <- blue_genes$protein_id
blue_domains <- domains[domains$Protein_id %in% blue_gene2domain, ]

darkred2degs <- rownames(darkred)
darkred_genes <- degs[degs$gene_id %in% darkred2degs, ]
darkred_gene2domain <- darkred_genes$protein_id
darkred_domains <- domains[domains$Protein_id %in% darkred_gene2domain, ]


#Export to cytoscape

adjMat <- as.matrix(TOM)

edgeFile <- "Input-edges.txt"
nodeFile <- "Input-nodes.txt"

exportNetworkToCytoscape(adjMat,
                         edgeFile = edgeFile,
                         nodeFile = nodeFile,
                         weighted = TRUE,
                         threshold = 0.25,
                         nodeNames = dimnames(adjMat)[[1]],
                         altNodeNames = NULL,
                         nodeAttr = bwnet$colors,
                         includeColNames = TRUE)


#Venn

degs2net <- degs$gene_id

venn_net_lists <- list(DEGs = degs2net, Turquoise = turquoise2degs, Red = red2degs, Tan = tan2degs)
venn_net <- Venn(venn_net_lists)
venn_net_data <- process_data(venn_net)
venn_plot_net <- ggplot() +
  geom_polygon(aes(X, Y, fill = name, group = id), 
               data = venn_regionedge(venn_net_data),
               show.legend = FALSE) +
  geom_path(aes(X, Y, color = id, group = id), 
            data = venn_setedge(venn_net_data), 
            linewidth = 1,
            show.legend = FALSE) +
  geom_text(aes(X, Y, label = name),
            size = 7,
            data = venn_setlabel(venn_net_data)) +
  geom_label(aes(X, Y, label = count),
             size = 8,
             data = venn_regionlabel(venn_net_data),
             alpha = 0) +
  scale_fill_manual(values = c("DEGs" = "#F8766D", "Turquoise" = "#40E0D0", "Red" = "#FF0000", "Tan" = "#D2B48C", "DEGs/Turquoise" = "#9cab9e", "DEGs/Red" = "#fb3b36", "DEGs/Tan" = "#e5957c", "Turquoise/Red" = "#9f7068", "Turquoise/Tan" = "#89caae", "Red/Tan" = "#e85a46", "DEGs/Turquoise/Red" = "#cd554f", "DEGs/Turquoise/Tan" = "#b7af95", "DEGs/Red/Tan" = "#e67761", "Turquoise/Red/Tan" = "#b8927a", "DEGs/Turquoise/Red/Tan" = "#c28272")) + 
  scale_color_manual(values = c("1"= "black", "2" = "black", "3" = "black", "4" = "black")) +
  coord_equal() +
  theme_void() 

venn_plot_net


venn_net_neg_lists <- list(DEGs = degs2net, Yellow = yellow2degs, Blue = blue2degs, DarkRed = darkred2degs)
venn_net_neg <- Venn(venn_net_neg_lists)
venn_net_data_neg <- process_data(venn_net_neg)
venn_plot_net_neg <- ggplot() +
  geom_polygon(aes(X, Y, fill = name, group = id), 
               data = venn_regionedge(venn_net_data_neg),
               show.legend = FALSE) +
  geom_path(aes(X, Y, color = id, group = id), 
            data = venn_setedge(venn_net_data_neg), 
            linewidth = 1,
            show.legend = FALSE) +
  geom_text(aes(X, Y, label = name),
            size = 7,
            data = venn_setlabel(venn_net_data_neg)) +
  geom_label(aes(X, Y, label = count),
             size = 8,
             data = venn_regionlabel(venn_net_data_neg),
             alpha = 0) +
  scale_fill_manual(values = c("DEGs" = "#F8766D", "Yellow" = "#FFFF00", "Blue" = "#0000FF", "DarkRed" = "#8b0000", "DEGs/Yellow" = "#fbba36", "DEGs/Blue" = "#7c3bb6", "DEGs/DarkRed" = "#c13b36", "Yellow/Blue" = "#7f7f7f", "Yellow/DarkRed" = "#c57f00", "Blue/DarkRed" = "#45007f", "DEGs/Yellow/Blue" = "#7d5d9a", "DEGs/Yellow/DarkRed" = "#c35d1b", "DEGs/Blue/DarkRed" = "#831d5b", "Yellow/Blue/DarkRed" = "#853f3f", "DEGs/Yellow/Blue/DarkRed" = "#a05d5a")) + 
  scale_color_manual(values = c("1"= "black", "2" = "black", "3" = "black", "4" = "black")) +
  coord_equal() +
  theme_void() 

venn_plot_net_neg


turquoise_red <- union(turquoise2degs, red2degs)
turquoise_red <- unique(turquoise_red)
related_degs <- union(turquoise_red, tan2degs)
related_degs <- unique(related_degs)

yellow_blue <- union(yellow2degs, blue2degs)
yellow_blue <- unique(yellow_blue)
unrelated_degs <- union(yellow_blue, darkred2degs)
unrelated_degs <- unique(unrelated_degs)

venn_lists_rela <- list(Related = related_degs, Unrelated = unrelated_degs)
venn_rela <- Venn(venn_lists_rela)
venn_data_rela <- process_data(venn_rela)
venn_plot_rela <- ggplot() +
  geom_polygon(aes(Y, X, fill = name, group = id), 
               data = venn_regionedge(venn_data_rela),
               show.legend = FALSE) +
  geom_path(aes(Y, X, color = id, group = id), 
            data = venn_setedge(venn_data_rela), 
            linewidth = 1,
            show.legend = FALSE) +
  geom_text(aes(Y, X, label = name),
            size = 7,
            data = venn_setlabel(venn_data_rela)) +
  geom_label(aes(Y, X, label = count),
             size = 8,
             data = venn_regionlabel(venn_data_rela),
             alpha = 0) +
  scale_fill_manual(values = c("Related" = "#9bf725", "Unrelated" = "#f5443b", "Related/Unrelated" = "#c89d30")) + 
  scale_color_manual(values = c("1"= "black", "2" = "black")) +
  coord_equal() +
  theme_void() 

venn_plot_rela



#Network + Domains


not_related2IDCP <- degs[degs$gene_id %in% unrelated_degs, ]
not_related_IDCP <- not_related2IDCP[not_related2IDCP$protein_id %in% domains2seq, ]

related2IDCP <- degs[degs$gene_id %in% related_degs, ]
related_IDCP <- related2IDCP[related2IDCP$protein_id %in% domains2seq, ]

venn_lists_dom <- list(RDEGs = related2IDCP$protein_id, UDEGs = not_related2IDCP$protein_id, IDCPs = domains2seq)
venn_dom <- Venn(venn_lists_dom)
venn_data_dom <- process_data(venn_dom)
venn_plot_dom <- ggplot() +
  geom_polygon(aes(Y, X, fill = name, group = id), 
               data = venn_regionedge(venn_data_dom),
               show.legend = FALSE) +
  geom_path(aes(Y, X, color = id, group = id), 
            data = venn_setedge(venn_data_dom), 
            linewidth = 1,
            show.legend = FALSE) +
  geom_text(aes(Y, X, label = name),
            size = 7,
            data = venn_setlabel(venn_data_dom)) +
  geom_label(aes(Y, X, label = count),
             size = 8,
             data = venn_regionlabel(venn_data_dom),
             alpha = 0) +
  scale_fill_manual(values = c("RDEGs" = "#17b700", "UDEGs" = "#b71900", "IDCPs" = "#ee8037", "RDEGs/UDEGs" = "#676800", "RDEGs/IDCPs" = "#829b1b", "UDEGs/IDCPs" = "#d24c1b", "RDEGs/UDEGs/IDCPs" = "#9c5a0d")) + 
  scale_color_manual(values = c("1"= "black", "2" = "black", "3" = "black")) +
  coord_equal() +
  theme_void() 

venn_plot_dom

related_not_related <- not_related2IDCP[not_related2IDCP$protein_id %in% related2IDCP$protein_id, ]

t_mod_cot <- t(module_membership_measure)
t_mod_p <- t(module_membership_measure_pvals)

related_not_related_mc <- t_mod_cot[related_not_related$gene_id, , drop = FALSE ]
related_not_related_p <- t_mod_p[related_not_related$gene_id, , drop = FALSE ]

rnr_mc <- related_not_related_mc[, c(6, 7, 21)]
rnr_p <- related_not_related_p[, c(6, 7, 21)]

turquoise_idcp <- turquoise_genes[turquoise_genes$gene_id %in% ids2seq$gene_id, ] 
red_idcp <- red_genes[red_genes$gene_id %in% ids2seq$gene_id, ] 
tan_idcp <- tan_genes[tan_genes$gene_id %in% ids2seq$gene_id, ] 
yellow_idcp <- yellow_genes[yellow_genes$gene_id %in% ids2seq$gene_id, ] 
blue_idcp <- blue_genes[blue_genes$gene_id %in% ids2seq$gene_id, ] 
darkred_idcp <- darkred_genes[darkred_genes$gene_id %in% ids2seq$gene_id, ] 

#Related threes

alignment_XP_057292180 <- read.phyDat('~/IEI/Maestria/Tesis/Tesis tras burnout/Tesis/Final Tesis/Align/XP_057292180.fas', format = 'fasta', type = 'AA')
dist_XP_057292180 <- dist.ml(alignment_XP_057292180)
tree_XP_057292180 <- NJ(dist_XP_057292180)
fit_XP_057292180 <- pml(tree_XP_057292180, data = alignment_XP_057292180)

alignment_XP_057294781 <- read.phyDat('~/IEI/Maestria/Tesis/Tesis tras burnout/Tesis/Final Tesis/Align/XP_057294781.fas', format = 'fasta', type = 'AA')
dist_XP_057294781 <- dist.ml(alignment_XP_057294781)
tree_XP_057294781 <- NJ(dist_XP_057294781)
fit_XP_057294781 <- pml(tree_XP_057294781, data = alignment_XP_057294781)

alignment_XP_057303477 <- read.phyDat('~/IEI/Maestria/Tesis/Tesis tras burnout/Tesis/Final Tesis/Align/XP_057303477.fas', format = 'fasta', type = 'AA')
dist_XP_057303477 <- dist.ml(alignment_XP_057303477)
tree_XP_057303477 <- NJ(dist_XP_057303477)
fit_XP_057303477 <- pml(tree_XP_057303477, data = alignment_XP_057303477)

alignment_XP_057304961 <- read.phyDat('~/IEI/Maestria/Tesis/Tesis tras burnout/Tesis/Final Tesis/Align/XP_057304961.fas', format = 'fasta', type = 'AA')
dist_XP_057304961 <- dist.ml(alignment_XP_057304961)
tree_XP_057304961 <- NJ(dist_XP_057304961)
fit_XP_057304961 <- pml(tree_XP_057304961, data = alignment_XP_057304961)

alignment_XP_057308464 <- read.phyDat('~/IEI/Maestria/Tesis/Tesis tras burnout/Tesis/Final Tesis/Align/XP_057308464.fas', format = 'fasta', type = 'AA')
dist_XP_057308464 <- dist.ml(alignment_XP_057308464)
tree_XP_057308464 <- NJ(dist_XP_057308464)
fit_XP_057308464 <- pml(tree_XP_057308464, data = alignment_XP_057308464)

alignment_XP_057312802 <- read.phyDat('~/IEI/Maestria/Tesis/Tesis tras burnout/Tesis/Final Tesis/Align/XP_057312802.fas', format = 'fasta', type = 'AA')
dist_XP_057312802 <- dist.ml(alignment_XP_057312802)
tree_XP_057312802 <- NJ(dist_XP_057312802)
fit_XP_057312802 <- pml(tree_XP_057312802, data = alignment_XP_057312802)

alignment_XP_057313020 <- read.phyDat('~/IEI/Maestria/Tesis/Tesis tras burnout/Tesis/Final Tesis/Align/XP_057313020.fas', format = 'fasta', type = 'AA')
dist_XP_057313020 <- dist.ml(alignment_XP_057313020)
tree_XP_057313020 <- NJ(dist_XP_057313020)
fit_XP_057313020 <- pml(tree_XP_057313020, data = alignment_XP_057313020)


mt_XP_057292180 <- modelTest(fit_XP_057292180, model=c("WAG", "JTT", "LG"), control = pml.control(trace = 0))
mt_XP_057294781 <- modelTest(fit_XP_057294781, model=c("WAG", "JTT", "LG"), control = pml.control(trace = 0))
mt_XP_057303477 <- modelTest(fit_XP_057303477, model=c("WAG", "JTT", "LG"), control = pml.control(trace = 0))
mt_XP_057304961 <- modelTest(fit_XP_057304961, model=c("WAG", "JTT", "LG"), control = pml.control(trace = 0))
mt_XP_057308464 <- modelTest(fit_XP_057308464, model=c("WAG", "JTT", "LG"), control = pml.control(trace = 0))
mt_XP_057312802 <- modelTest(fit_XP_057312802, model=c("WAG", "JTT", "LG"), control = pml.control(trace = 0))
mt_XP_057313020 <- modelTest(fit_XP_057313020, model=c("WAG", "JTT", "LG"), control = pml.control(trace = 0))


#Non - Related threes

alignment_XP_057290626 <- read.phyDat('~/IEI/Maestria/Tesis/Tesis tras burnout/Tesis/Final Tesis/Align/XP_057290626.fas', format = 'fasta', type = 'AA')
dist_XP_057290626 <- dist.ml(alignment_XP_057290626)
tree_XP_057290626 <- NJ(dist_XP_057290626)
fit_XP_057290626 <- pml(tree_XP_057290626, data = alignment_XP_057290626)

alignment_XP_057293581 <- read.phyDat('~/IEI/Maestria/Tesis/Tesis tras burnout/Tesis/Final Tesis/Align/XP_057293581.fas', format = 'fasta', type = 'AA')
dist_XP_057293581 <- dist.ml(alignment_XP_057293581)
tree_XP_057293581 <- NJ(dist_XP_057293581)
fit_XP_057293581 <- pml(tree_XP_057293581, data = alignment_XP_057293581)

alignment_XP_057297702 <- read.phyDat('~/IEI/Maestria/Tesis/Tesis tras burnout/Tesis/Final Tesis/Align/XP_057297702.fas', format = 'fasta', type = 'AA')
dist_XP_057297702 <- dist.ml(alignment_XP_057297702)
tree_XP_057297702 <- NJ(dist_XP_057297702)
fit_XP_057297702 <- pml(tree_XP_057297702, data = alignment_XP_057297702)

alignment_XP_057298019 <- read.phyDat('~/IEI/Maestria/Tesis/Tesis tras burnout/Tesis/Final Tesis/Align/XP_057298019.fas', format = 'fasta', type = 'AA')
dist_XP_057298019 <- dist.ml(alignment_XP_057298019)
tree_XP_057298019 <- NJ(dist_XP_057298019)
fit_XP_057298019 <- pml(tree_XP_057298019, data = alignment_XP_057298019)

alignment_XP_057298407 <- read.phyDat('~/IEI/Maestria/Tesis/Tesis tras burnout/Tesis/Final Tesis/Align/XP_057298407.fas', format = 'fasta', type = 'AA')
dist_XP_057298407 <- dist.ml(alignment_XP_057298407)
tree_XP_057298407 <- NJ(dist_XP_057298407)
fit_XP_057298407 <- pml(tree_XP_057298407, data = alignment_XP_057298407)

alignment_XP_057308713 <- read.phyDat('~/IEI/Maestria/Tesis/Tesis tras burnout/Tesis/Final Tesis/Align/XP_057308713.fas', format = 'fasta', type = 'AA')
dist_XP_057308713 <- dist.ml(alignment_XP_057308713)
tree_XP_057308713 <- NJ(dist_XP_057308713)
fit_XP_057308713 <- pml(tree_XP_057308713, data = alignment_XP_057308713)

alignment_XP_057311725 <- read.phyDat('~/IEI/Maestria/Tesis/Tesis tras burnout/Tesis/Final Tesis/Align/XP_057311725.fas', format = 'fasta', type = 'AA')
dist_XP_057311725 <- dist.ml(alignment_XP_057311725)
tree_XP_057311725 <- NJ(dist_XP_057311725)
fit_XP_057311725 <- pml(tree_XP_057311725, data = alignment_XP_057311725)

alignment_XP_057315269 <- read.phyDat('~/IEI/Maestria/Tesis/Tesis tras burnout/Tesis/Final Tesis/Align/XP_057315269.fas', format = 'fasta', type = 'AA')
dist_XP_057315269 <- dist.ml(alignment_XP_057315269)
tree_XP_057315269 <- NJ(dist_XP_057315269)
fit_XP_057315269 <- pml(tree_XP_057315269, data = alignment_XP_057315269)

alignment_XP_057316842 <- read.phyDat('~/IEI/Maestria/Tesis/Tesis tras burnout/Tesis/Final Tesis/Align/XP_057316842.fas', format = 'fasta', type = 'AA')
dist_XP_057316842 <- dist.ml(alignment_XP_057316842)
tree_XP_057316842 <- NJ(dist_XP_057316842)
fit_XP_057316842 <- pml(tree_XP_057316842, data = alignment_XP_057316842)

mt_XP_057290626 <- modelTest(fit_XP_057290626, model=c("WAG", "JTT", "LG"), control = pml.control(trace = 0))
mt_XP_057293581 <- modelTest(fit_XP_057293581, model=c("WAG", "JTT", "LG"), control = pml.control(trace = 0))
mt_XP_057297702 <- modelTest(fit_XP_057297702, model=c("WAG", "JTT", "LG"), control = pml.control(trace = 0))
mt_XP_057298019 <- modelTest(fit_XP_057298019, model=c("WAG", "JTT", "LG"), control = pml.control(trace = 0))
mt_XP_057298407 <- modelTest(fit_XP_057298407, model=c("WAG", "JTT", "LG"), control = pml.control(trace = 0))
mt_XP_057308713 <- modelTest(fit_XP_057308713, model=c("WAG", "JTT", "LG"), control = pml.control(trace = 0))
mt_XP_057311725 <- modelTest(fit_XP_057311725, model=c("WAG", "JTT", "LG"), control = pml.control(trace = 0))
mt_XP_057315269 <- modelTest(fit_XP_057315269, model=c("WAG", "JTT", "LG"), control = pml.control(trace = 0))
mt_XP_057316842 <- modelTest(fit_XP_057316842, model=c("WAG", "JTT", "LG"), control = pml.control(trace = 0))

#Model test
combined_mt <- bind_rows(mt_XP_057290626, mt_XP_057292180, mt_XP_057293581, mt_XP_057294781, mt_XP_057297702, mt_XP_057298019, mt_XP_057298407, mt_XP_057303477, mt_XP_057304961, mt_XP_057308464, mt_XP_057308713, mt_XP_057311725, mt_XP_057312802, mt_XP_057313020, mt_XP_057315269, mt_XP_057316842, .id = "source")

result_mt <- combined_mt %>%
  group_by(Model) %>%
  summarize(across(where(is.numeric), mean))

#Related

fit_XP_057292180_tree <- optim.pml(fit_XP_057292180, model = "WAG")
fit_XP_057294781_tree <- optim.pml(fit_XP_057294781, model = "WAG")
fit_XP_057303477_tree <- optim.pml(fit_XP_057303477, model = "WAG")
fit_XP_057304961_tree <- optim.pml(fit_XP_057304961, model = "WAG")
fit_XP_057308464_tree <- optim.pml(fit_XP_057308464, model = "WAG")
fit_XP_057312802_tree <- optim.pml(fit_XP_057312802, model = "WAG")
fit_XP_057313020_tree <- optim.pml(fit_XP_057313020, model = "WAG")

#Non - related

fit_XP_057290626_tree <- optim.pml(fit_XP_057290626, model = "WAG")
fit_XP_057293581_tree <- optim.pml(fit_XP_057293581, model = "WAG")
fit_XP_057297702_tree <- optim.pml(fit_XP_057297702, model = "WAG")
fit_XP_057298019_tree <- optim.pml(fit_XP_057298019, model = "WAG")
fit_XP_057298407_tree <- optim.pml(fit_XP_057298407, model = "WAG")
fit_XP_057308713_tree <- optim.pml(fit_XP_057308713, model = "WAG")
fit_XP_057311725_tree <- optim.pml(fit_XP_057311725, model = "WAG")
fit_XP_057315269_tree <- optim.pml(fit_XP_057315269, model = "WAG")
fit_XP_057316842_tree <- optim.pml(fit_XP_057316842, model = "WAG")

#Related

phylo_tree_XP_057292180 <- as.treedata(fit_XP_057292180_tree)
phylo_tree_XP_057294781 <- as.treedata(fit_XP_057294781_tree)
phylo_tree_XP_057303477 <- as.treedata(fit_XP_057303477_tree)
phylo_tree_XP_057304961 <- as.treedata(fit_XP_057304961_tree)
phylo_tree_XP_057308464 <- as.treedata(fit_XP_057308464_tree)
phylo_tree_XP_057312802 <- as.treedata(fit_XP_057312802_tree)
phylo_tree_XP_057313020 <- as.treedata(fit_XP_057313020_tree)

#Non - related

phylo_tree_XP_057290626 <- as.treedata(fit_XP_057290626_tree)
phylo_tree_XP_057293581 <- as.treedata(fit_XP_057293581_tree)
phylo_tree_XP_057297702 <- as.treedata(fit_XP_057297702_tree)
phylo_tree_XP_057298019 <- as.treedata(fit_XP_057298019_tree)
phylo_tree_XP_057298407 <- as.treedata(fit_XP_057298407_tree)
phylo_tree_XP_057308713 <- as.treedata(fit_XP_057308713_tree)
phylo_tree_XP_057311725 <- as.treedata(fit_XP_057311725_tree)
phylo_tree_XP_057315269 <- as.treedata(fit_XP_057315269_tree)
phylo_tree_XP_057316842 <- as.treedata(fit_XP_057316842_tree)



#ggtree

labels <- c("uncharacterized protein LOC130614748 Hydractinia symbiolongicarpus", "uncharacterized protein LOC130613291 isoform X2 Hydractinia symbiolongicarpus", "uncharacterized protein LOC130622166 Hydractinia symbiolongicarpus", "heparanase-like Hydractinia symbiolongicarpus", "basic leucine zipper transcriptional factor ATF-like isoform X2 Hydractinia symbiolongicarpus", "SCO-spondin-like isoform X2 Hydractinia symbiolongicarpus(2)", "matrix metalloproteinase-17-like Hydractinia symbiolongicarpus(2)", "hatching enzyme-like Hydractinia symbiolongicarpus", "zinc metalloproteinase nas-6-like Hydractinia symbiolongicarpus(2)", "uncharacterized protein LOC130646309 Hydractinia symbiolongicarpus", "uncharacterized protein LOC130647021 Hydractinia symbiolongicarpus", "netrin receptor DCC-like isoform X1 Hydractinia symbiolongicarpus", "uncharacterized protein LOC130654269 Hydractinia symbiolongicarpus", "interferon-gamma-inducible GTPase 10-like Hydractinia symbiolongicarpus", "hydralysin-like Hydractinia symbiolongicarpus(2)", "uncharacterized protein LOC130657866 Hydractinia symbiolongicarpus")

#Related

ggtree(phylo_tree_XP_057292180, ladderize = "right") +
  geom_tiplab(aes(color = label %in% labels), align = FALSE) +
  scale_color_manual(values=c("black", "red")) +
  theme_tree2() +
  theme(legend.position = "none") +
  xlim_tree(5)

ggtree(phylo_tree_XP_057294781, ladderize = "right") +
  geom_tiplab(aes(color = label %in% labels), align = FALSE) +
  scale_color_manual(values=c("black", "red")) +
  theme_tree2() +
  theme(legend.position = "none") +
  xlim_tree(5)

ggtree(phylo_tree_XP_057303477, ladderize = "right") +
  geom_tiplab(aes(color = label %in% labels), align = FALSE) +
  scale_color_manual(values=c("black", "red")) +
  theme_tree2() +
  theme(legend.position = "none") +
  xlim_tree(5)

ggtree(phylo_tree_XP_057304961, ladderize = "right") +
  geom_tiplab(aes(color = label %in% labels), align = FALSE) +
  scale_color_manual(values=c("black", "red")) +
  theme_tree2() +
  theme(legend.position = "none") +
  xlim_tree(5)

ggtree(phylo_tree_XP_057308464, ladderize = "right") +
  geom_tiplab(aes(color = label %in% labels), align = FALSE) +
  scale_color_manual(values=c("black", "red")) +
  theme_tree2() +
  theme(legend.position = "none") +
  xlim_tree(5)

ggtree(phylo_tree_XP_057312802, ladderize = "right") +
  geom_tiplab(aes(color = label %in% labels), align = FALSE) +
  scale_color_manual(values=c("black", "red")) +
  theme_tree2() +
  theme(legend.position = "none") +
  xlim_tree(5)

ggtree(phylo_tree_XP_057313020, ladderize = "right") +
  geom_tiplab(aes(color = label %in% labels), align = FALSE) +
  scale_color_manual(values=c("black", "red")) +
  theme_tree2() +
  theme(legend.position = "none") +
  xlim_tree(5)

#Non - related

ggtree(phylo_tree_XP_057290626, ladderize = "right") +
  geom_tiplab(aes(color = label %in% labels), align = FALSE) +
  scale_color_manual(values=c("black", "red")) +
  theme_tree2() +
  theme(legend.position = "none") +
  xlim_tree(5)

ggtree(phylo_tree_XP_057293581, ladderize = "right") +
  geom_tiplab(aes(color = label %in% labels), align = FALSE) +
  scale_color_manual(values=c("black", "red")) +
  theme_tree2() +
  theme(legend.position = "none") +
  xlim_tree(5)

ggtree(phylo_tree_XP_057297702, ladderize = "right") +
  geom_tiplab(aes(color = label %in% labels), align = FALSE) +
  scale_color_manual(values=c("black", "red")) +
  theme_tree2() +
  theme(legend.position = "none") +
  xlim_tree(5)

ggtree(phylo_tree_XP_057298019, ladderize = "right") +
  geom_tiplab(aes(color = label %in% labels), align = FALSE) +
  scale_color_manual(values=c("black", "red")) +
  theme_tree2() +
  theme(legend.position = "none") +
  xlim_tree(5)

ggtree(phylo_tree_XP_057298407, ladderize = "right") +
  geom_tiplab(aes(color = label %in% labels), align = FALSE) +
  scale_color_manual(values=c("black", "red")) +
  theme_tree2() +
  theme(legend.position = "none") +
  xlim_tree(5)

ggtree(phylo_tree_XP_057308713, ladderize = "right") +
  geom_tiplab(aes(color = label %in% labels), align = FALSE) +
  scale_color_manual(values=c("black", "red")) +
  theme_tree2() +
  theme(legend.position = "none") +
  xlim_tree(5)

ggtree(phylo_tree_XP_057311725, ladderize = "right") +
  geom_tiplab(aes(color = label %in% labels), align = FALSE) +
  scale_color_manual(values=c("black", "red")) +
  theme_tree2() +
  theme(legend.position = "none") +
  xlim_tree(5)

ggtree(phylo_tree_XP_057315269, ladderize = "right") +
  geom_tiplab(aes(color = label %in% labels), align = FALSE) +
  scale_color_manual(values=c("black", "red")) +
  theme_tree2() +
  theme(legend.position = "none") +
  xlim_tree(5)

ggtree(phylo_tree_XP_057316842, ladderize = "right") +
  geom_tiplab(aes(color = label %in% labels), align = FALSE) +
  scale_color_manual(values=c("black", "red")) +
  theme_tree2() +
  theme(legend.position = "none") +
  xlim_tree(5)


#Save as excel

wb <- createWorkbook()

GSEA_GO <- yy_df
GSEA_InterPro <- IP2

addWorksheet(wb, "Sheet1")
writeData(wb, sheet = "Sheet1", GSEA_GO)
addWorksheet(wb, "Sheet2")
writeData(wb, sheet = "Sheet2", GSEA_InterPro)

saveWorkbook(wb, "Hysym_GSEA.xlsx")

#Domains Master Table

wb <- createWorkbook()

for (i in seq_along(table_list)) {
  
  addWorksheet(wb, sheetName = Sheet_Names[i])
  
  writeData(wb, sheet = i, x = table_list[[i]])
}

saveWorkbook(wb, "Domains_Master_Table.xlsx")

#ALT 

wb_05 <- createWorkbook()

for (i in seq_along(table_list_05)) {
  
  addWorksheet(wb_05, sheetName = Sheet_Names_05[i])
  
  writeData(wb_05, sheet = i, x = table_list_05[[i]])
}

saveWorkbook(wb_05, "Domains_Master_Table_05.xlsx")
