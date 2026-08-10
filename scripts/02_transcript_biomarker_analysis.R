rm(list=ls())
gc()

library(edgeR)
library(limma)
library(IsoformSwitchAnalyzeR)
library(glmnet)
library(caret)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(rtracklayer)
library(reshape2)
library(biomaRt)

gene_counts <- read.csv("results/stringtie/gene_count_matrix.csv",
                        row.names = 1, check.names = FALSE)

tx_counts <- read.csv("results/stringtie/transcript_count_matrix.csv",
                      row.names = 1, check.names = FALSE)

group <- factor(c(rep("Viral",5),
                  rep("Healthy",5),
                  rep("Bacterial",5)))
dge_gene <- DGEList(counts = gene_counts, group = group)
keep_gene <- filterByExpr(dge_gene)
dge_gene <- dge_gene[keep_gene,, keep.lib.sizes=FALSE]
dge_gene <- calcNormFactors(dge_gene)

design <- model.matrix(~group)
v_gene <- voom(dge_gene, design)

fit_gene <- lmFit(v_gene, design)
fit_gene <- eBayes(fit_gene)
gene_res <- topTable(fit_gene, number=Inf)

write.csv(
  gene_res,
  "results/differential_expression/Gene_Level_DE_Results.csv"
)

dge_tx <- DGEList(counts = tx_counts, group = group)
keep_tx <- filterByExpr(dge_tx)
dge_tx <- dge_tx[keep_tx,, keep.lib.sizes=FALSE]
dge_tx <- calcNormFactors(dge_tx)

v_tx <- voom(dge_tx, design)
fit_tx <- lmFit(v_tx, design)
fit_tx <- eBayes(fit_tx)

tx_res <- topTable(fit_tx, number=Inf)
write.csv(tx_res,"results/differential_expression/Transcript_Level_DE_Results.csv")
cpm_tx <- cpm(dge_tx)

gtf <- import("data/stringtie/merged_transcripts.gtf")
gtf_tx <- gtf[gtf$type == "transcript"]
valid_ids <- unique(gtf_tx$transcript_id)

common_ids <- intersect(rownames(tx_counts), valid_ids)
tx_counts_sync <- tx_counts[common_ids, ]

design_iso <- data.frame(
  sampleID = colnames(tx_counts_sync),
  condition = group
)

switchList <- importRdata(
  isoformCountMatrix   = tx_counts_sync,
  designMatrix         = design_iso,
  isoformExonAnnoation = "data/stringtie/merged_transcripts.gtf",
  ignoreAfterPeriod    = FALSE,
  removeNonConvensionalChr = TRUE
)

switchList <- preFilter(
  switchList,
  geneExpressionCutoff = 1,
  isoformExpressionCutoff = 0.5,
  IFcutoff = 0.01
)

switchList <- isoformSwitchTestDEXSeq(
  switchList,
  reduceToSwitchingGenes = FALSE
)

switch_res <- extractTopSwitches(
  switchList,
  filterForConsequences = FALSE,
  n = Inf
)

write.csv(
  switch_res,
  "results/isoform_switching/Isoform_Switch_Results.csv",
  row.names = FALSE
)

save(
  switchList,
  file = "results/isoform_switching/switchList_processed.RData"
)

ml_features <- rownames(tx_res)[1:100]

x <- as.matrix(t(log2(cpm_tx[ml_features, ] + 1)))
y <- factor(group)

n_bootstrap <- 100
selection_count <- numeric(ncol(x))
names(selection_count) <- colnames(x)

set.seed(123)

for(i in 1:n_bootstrap){
  
  idx <- createDataPartition(y, p = 0.8, list = FALSE)
  
  x_sub <- x[idx, ]
  y_sub <- y[idx]
  
  fit <- tryCatch({
    cv.glmnet(x_sub, y_sub, family="multinomial", alpha=0.5)
  }, error=function(e) NULL)
  
  if(!is.null(fit)){
    co <- coef(fit, s="lambda.min")
    
    selected <- unique(unlist(lapply(co, function(m){
      rownames(m)[as.numeric(m) != 0]
    })))
    
    selection_count[names(selection_count) %in% selected] <-
      selection_count[names(selection_count) %in% selected] + 1
  }
}

stability_results <- data.frame(
  isoform_id = names(selection_count),
  Frequency = (selection_count/n_bootstrap)*100
) %>%
  arrange(desc(Frequency)) %>%
  filter(isoform_id != "(Intercept)")

final_panel <- stability_results %>%
  slice(1:20) %>%
  left_join(unique(switchList$isoformFeatures[, c("isoform_id","gene_name")]),
            by="isoform_id") %>%
  mutate(gene_name = ifelse(is.na(gene_name)|gene_name=="",
                            isoform_id, gene_name))

write.csv(final_panel, "results/biomarkers/Stable_Top20_Biomarkers.csv", row.names=FALSE)

top_ids <- final_panel$isoform_id
heat_matrix <- log2(cpm_tx[top_ids, ] + 1)

annotation_col <- data.frame(Group=group)
rownames(annotation_col) <- colnames(heat_matrix)

png("figures/Top20_Biomarker_Heatmap.png", width=2000, height=2500, res=300)

pheatmap(heat_matrix,
         scale="row",
         annotation_col=annotation_col,
         clustering_method="complete",
         show_rownames=TRUE,
         fontsize_row=8,
         main="Top 20 Isoform Biomarkers")

dev.off()

#Correlation
cor_matrix <- cor(t(heat_matrix), method = "pearson")
cor_matrix[upper.tri(cor_matrix)] <- NA
cor_df <- melt(cor_matrix, na.rm = TRUE)
n_vars <- nrow(cor_matrix)

text_size <- ifelse(n_vars >= 20, 2.8,
                    ifelse(n_vars >= 15, 3.2, 4))

p <- ggplot(cor_df, aes(Var1, Var2, fill = value)) +
  
  geom_tile(color = "white", linewidth = 0.4) +
  
  geom_text(aes(label = sprintf("%.2f", value)),
            size = text_size,
            fontface = "bold") +
  
  scale_fill_gradient2(
    low = "#2C7BB6",
    mid = "white",
    high = "#D7191C",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Pearson\nCorrelation"
  ) +
  
  coord_fixed() + 
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    panel.grid = element_blank(),
    axis.title = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  
  labs(title = "Correlation Between Top 20 Biomarkers")

ggsave("figures/Top20_Correlation_Triangular.png",
       plot = p,
       width = 10,  
       height = 10,  
       dpi = 300)

pca_final <- prcomp(t(heat_matrix), scale.=TRUE)

pca_df <- data.frame(pca_final$x, Group=group)

png("figures/Top20_PCA_Separation.png", width=2100, height=1800, res=300)

ggplot(pca_df, aes(x=PC1, y=PC2, color=Group, fill=Group)) +
  geom_point(size=6, shape=21, color="black", stroke=1.5) +
  stat_ellipse(geom="polygon", alpha=0.2) +
  theme_minimal(base_size=16) +
  labs(title="PCA Separation Using Top 20 Isoforms",
       x=paste0("PC1 (",
                round(summary(pca_final)$importance[2,1]*100,1),"%)"),
       y=paste0("PC2 (",
                round(summary(pca_final)$importance[2,2]*100,1),"%)"))

dev.off()

#Connect to Ensembl BioMart (Human)
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

final_panel <- read.csv("results/biomarkers/Stable_Top20_Biomarkers.csv", stringsAsFactors = FALSE)
final_panel$ensembl_transcript_id <- gsub("\\..*", "", final_panel$isoform_id)

annotations <- getBM(
  attributes = c("ensembl_transcript_id", "external_gene_name", 
                 "transcript_biotype", "description", "entrezgene_id"),
  filters = "ensembl_transcript_id",
  values = final_panel$ensembl_transcript_id,
  mart = mart
)

final_annotated_file <- final_panel %>%
  left_join(annotations, by = "ensembl_transcript_id") %>%
  mutate(
    external_gene_name = ifelse(is.na(external_gene_name) | external_gene_name == "",
                                isoform_id,
                                external_gene_name)
  )

write.csv(final_annotated_file, "results/biomarkers/Final_Top20_Biomarker_Database_Correlation.csv", row.names = FALSE)

functional_info <- getBM(
  attributes = c("external_gene_name", 
                 "description", 
                 "uniprotswissprot", 
                 "entrezgene_description"),
  filters = "external_gene_name",
  values = final_annotated_file$external_gene_name,
  mart = mart
)

functional_info_clean <- functional_info %>%
  group_by(external_gene_name) %>%
  summarize(
    gene_description = first(description),
    uniprot = paste(unique(uniprotswissprot), collapse=", "),
    entrez_desc = first(entrezgene_description)
  )

final_functional_panel <- final_annotated_file %>%
  left_join(functional_info_clean, by = "external_gene_name") %>%
  mutate(
    derived_function = case_when(
      !is.na(entrez_desc) ~ as.character(entrez_desc),
      !is.na(gene_description) ~ as.character(gene_description),
      grepl("MSTRG", isoform_id) ~ "Novel StringTie-assembled isoform (Uncharacterized)",
      TRUE ~ "Hypothetical Protein / Novel Regulatory RNA"
    )
  )

#Save Functional Landscape for Thesis / Supplement
write.csv(final_functional_panel, "results/biomarkers/Biomarker_Functional_Landscape.csv", row.names = FALSE)
