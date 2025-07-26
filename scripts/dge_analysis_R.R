setwd("~/Desktop/sample_counts")
count_files <- c("bonemarrow_5a.s.count.txt","bonemarrow_6a.s.count.txt", "bonemarrow_6b.s.count.txt","placenta_3a.s.count.txt", "placenta_6a.s.count.txt", "placenta_6b.s.count.txt")
count_list <- lapply(file_paths, function(file){read.table(file, header=TRUE, row.names = 1, sep = "\t")})
counts <- do.call(cbind, count_list)
colnames(counts) <- c("bonemarrow_5a","bonemarrow_6a","bonemarrow_6b","placenta_3a", "placenta_6a","placenta_6b")
head(counts)
group <- factor(c("BoneMarrow", "BoneMarrow", "BoneMarrow","Placenta", "Placenta", "Placenta"))
y <- DGEList(counts = counts, group = group)
y
y <- calcNormFactors(y)
plotMDS(y)
y <- estimateDisp(y)
et <- exactTest(y, pair = c("BoneMarrow", "Placenta"))
de <- decideTests(et)
summary(de)
detags <- rownames(y)[as.logical(de)]
plotSmear(et, de.tags=detags)
abline(h=0, col="blue")
topTags(et)
all_genes_df <- as.data.frame(all_genes)
all_genes_df
diffExpGenes <- topTags(et, n=1000, p.value = 0.05)
head(diffExpGenes$table)
sample_data <- data.frame(row.names = colnames(counts), Group = factor(c("BoneMarrow", "BoneMarrow", "BoneMarrow", "Placenta","Placenta", "Placenta")))
data_deseq <- DESeqDataSetFromMatrix(countData = counts, colData = sample_data, design = ~ 1)
data_deseq
nrow(data_deseq)    
data_deseq <- data_deseq[ rowSums(counts(data_deseq)) > 1, ]
data_deseq
rld <- rlog(data_deseq, blind=FALSE)
sampleDists <- dist(t(assay(rld)))
sampleDists
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste( rld$cell_type, rld$dev_stage, rld$replicate, sep="-" )
rownames(sampleDistMatrix) <- sampleLabels
colnames(sampleDistMatrix) <- sampleLabels
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pheatmap(sampleDistMatrix, clustering_distance_rows = sampleDists, clustering_distance_cols = sampleDists, col = colors)
plotPCA(rld, intgroup = "Group")
write.csv(all_genes_df, file = "dge_results_all_genes.csv", row.names = TRUE)
write.csv(all_genes_df, file = "~/Desktop/dge_results_all_genes.csv", row.names = TRUE)
top50_genes <- rownames(all_genes_df[order(all_genes_df$FDR), ])[1:50]
top50_counts <- y$counts[top50_genes, ]
log_cpm <- cpm(y, log=TRUE)[top50_genes, ]
log_cpm_scaled <- t(scale(t(log_cpm)))  
annotation_col <- data.frame(Tissue_Type = group)
rownames(annotation_col) <- colnames(log_cpm_scaled)
upregulated_placenta <- all_genes_df[all_genes_df$logFC > 0 & all_genes_df$FDR < 0.05, ]
placenta_genes <- rownames(upregulated_placenta)
upregulated_bonemarrow <- all_genes_df[all_genes_df$logFC < 0 & all_genes_df$FDR < 0.05, ]
bonemarrow_genes <- rownames(upregulated_bonemarrow)
bonemarrow_genes <- rownames(all_genes_df[all_genes_df$logFC < 0 & all_genes_df$PValue < 0.05, ])
placenta_genes <- rownames(all_genes_df[all_genes_df$logFC > 0 & all_genes_df$PValue < 0.05, ])
> write.table(placenta_genes, file = "~/Desktop/upregulated_placenta_gene_list_final.txt",row.names = FALSE, col.names = FALSE, quote = FALSE, sep = ",")
> write.table(bonemarrow_genes, file = "~/Desktop/upregulated_bonemarrow_gene_list_final.txt",row.names = FALSE, col.names = FALSE, quote = FALSE, sep = ",")
top50_genes <- rownames(all_genes_df[order(all_genes_df$FDR), ])[1:50]
top50_counts <- y$counts[top50_genes, ]
log_cpm <- cpm(y, log=TRUE)[top50_genes, ]
log_cpm_scaled <- t(scale(t(log_cpm)))
annotation_col <- data.frame(Tissue_Type = group)
rownames(annotation_col) <- colnames(log_cpm_scaled)
pheatmap(log_cpm_scaled,annotation_col = annotation_col,show_rownames = TRUE,show_colnames = TRUE,fontsize_row = 8,scale = "none",clustering_distance_rows = "euclidean",clustering_distance_cols = "euclidean",clustering_method = "complete",color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),main = "Top 50 Most Differentially Expressed Genes")
