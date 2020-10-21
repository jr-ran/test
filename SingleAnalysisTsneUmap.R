#!/usr/bin/Rscript

# parameter
suppressMessages(library(getopt))

spec <- matrix(c(
  'qc.rda'     , 'R', 2, 'character', 'Rdata after samples qc, example: single_qc.Rda',
  'od'         , 'o', 2, 'character', 'out dir name',
  'sample'     , 'a', 1, 'character', 'sample name',
  'FDR'        , 'q', 1, 'numeric'  , 'FDR threshold',
  'pvalue'     , 'p', 1, 'numeric'  , 'pvalue', 
  'fold'       , 'f', 1, 'numeric'  , 'fold change threshold',
  'resolution' , 'r', 1, 'numeric'  , 'resolution',
  'min.pct'    , 'm', 1, 'numeric'  , 'gene expression min percent',
  'id'         , 'i', 2, 'character', 'gene id list, example: id_name.list',
  'color'      , 'c', 1, 'character', 'color',
  'size'       , 's', 1, 'numeric'  , 'font size',
  'top.n'      , 'n', 1, 'numeric'  , 'top.n',
  'help'       , 'h', 0, 'logical'  , 'print usage'
), byrow = TRUE, ncol =5)
opt <- getopt(spec)

PrintUsage <- function(){
  cat("ProgramName:AnalysisTsneUmap.R
Version: Version v1.0
Contact: Ranjr <ranjr@biomarker.com.cn> 
Program Date:   2020.9.3 
Description: this program is used to find cluster marker gene ;\n")
  cat(paste(getopt(spec, usage = TRUE), sep = "\n"))
  q()
}

if(!is.null(opt$help)) {PrintUsage()}
if(is.null(opt$qc.rda) && is.null(opt$od) && is.null(opt$id)) {PrintUsage()}
if(is.null(opt$sample)) {PrintUsage()}
# if(is.null(opt$FDR)) {opt$FDR <- 0.1}
# if(is.null(opt$fold)) {opt$fold <- 2}
if(is.null(opt$resolution)) {opt$resolution <- 0.2}
if(is.null(opt$min.pct)) {opt$min.pct <- 0.25}
if(is.null(opt$color)) {opt$color <- 0.25}
if(is.null(opt$size)) {opt$size <- 10}
if(is.null(opt$top.n)) {opt$top.n <- 10}
if(!dir.exists(opt$od)) {dir.create(opt$od)}

# load packages
suppressMessages(library(SingleCellExperiment))
suppressMessages(library(Seurat))
suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(Matrix))
suppressMessages(library(patchwork))
suppressMessages(library(ggplot2))
suppressMessages(library(scales))
suppressMessages(library(cowplot))
suppressMessages(library(circlize))
suppressMessages(library(parallel))
suppressMessages(library(ComplexHeatmap, lib.loc = "/share/nas1/ranjr/packages/3.6"))

time1 <- proc.time()
# parameter set
rda <- opt$qc.rda
od <- opt$od
sample <- opt$sample
FDR <- opt$FDR
pvalue <- opt$pvalue
fold <- opt$fold
res <- opt$resolution
min.pct <- opt$min.pct
ensembleID_symbol <- opt$id
color <- opt$color
size <- opt$size
top.n <- opt$top.n

# function
SavePlot <- function(od, filename, data, width = 6, height = 4, scale = 1.3){
  file.png <- paste(sample, filename, "png", sep = ".")
  file.pdf <- paste(sample, filename, "pdf", sep = ".")
  ggsave(filename = file.path(od, file.png), plot = data, width = width, height = height, scale = scale)
  ggsave(filename = file.path(od, file.pdf), plot = data, width = width, height = height, scale = scale)
}

AddId <- function(data, id_list){
  data.re <- inner_join(data, id_list, by = c("gene" = "Symbol")) %>% 
    dplyr::select(EnsembleId, gene, dplyr::everything())
  return(data.re)
}


load(rda)
ensemble.id.symbol <- read.table(file = ensembleID_symbol, header = FALSE, col.names = c("EnsembleId", "Symbol"))

# Findvariabls, Normalized, scaledata and remove batch effect
single.seurat <- FindVariableFeatures(single.seurat, selection.method = "vst", nfeatures = 2000)
single.seurat <- SCTransform(object = single.seurat, do.scale = TRUE, verbose = FALSE, vars.to.regress = c("nCount_RNA"))
single.seurat <- RunPCA(single.seurat, features = VariableFeatures(single.seurat))
single.seurat <- RunUMAP(single.seurat, features = VariableFeatures(single.seurat))
single.seurat <- RunTSNE(single.seurat, features = VariableFeatures(single.seurat))
single.seurat <- FindNeighbors(single.seurat, dims = 1:20)
single.seurat <- FindClusters(single.seurat, resolution = res) # resolution 
single.seurat[['cellcluster']] <- Idents(single.seurat)[rownames(single.seurat@meta.data)]

# plot
outputdir <- paste(od, "reduction", sep = "/")
if(!dir.exists(outputdir)) {dir.create(outputdir)}
P.umap <- UMAPPlot(single.seurat, label = TRUE) 
P.tsne <- TSNEPlot(single.seurat, label = TRUE) 
SavePlot(od = outputdir, filename = "umap_all", data = P.umap)
SavePlot(od = outputdir, filename = "tsne_all", data = P.tsne)

pca.projection <- single.seurat@reductions$pca@cell.embeddings
cells <- rownames(pca.projection)
pca.out <- data.frame(Cells = cells, pca.projection)
write.table(pca.out, file = file.path(outputdir, paste(sample, "pca_components.xls", sep = ".")), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

tsne.projection <- single.seurat@reductions$tsne@cell.embeddings
cells <- rownames(tsne.projection)
tsne.out <- data.frame(Cells = cells, tsne.projection)
write.table(tsne.out, file = file.path(outputdir, paste(sample, "tsne_components.xls", sep = ".")), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

umap.projection <- single.seurat@reductions$umap@cell.embeddings
cells <- rownames(umap.projection)
umap.out <- data.frame(Cells = cells, umap.projection)
write.table(umap.out, file = file.path(outputdir, paste(sample, "umap_components.xls", sep = ".")), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# diff marker gene in clusters
logfc.threshold <- log2(fold)
diff.exp.all <- FindAllMarkers(object = single.seurat, assay = "SCT", slot = "counts", min.pct = 0.25, logfc.threshold = 0)
diff.exp.all$avg_logFC <- diff.exp.all$avg_logFC/log(2) 
diff.exp.all <- AddId(diff.exp.all, ensemble.id.symbol)
print("diff threshold is:")
print(FDR)
print(pvalue)
print(fold)
if(!is.null(FDR)){
  if(!is.null(fold)){
    diff.exp.all.filter <- diff.exp.all %>% filter(p_val_adj < FDR & avg_logFC > logfc.threshold)
  }else{
    diff.exp.all.filter <- diff.exp.all %>% filter(p_val_adj < FDR)
  }
}
if(!is.null(pvalue)){
  if(!is.null(fold)){
    diff.exp.all.filter <- diff.exp.all %>% filter(p_val < pvalue & avg_logFC > logfc.threshold)
  }else{
    diff.exp.all.filter <- diff.exp.all %>% filter(p_val < pvalue)
  }
}

# top10 marker gene plot
outputdir <- paste(od, paste("top", top.n, sep = ""), sep = "/")
if(!dir.exists(outputdir)) {dir.create(outputdir)}

# plot vln tsane umap
GlobPlotTheme <- function (font.size = size){
  Theme <- theme(title = element_text(size = font.size), 
                   axis.ticks = element_blank(), 
                   axis.line = element_blank(), 
                   axis.title = element_blank(), 
                   axis.text = element_blank())
  return(Theme)
}

parallel::mclapply(sort(unique(Idents(single.seurat))), function(x){
  top10.genes <- diff.exp.all.filter %>% 
    group_by(cluster) %>% filter(., cluster == x) %>%
    top_n(., n = top.n, wt = avg_logFC) %>%
    .[, "gene", drop = TRUE] %>% unique()
  p.tsne <- FeaturePlot(single.seurat, features = top10.genes, reduction = "tsne", cols = c("lightgrey", "purple"), ncol = 5, pt.size = 0.1) & NoLegend() & GlobPlotTheme()
  outfile <- paste(paste("cluster", x, sep = ""), paste("top", top.n, sep = ""), "markerTsne", sep = "_")
  SavePlot(od = outputdir, filename = outfile, data = p.tsne)
  p.umap <- FeaturePlot(single.seurat, features = top10.genes, reduction = "umap", cols = c("lightgrey", "purple"), ncol = 5, pt.size = 0.1) & NoLegend() & GlobPlotTheme()
  outfile <- paste(paste("cluster", x, sep = ""), paste("top", top.n, sep = ""), "markerUmap", sep = "_")
  SavePlot(od = outputdir, filename = outfile, data = p.umap)
  p.vln <- VlnPlot(single.seurat, features = top10.genes, ncol = 5, pt.size = 0) * theme(title = element_text(size = 7), 
      axis.ticks = element_blank(), 
      axis.title = element_blank(), 
      axis.text.y = element_text(size = 7),
	  axis.text.x  = element_text(size = 6, angle = 90))
  outfile <- paste(paste("cluster", x, sep = ""), paste("top", top.n, sep = ""), "markerVln", sep = "_")
  SavePlot(od = outputdir, filename = outfile, data = p.vln)
  return("done")
}, mc.cores = length(unique(Idents(single.seurat))))

# marker gene expression
marker.gene.avgExp.data <- lapply(sort(unique(Idents(single.seurat))), function(x) {
  cells <- WhichCells(object = single.seurat, idents = x)
  t <- as.matrix(apply(single.seurat[['SCT']]@counts[match(diff.exp.all.filter$gene,rownames(single.seurat[['SCT']]@counts)), cells], 1, median), ncol = 1)
  colnames(t) <- x
  return(t)
})
marker.gene.avgExp <- do.call(cbind,marker.gene.avgExp.data)
marker.gene.avgExp %<>% as.data.frame %>% mutate(gene = rownames(.)) %>% AddId(., ensemble.id.symbol)

# output file
outputdir <- paste(od, "clusterDiff", sep = "/")
if(!dir.exists(outputdir)) {dir.create(outputdir)}

# plot heatmap
# P.heatmap <- DoHeatmap(single.seurat, features = diff.exp.all.filter %>% 
                       # group_by(cluster) %>% 
                       # top_n(., n = top.n, wt = avg_logFC) %>%
                       # .[, "gene", drop = TRUE] %>% unique(), label = F, draw.lines = F) + theme(axis.text = element_text(size = 3))			
topmarker <- diff.exp.all.filter %>% group_by(cluster) %>% top_n(., n = top.n, wt = avg_logFC) %>% .[, "gene", drop = TRUE] %>% unique()
cluster_info <- sort(single.seurat$seurat_clusters)
exp_data <- GetAssayData(single.seurat, assay = "SCT", slot = "scale.data")
exp_data <- as.matrix(exp_data[intersect(topmarker, rownames(exp_data)), names(cluster_info)])
top_anno <- HeatmapAnnotation(cluster = anno_block(gp = gpar(fill = rainbow(length(levels(Idents(single.seurat)))), col = "white"),
                       labels = levels(cluster_info), 
                       labels_gp = gpar(cex = 0.4, col = "black"))) 
mark_gene <- diff.exp.all.filter %>% group_by(cluster) %>% top_n(., n = 2, wt = avg_logFC) %>% .[, "gene", drop = TRUE] %>% unique()
gene_pos <- match(mark_gene, rownames(exp_data))
grid.draw.Heatmap <- function(x){
  print(x)
}
row_anno <-  rowAnnotation(mark_gene = anno_mark(at = gene_pos, labels = mark_gene, labels_gp = gpar(fontsize = 8)))
P.heatmap <- Heatmap(exp_data,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_column_names = FALSE,
        show_row_names = FALSE,
        column_split = cluster_info,
        top_annotation = top_anno,
        right_annotation = row_anno,
        column_title = NULL,
		col = colorRamp2(c(-5, 1, 5), c("blue","white", "red")),
        heatmap_legend_param = list(
          title = ""
        ))
# plot dotplot					   
P.dotplot <- DotPlot(single.seurat, features = diff.exp.all %>% 
                       group_by(cluster) %>% 
                       top_n(., n = 2, wt = avg_logFC) %>%
                       .[, "gene", drop = TRUE] %>% unique(), dot.min = 0.3) + labs(x = "marker gene", y = "cluster") + theme(axis.text = element_text(size = 10)) + coord_flip()
SavePlot(od = outputdir, filename = "markergene_heatmap", data = P.heatmap)
SavePlot(od = outputdir, filename = "top2_markergene_dotplot", data = P.dotplot)

header <- c("ID", "symbol", "Pvalue", "log2FC", "pct.1", "pct.2", "Qvalue", "Cluster")
sapply(sort(unique(Idents(single.seurat))), function(x){
  clusterx <- diff.exp.all.filter %>% 
    group_by(cluster) %>% filter(., cluster == x)
  outfile <- paste(sample, paste("cluster", x, sep = ""),"diff_featuregene.xls", sep = ".")
  write.table(t(header), file = file.path(outputdir, outfile), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(clusterx, file = file.path(outputdir, outfile), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE, append = T)	
})
outputdir <- od
# write.table(t(header), file = file.path(outputdir, "nofilter_All_clusetr_featuregene.xls"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
# write.table(diff.exp.all, file = file.path(outputdir, "nofilter_All_clusetr_featuregene.xls"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE, append = T)
# write.table(t(header), file = file.path(outputdir, "filterd_All_clusetr_markergene.xls"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
# write.table(diff.exp.all.filter, file = file.path(outputdir, "filterd_All_clusetr_markergene.xls"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE, append = T)

# statistic marker gene in clusters
up.deg <- diff.exp.all.filter$cluster %>% table() %>% as.matrix() %>% t() %>% as.data.frame()
up.deg.stat <- data.frame(cluster = "deg_number", up.deg)
colnames(up.deg.stat) <- c("cluster", colnames(up.deg))
write.table(up.deg.stat, file = file.path(outputdir, paste(sample, "cluster_deg.stat.xls", sep = ".")), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)


# marker gene average expression in cluster
header <- c("ID", "symbol", paste("cluster", sort(unique(Idents(single.seurat))), sep = ""))
write.table(t(header), file = file.path(outputdir, paste(sample, "All_cluster_Markergene_avgExp.xls", sep = ".")), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(marker.gene.avgExp, file = file.path(outputdir, paste(sample, "All_cluster_Markergene_avgExp.xls", sep = ".")), sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE, append = T)

cluster.out <- data.frame(rownames(single.seurat@meta.data), single.seurat@meta.data$seurat_clusters)
names(cluster.out) <- c("Barcode", "Cluster")
write.table(cluster.out, file = file.path(outputdir, paste(sample, "clusters.xls", sep = ".")), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

counts <- single.seurat@assays$RNA@counts %>% as.data.frame %>% mutate(gene = rownames(.)) %>% AddId(., ensemble.id.symbol)
# scale.data <- single.seurat@assays$RNA@scale.data %>% as.data.frame %>% mutate(gene = rownames(.), .keep = "all") %>% AddId(., ensemble.id.symbol)
write.table(counts, file = file.path(outputdir, paste(sample, "All_cell_counts.xls", sep = ".")), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
# write.table(scale.data, file = file.path(outputdir, "All_cell_expression.tsv"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

save(single.seurat, file = file.path(outputdir, "single_analysis.Rda"))
saveRDS(single.seurat, file = file.path(outputdir, "single_seruat.Rds"))

time2 <- proc.time()
run.time <- time2 - time1
print(paste0('执行时间：',run.time[3][[1]],'秒'))
