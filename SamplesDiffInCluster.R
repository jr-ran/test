#!/usr/bin/Rscript

# parameter
suppressMessages(library(getopt))

spec <- matrix(c(
  'rda'        , 'R', 2, 'character', 'Rdata after step1 filter, example: single.Rda',
  'od'         , 'o', 2, 'character', 'out dir name',
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
if(is.null(opt$rda) && is.null(opt$od) && is.null(opt$id)) {PrintUsage()}
# if(is.null(opt$FDR)) {opt$FDR <- 0.05}
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
suppressMessages(library(parallel))
suppressMessages(library(circlize))
suppressMessages(library(ComplexHeatmap, lib.loc = "/share/nas1/ranjr/packages/3.6"))

time1 <- proc.time()
# parameter set
rda <- opt$rda
od <- opt$od
FDR <- opt$FDR
pvalue <- opt$pvalue
fold <- opt$fold
res <- opt$resolution
min.pct <- opt$min.pct
ensembleID_symbol <- opt$id
color <- opt$color
size <- opt$size
top.n <- opt$top.n

print(pvalue)
# function
SavePlot <- function(od, filename, data, width = 6, height = 4, scale = 1.3){
  file.png <- paste(filename, "png", sep = ".")
  file.pdf <- paste(filename, "pdf", sep = ".")
  ggsave(filename = file.path(od, file.png), plot = data, width = width, height = height, scale = scale)
  ggsave(filename = file.path(od, file.pdf), plot = data, width = width, height = height, scale = scale)
}

AddId <- function(data, id_list){
  data.re <- inner_join(data, id_list, by = c("gene" = "Symbol")) %>% 
    dplyr::select(EnsembleId, gene, dplyr::everything())
  return(data.re)
}

grid.draw.Heatmap <- function(x){
  print(x)
}


TopMarkerHeatmapAndDotplot <- function(seuratObject, diff.exp, topN, odir){
  topmarker <- diff.exp %>% group_by(cluster) %>% top_n(., n = top.n, wt = avg_logFC) 
  cluster_info <- sort(seuratObject$seurat_clusters)
  exp_data <- GetAssayData(seuratObject, assay = "SCT", slot = "scale.data")
  exp_data <- as.matrix(exp_data[intersect(topmarker$gene, rownames(exp_data)), names(cluster_info)])
  top_anno <- HeatmapAnnotation(cluster = anno_block(gp = gpar(fill = rainbow(length(levels(Idents(seuratObject)))), col = "white"),
                                                     labels = levels(cluster_info), 
                                                     labels_gp = gpar(cex = 0.8, col = "black"))) 
  mark_gene <- topmarker %>% group_by(cluster) %>% top_n(., n = 2, wt = avg_logFC) %>% .[, "gene", drop = TRUE] %>% unique()
  gene_pos <- match(mark_gene, rownames(exp_data))
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
  P.dotplot <- DotPlot(seuratObject, features = mark_gene, dot.min = 0.3) + labs(x = "marker gene", y = "cluster") + theme(axis.text = element_text(size = 10)) + coord_flip()
  SavePlot(od = odir, filename = paste("top", topN, "_", "markergene_heatmap", sep = ""), data = P.heatmap)
  SavePlot(od = odir, filename = "top2_markergene_dotplot", data = P.dotplot)
                       
}  
######################

load(rda)
ensemble.id.symbol <- read.table(file = ensembleID_symbol, header = FALSE, col.names = c("EnsembleId", "Symbol"))

# diff marker gene in clusters
print("diff marker gene in clusters start...")
logfc.threshold <- log2(fold)
diff.exp.all <- FindAllMarkers(object = single.integrated, assay = "SCT", slot = "counts", min.pct = 0.25, logfc.threshold = 0, verbose = FALSE)
diff.exp.all$avg_logFC <- diff.exp.all$avg_logFC/log(2) 
diff.exp.all <- AddId(diff.exp.all, ensemble.id.symbol)
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
print("diff marker gene in clusters done...")


# top10 marker gene plot
print("top 10 diff marker gene in clusters start...")
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

sapply(sort(unique(Idents(single.integrated))), function(x){
  
  diff.features <- diff.exp.all.filter[match(intersect(rownames(single.integrated), diff.exp.all.filter$gene), diff.exp.all.filter$gene),]
  top10.genes <- diff.features %>% 
    group_by(cluster) %>% filter(., cluster == x) %>%
    top_n(., n = top.n, wt = avg_logFC) %>%
    .[, "gene", drop = TRUE] %>% unique()
  len <- length(intersect(top10.genes, rownames(single.integrated)))
  print(len)
  #print(top10.genes)
  print("AAAAAAAAAAAAAAAAA")
  p.tsne <- FeaturePlot(single.integrated, features = top10.genes, slot = "scale.data", reduction = "tsne", cols = c("lightgrey", "purple"), ncol = 5, pt.size = 0.1) & NoLegend() & GlobPlotTheme()
  outfile <- paste(paste("cluster", x, sep = ""), paste("top", top.n, sep = ""), "markerTsne", sep = "_")
  SavePlot(od = outputdir, filename = outfile, data = p.tsne)
  p.umap <- FeaturePlot(single.integrated, features = top10.genes, slot = "scale.data", reduction = "umap", cols = c("lightgrey", "purple"), ncol = 5, pt.size = 0.1) & NoLegend() & GlobPlotTheme()
  outfile <- paste(paste("cluster", x, sep = ""), paste("top", top.n, sep = ""), "markerUmap", sep = "_")
  SavePlot(od = outputdir, filename = outfile, data = p.umap)
  p.vln <- VlnPlot(single.integrated, features = top10.genes, slot = "scale.data", ncol = 5, pt.size = 0) * theme(title = element_text(size = 7), 
      axis.ticks = element_blank(), 
      axis.title = element_blank(), 
      axis.text.y = element_text(size = 7),
	  axis.text.x  = element_text(size = 6, angle = 90))
  outfile <- paste(paste("cluster", x, sep = ""), paste("top", top.n, sep = ""), "markerVln", sep = "_")
  SavePlot(od = outputdir, filename = outfile, data = p.vln)
  finish <- paste("cluster", x, ":", top.n, "gene display is","done")
  return(finish)
})
print("top 10 diff marker gene in clusters done...")


# output cluster diff analysis file
outputdir <- paste(od, "statistic", sep = "/")
if(!dir.exists(outputdir)) {dir.create(outputdir)}

# marker gene expression
marker.gene.avgExp.data <- lapply(sort(unique(Idents(single.integrated))), function(x) {
  cells <- WhichCells(object = single.integrated, idents = x)
  t <- as.matrix(apply(single.integrated[['SCT']]@counts[match(diff.exp.all.filter$gene,rownames(single.integrated[['SCT']]@counts)), cells], 1, median), ncol = 1)
  colnames(t) <- x
  return(t)
})
marker.gene.avgExp <- do.call(cbind, marker.gene.avgExp.data)
marker.gene.avgExp %<>% as.data.frame %>% mutate(gene = rownames(.)) %>% AddId(., ensemble.id.symbol)
print("marker gene expression in clusters done...")
# marker gene average expression in cluster
header <- c("ID", "symbol", paste("cluster", sort(unique(Idents(single.integrated))), sep = ""))
write.table(t(header), file = file.path(outputdir, "All_cluster_Markergene_avgExp.xls"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(marker.gene.avgExp, file = file.path(outputdir, "All_cluster_Markergene_avgExp.xls"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE, append = T)

# plot heatmap
TopMarkerHeatmapAndDotplot(single.integrated, diff.exp.all.filter, top.n, outputdir)

header <- c("ID", "symbol", "Pvalue", "log2FC", "pct.1", "pct.2", "Qvalue", "Cluster")
sapply(sort(unique(Idents(single.integrated))), function(x){
  clusterx <- diff.exp.all.filter %>% 
    group_by(cluster) %>% filter(., cluster == x)
  outfile <- paste(paste("cluster", x, sep = ""),"diff_featuregene.xls", sep = ".")
  write.table(t(header), file = file.path(outputdir, outfile), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(clusterx, file = file.path(outputdir, outfile), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE, append = T)	
})
# statistic marker gene in clusters
# up.deg.stat <- diff.exp.all.filter$cluster %>% table() %>% as.matrix() %>% t() %>% as.data.frame() %>% mutate(cluster = "up_deg_number", .before = 1)
up.deg <- diff.exp.all.filter$cluster %>% table() %>% as.matrix() %>% t() %>% as.data.frame()
up.deg.stat <- data.frame(cluster = "deg_number", up.deg)
colnames(up.deg.stat) <- c("cluster", colnames(up.deg))
write.table(up.deg.stat, file = file.path(outputdir, "cluster_deg.stat.xls"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)

time2 <- proc.time()
run.time <- time2 - time1
print(paste0('执行时间：',run.time[3][[1]],'秒'))
