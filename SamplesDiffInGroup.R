#!/usr/bin/Rscript

# parameter
suppressMessages(library(getopt))

spec <- matrix(c(
  'rda'        , 'R', 2, 'character', 'Rdata after step1 filter, example: upload.Rdata',
  'group'      , 'g', 1, 'character', 'sample group for integrate',
  'od'         , 'o', 2, 'character', 'out dir name',
  'FDR'        , 'q', 1, 'numeric'  , 'FDR threshold',
  'pvalue'     , 'P', 1, 'numeric'  , 'pvalue',
  'fold'       , 'f', 1, 'numeric'  , 'fold change threshold',
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
if(is.null(opt$group)) {PrintUsage()}
# if(is.null(opt$FDR)) {opt$FDR <- 0.05}
# if(is.null(opt$fold)) {opt$fold <- 2}
if(is.null(opt$min.pct)) {opt$min.pct <- 0.25}
if(is.null(opt$color)) {opt$color <- 0.25}
if(is.null(opt$size)) {opt$size <- 10}
if(is.null(opt$top.n)) {opt$top.n <- 3}
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
rda <- opt$rda
group <- opt$group
od <- opt$od
FDR <- opt$FDR
pvalue <- opt$pvalue
fold <- opt$fold
min.pct <- opt$min.pct
ensembleID_symbol <- opt$id
color <- opt$color
size <- opt$size
top.n <- opt$top.n
logfc.threshold <- log2(fold)

samples.A <- strsplit(unlist(strsplit(group, split = "_vs_")), split = "-")[[1]]
samples.B <- strsplit(unlist(strsplit(group, split = "_vs_")), split = "-")[[2]]

# function
SavePlot <- function(od, filename, data, width = 6, height = 4, scale = 1.3){
  file.png <- paste(group, filename, "png", sep = ".")
  file.pdf <- paste(group, filename, "pdf", sep = ".")
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
  
GlobPlotTheme <- function (font.size = size){
  Theme <- theme(title = element_text(size = font.size), 
                   axis.ticks = element_blank(), 
                   axis.line = element_blank(), 
                   axis.title = element_blank(), 
                   axis.text = element_blank())
  return(Theme)
}

ExpHeatmapAndDotplot <- function(seurat.object, diff.exp, odir, cluster){
  # plot heatmap
  cluster_info <- WhichCells(object = seurat.object, idents = cluster)
  data <- seurat.object@meta.data[cluster_info,] %>% group_by(sample)
  sample_info = as.factor(data$sample)
  # exp_data <- GetAssayData(seurat.object, assay = "integrated", slot = "scale.data")
  exp_data <- GetAssayData(seurat.object, assay = "RNA", slot = "counts")
  exp_data <- as.matrix(exp_data[intersect(diff.exp$gene, rownames(exp_data)), data$cell])
  exp_data <- log10(exp_data+1)
  top_anno <- HeatmapAnnotation(
    cluster = anno_block(gp = gpar(fill = rainbow(length(levels(sample_info))), col = "white"), # 设置填充色
                         labels = levels(sample_info), 
                         labels_gp = gpar(cex = 0.7, col = "black"))) # 设置字体
  P.heatmap <- Heatmap(exp_data,
                       cluster_rows = FALSE,
                       cluster_columns = FALSE,
                       show_column_names = FALSE,
                       show_row_names = FALSE,
					   column_split = sample_info,
                       top_annotation = top_anno,
                       column_title = NULL,
					   col = colorRamp2(c(-3, 1, 3), c("blue","white", "red")),
                       heatmap_legend_param = list(
                         title = ""
                       ))
					  
  clusterName <- paste("cluster", cluster, sep = "")
  SavePlot(od = odir, filename = paste(clusterName, "markergene_heatmap", sep = "."), data = P.heatmap)
}



TopNPlot <- function(seurat.data, top.gene, odir, top = 3, x){
  n.col <- top
  print(top.gene)

  p.tsne <- FeaturePlot(seurat.data, features = top.gene, reduction = "tsne", cols = c("lightgrey", "red"), ncol = n.col, pt.size = 0.1) & NoLegend() & GlobPlotTheme()
  outfile <- paste(paste("cluster", x, sep = ""), paste("top", top, sep = ""), "markerTsne", sep = "_")
  SavePlot(od = odir, filename = outfile, data = p.tsne)
  p.umap <- FeaturePlot(single.integrated, features = top.gene, reduction = "umap", cols = c("lightgrey", "red"), ncol = n.col, pt.size = 0.1) & NoLegend() & GlobPlotTheme()
  outfile <- paste(paste("cluster", x, sep = ""), paste("top", top, sep = ""), "markerUmap", sep = "_")
  SavePlot(od = odir, filename = outfile, data = p.umap)
  p.vln <- VlnPlot(seurat.data, features = top.gene, ncol = n.col, pt.size = 0, group.by = "sample", idents = x) * theme(title = element_text(size = 7), 
                                                                                               axis.ticks = element_blank(), 
                                                                                               axis.title = element_blank(), 
                                                                                               axis.text.y = element_text(size = 7),
																							   axis.text.x  = element_text(size = 6, angle = 90))
  outfile <- paste(paste("cluster", x, sep = ""), paste("top", top, sep = ""), "markerVln", sep = "_")
  SavePlot(od = odir, filename = outfile, data = p.vln)
  
}

load(rda)
ensemble.id.symbol <- read.table(file = ensembleID_symbol, header = FALSE, col.names = c("EnsembleId", "Symbol"))

# diff analysis in samples
# samples.A <- c("T01", "T02")
# samples.B <- c("T03", "T04")
outputdir <- paste(od, "statistic", sep = "/")
if(!dir.exists(outputdir)) {dir.create(outputdir)}
odir <- paste(od, "topmarker", sep = "/")
if(!dir.exists(odir)) {dir.create(odir)}
print("diff analysis in samples start...")
diffList <- list()
for( x in 0:(length(unique(Idents(single.integrated)))-1)){
  print(paste("group diff analysis in", "cluster", x, "start..."))
  single.integrated@meta.data$cell <- rownames(single.integrated@meta.data)
  group.A <- single.integrated@meta.data %>% filter(seurat_clusters == x & sample == samples.A)
  group.B <- single.integrated@meta.data %>% filter(seurat_clusters == x & sample == samples.B)
  i <- x+1
  if (nrow(group.A) <= 2 | nrow(group.B) <= 2 ) {
    # return(NULL)
	diffList[i] <- NULL
	next
  }else{
	print("FindMarkers start...")
    # diff.all <- FindMarkers(object = single.integrated, assay = "SCT", slot = "counts", ident.1 = group.A$cell, ident.2 = group.B$cell, min.pct = 0.25, logfc.threshold = 0) %>% rownames_to_column("gene") 
	diff.all <- FindMarkers(object = single.integrated, assay = "RNA", slot = "data", ident.1 = group.A$cell, ident.2 = group.B$cell, min.pct = 0.25, logfc.threshold = 0) %>% rownames_to_column("gene")
	print("FindMarkers done...")
	print("diff start...")
	diff.all$avg_logFC <- diff.all$avg_logFC/log(2) 
	diff.all <- AddId(diff.all, ensemble.id.symbol)
	print("Add symbol done")
	if(!is.null(FDR)){
      if(!is.null(fold)){
        diff.filter <- diff.all %>% filter(p_val_adj < FDR & avg_logFC > logfc.threshold)
      }else{
        diff.filter <- diff.all %>% filter(p_val_adj < FDR)
      }
    }
    if(!is.null(pvalue)){
      if(!is.null(fold)){
        diff.filter <- diff.all %>% filter(p_val < pvalue & avg_logFC > logfc.threshold)
      }else{
        diff.filter <- diff.all %>% filter(p_val < pvalue)
      }
    }
	if (nrow(diff.filter) < 1){
	  # return(NULL)
	  diffList[i] <- NULL
	  next
	}
	print("diff done...")
	outfile <- paste(group, paste("cluster", x, sep = ""), "diff_featuregene.xls", sep = ".")
	header <- c("ID", "symbol", "Pvalue", "log2FC", "pct.1", "pct.2", "Qvalue")
	write.table(t(header), file = file.path(outputdir, outfile), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
	write.table(diff.filter, file = file.path(outputdir, outfile), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE, append = T)
	diff.features <- diff.filter[match(intersect(rownames(single.integrated), diff.filter$gene), diff.filter$gene),]
	if (nrow(diff.filter) < top.n){
	  top <- nrow(diff.filter)
	  topmarker <- diff.features %>% top_n(., n = top, wt = avg_logFC) %>% .[, "gene", drop = TRUE] %>% unique()
	  ExpHeatmapAndDotplot(single.integrated, diff.filter, outputdir, x)
	  TopNPlot(single.integrated, topmarker, odir, top, x)
	  
	}else{
	  top <- top.n
	  topmarker <- diff.features %>% top_n(., n = top, wt = avg_logFC) %>% .[, "gene", drop = TRUE] %>% unique()
	  ExpHeatmapAndDotplot(single.integrated, diff.filter, outputdir, x)
	  TopNPlot(single.integrated, topmarker, odir, top, x)
	}
	diffList[i] <- nrow(diff.filter)
  }
  print("diff analysis in samples done...")
}
print(diffList)
diffList <- as.list(diffList)
names(diffList) <- sort(unique(Idents(single.integrated)))
diffList <- diffList[which(sapply(diffList, is.null)==FALSE)]
deg.stat <- do.call(data.frame, diffList)
deg.stat <-  data.frame(cluster = group, diffList)
print(deg.stat)
colnames(deg.stat) <- c("cluster", names(diffList))
write.table(deg.stat, file = file.path(outputdir, "sample_deg.stat.xls"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
print("deg statistic in samples done...")

time2 <- proc.time()
run.time <- time2 - time1
print(paste0('执行时间：',run.time[3][[1]],'秒'))
