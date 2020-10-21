#!/usr/bin/Rscript

# parameter
suppressMessages(library(getopt))
spec <- matrix(c(
  'input.dirs'          , 'i', 2, 'character', 'input dir, Comma separated',
  'od'                  , 'o', 2, 'character', 'output dir name',
  'nUMI.min'            , 'u', 1, 'integer'  , 'min UMI number, default 100',
  'nGene.min'           , 'g', 1, 'integer'  , 'min gene number, default 500',
  'percent.mt.max'      , 'P', 1, 'numeric'  , 'mt percent, default 0.2',
  'min.cells'           , 'c', 1, 'integer'  , 'gene expression\'s min cells',
  'help'                , 'h', 0, 'logical'  , 'print usage'
), byrow = TRUE, ncol = 5)
opt <- getopt(spec)

PrintUsage <- function(){
  cat("ProgramName:filterAndMergeData.R
Version: Version v1.0
Contact: Ranjr <ranjr@biomarker.com.cn> 
Program Date:   2020.9.2 
Description: this program is used to filter single cell data ;\n")
  cat(paste(getopt(spec, usage = TRUE), sep = "\n"))
  q(status = 1)		
}

if(!is.null(opt$help)) {PrintUsage()}
if(is.null(opt$input.dirs) & is.null(opt$od)) {PrintUsage()}
if(is.null(opt$nUMI.min)) {opt$nUMI.min <- 100}
if(is.null(opt$nGene.min)) {opt$nGene.min <- 500}
if(is.null(opt$percent.mt.max)) {opt$percent.mt.max <- 0.2}
if(is.null(opt$min.cells)) {opt$min.cells <- 10}
if(!dir.exists(opt$od)) {dir.create(opt$od)}

# load packages
suppressMessages(library(SingleCellExperiment))
suppressMessages(library(Seurat))
suppressMessages(library(tidyverse))
suppressMessages(library(Matrix))
suppressMessages(library(scales))
suppressMessages(library(cowplot))
suppressMessages(library(RCurl))
suppressMessages(library(parallel))

# parameter set
nUMI.min <- opt$nUMI.min
nGene.min <- opt$nGene.min
pct.mt.max <- opt$percent.mt.max
cells.min <- opt$min.cells
od <- opt$od
object <- strsplit(opt$input.dirs, split = ",")[[1]]
sample <- basename(unlist(strsplit(object, split = "/outs"))[c(1,3)])

# function
SavePlot <- function(od, filename, data, width = 6, height = 4, scale = 1.3){
  file.png <- paste(filename, "png", sep = ".")
  file.pdf <- paste(filename, "pdf", sep = ".")
  ggsave(filename = file.path(od, file.png), plot = data, width = width, height = height, scale = scale)
  ggsave(filename = file.path(od, file.pdf), plot = data, width = width, height = height, scale = scale)
}

# parallel computing
outputdir <- paste(od, "singleSample", sep = "/")
if(!dir.exists(outputdir)) {dir.create(outputdir)}
print(length(object))
singleList <- parallel::mclapply(object, function(x) {
  single.matrix <- Read10X(data.dir = x)
  sample.name <- basename(strsplit(x, split = "/outs")[[1]][1])
  singledir <- paste(outputdir, sample.name, sep = "/")
  if(!dir.exists(singledir)) {dir.create(singledir)}
  seurat.object <- CreateSeuratObject(counts = single.matrix, project = sample.name, min.cells = cells.min)
  seurat.object[['sample']] <- sample.name
  seurat.object$"percent.mt" <- PercentageFeatureSet(object = seurat.object, pattern = "^MT-|^mt-")
  seurat.object$"percent.mt" <- seurat.object$"percent.mt" /100
  
  seurat.object@meta.data$nUMI <- seurat.object@meta.data$nCount_RNA 
  seurat.object@meta.data$nGene <- seurat.object@meta.data$nFeature_RNA 
  P.vlnplot.pre <- VlnPlot(seurat.object,
                    features = c("nGene", "nUMI", "percent.mt"), 
                    ncol = 3, pt.size = 0.1) * theme(axis.title.x = element_blank(), axis.text = element_text(size = 10))
  SavePlot(od = singledir, filename = paste(sample.name, "nofilter.vln", sep = "_"), data = P.vlnplot.pre)
  
  P1.scatter.pre <- FeatureScatter(seurat.object, feature1 = "nUMI", feature2 = "percent.mt", pt.size = 0.5) + theme(legend.position = "none", axis.text = element_text(size = 10))
  P2.scatter.pre <- FeatureScatter(seurat.object, feature1 = "nUMI", feature2 = "nGene", pt.size = 0.5) + theme(axis.text = element_text(size = 10))
  P.scatter.pre <- P1.scatter.pre + P2.scatter.pre
  SavePlot(od = singledir, filename = paste(sample.name, "nofilter.scatter", sep = "_"), data = P.scatter.pre)
  single.seurat <- subset(x = seurat.object, 
                        subset= ((nFeature_RNA >= nGene.min) & 
                                   percent.mt < pct.mt.max) & 
								   nCount_RNA >= nUMI.min)
	
  P.vlnplot <- VlnPlot(single.seurat,
                    features = c("nGene", "nUMI", "percent.mt"), 
                    ncol = 3, pt.size = 0.1) * theme(axis.title.x = element_blank(), axis.text = element_text(size = 10))
  SavePlot(od = singledir, filename = paste(sample.name, "vln", sep = "_"), data = P.vlnplot)

  P1.scatter <- FeatureScatter(single.seurat, feature1 = "nUMI", feature2 = "percent.mt", pt.size = 0.5) + theme(legend.position = "none", axis.text = element_text(size = 10))
  P2.scatter <- FeatureScatter(single.seurat, feature1 = "nUMI", feature2 = "nGene", pt.size = 0.5) + theme(axis.text = element_text(size = 10))
  P.scatter <- P1.scatter + P2.scatter
  SavePlot(od = singledir, filename = paste(sample.name, "scatter", sep = "_"), data = P.scatter)
  save(single.seurat, file = file.path(singledir, "single_qc.Rda"))
  saveRDS(single.seurat, file = file.path(singledir, "single_qc.Rds"))
  return(single.seurat)
}, mc.cores = length(object))

names(singleList) <- sample
save(singleList,file = file.path(od, 'upload.Rdata'))



if(length(object) > 1){
  # merge data
  outputdir <- paste(od, "mergedAllSample", sep = "/")
  if(!dir.exists(outputdir)) {dir.create(outputdir)}
  # Created seruat merge object with y
  merge.y <- function(object.list) {
    y = c(object.list[[2]])
    if (length(object.list) > 2L) {
      for (i in 3L:length(object.list)) {
        y <- c(y, object.list[[i]])
      }
    }
    return(y)
  }
  # begind to merge seurat object data, reuslts is mergedAll.seurat
  mergedAll.seurat <- merge(x = singleList[[1]], y = merge.y(singleList), add.cell.ids = sample)
  
  # mitochondrion 
  mergedAll.seurat$"percent.mt" <- PercentageFeatureSet(object = mergedAll.seurat, pattern = "^MT-|^mt-")
  mergedAll.seurat$"percent.mt" <- mergedAll.seurat$"percent.mt" /100
  # mergedAll.seurat@meta.data <- mergedAll.seurat@meta.data %>% dplyr::rename(nUMI = nCount_RNA, nGene = nFeature_RNA)
  
  # QC
  P.vlnplot.pre <- VlnPlot(mergedAll.seurat,
                           features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                           ncol = 3, pt.size = 0.1) * theme(axis.title.x = element_blank(), axis.text = element_text(size = 10))
  SavePlot(od = outputdir, filename = "nonfilter_cells_qc_vln", data = P.vlnplot.pre)
  
  P1.scatter.pre <- FeatureScatter(mergedAll.seurat, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size = 0.5) + theme(legend.position = "none", axis.text = element_text(size = 10)) + labs(x = "nUMI")
  P2.scatter.pre <- FeatureScatter(mergedAll.seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size = 0.5) + theme(axis.text = element_text(size = 10)) + labs(x = "uUMI", y = "nGene")
  P.scatter.pre <- P1.scatter.pre + P2.scatter.pre
  SavePlot(od = outputdir, filename = "nonfilter_cells_qc_scatter", data = P.scatter.pre)
  
  # filter QC
  mergedAll.seurat <- subset(x = mergedAll.seurat, 
                             subset= ((nFeature_RNA >= nGene.min) & 
                                        percent.mt < pct.mt.max) & 
                               nCount_RNA >= nUMI.min)
  
  P.vlnplot <- VlnPlot(mergedAll.seurat,
                       features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                       ncol = 3, pt.size = 0.1) * theme(axis.title.x = element_blank(), axis.text = element_text(size = 10)) 
  SavePlot(od = outputdir, filename = "filter_cells_qc_vln", data = P.vlnplot)
  
  P1.scatter <- FeatureScatter(mergedAll.seurat, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size = 0.5) + theme(legend.position = "none", axis.text = element_text(size = 10)) + labs(x = "nUMI")
  P2.scatter <- FeatureScatter(mergedAll.seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size = 0.5) + theme(axis.text = element_text(size = 10)) + labs(x = "uUMI", y = "nGene")
  P.scatter <- P1.scatter + P2.scatter
  SavePlot(od = outputdir, filename = "filter_cells_qc_scatter", data = P.scatter)
  
  save(mergedAll.seurat, file = file.path(outputdir, "AllSample_qc.Rda"))
  saveRDS(mergedAll.seurat, file = file.path(outputdir, "AllSample_qc.Rds"))
  
# integrate data 
  outputdir <- paste(od, "integratedAllSample", sep = "/")
  if(!dir.exists(outputdir)) {dir.create(outputdir)}
  for (i in 1:length(singleList)){
    singleList[[i]] <- SCTransform(object = singleList[[i]], do.scale = TRUE, verbose = FALSE, vars.to.regress = c("nCount_RNA"))
  }
  names(single.list) <- sample
  pancreas.features <- SelectIntegrationFeatures(object.list = singleList, nfeatures = 3000)
  pancreas.list <- PrepSCTIntegration(object.list = singleList, anchor.features = pancreas.features)
  print("Done")
  single.anchors <- FindIntegrationAnchors(object.list = pancreas.list, anchor.features = pancreas.features, dims = 1:30, normalization.method = "SCT")
  single.integrated <- IntegrateData(anchorset = single.anchors, dims = 1:30)
  print("data integrated done...")
  
  # VlnPlot scatter 
  outputdir <- paste(od, "base", sep = "/")
  if(!dir.exists(outputdir)) {dir.create(outputdir)}
  single.integrated@meta.data$nUMI <- single.integrated@meta.data$nCount_RNA 
  single.integrated@meta.data$nGene <- single.integrated@meta.data$nFeature_RNA 
  P.vlnplot <- VlnPlot(single.integrated,
                       features = c("nGene", "nUMI", "percent.mt"), 
                       ncol = 3, pt.size = 0.1) * theme(axis.title.x = element_blank(), axis.text = element_text(size = 10))
  SavePlot(od = outputdir, filename = "All_vln", data = P.vlnplot)
  
  P1.scatter <- FeatureScatter(single.integrated, feature1 = "nUMI", feature2 = "percent.mt", pt.size = 0.5) + theme(legend.position = "none", axis.text = element_text(size = 10))
  P2.scatter <- FeatureScatter(single.integrated, feature1 = "nUMI", feature2 = "nGene", pt.size = 0.5) + theme(axis.text = element_text(size = 10))
  P.scatter <- P1.scatter + P2.scatter
  SavePlot(od = outputdir, filename = "All_scatter", data = P.scatter)
  print("VlnPlot and scatter done...")
  
  single.integrated <- ScaleData(single.integrated, verbose = FALSE)
  single.integrated <- RunPCA(single.integrated, features = VariableFeatures(single.integrated))
  
  # PCA plot
  outputdir <- paste(od, "reduction", sep = "/")
  if(!dir.exists(outputdir)) {dir.create(outputdir)}
  single.integrated <- JackStraw(single.integrated, num.replicate = 100)
  single.integrated <- ScoreJackStraw(single.integrated, dims = 1:20)
  P.Jack <- JackStrawPlot(single.integrated, dims = 1:20)
  P.Elbow <- ElbowPlot(single.integrated, ndims = 30)
  PCA.heatmap <- DimHeatmap(single.integrated, dims = 1:9)
  P.pca <- DimPlot(single.integrated, reduction  = "pca")
  SavePlot(od = outputdir, filename = "JackStrawPlot", data = P.Jack)
  SavePlot(od = outputdir, filename = "ElbowPlot", data = P.Elbow)
  SavePlot(od = outputdir, filename = "PCs_heatmap", data = PCA.heatmap)
  SavePlot(od = outputdir, filename = "PCs_heatmap", data = P.pca)
  
  single.integrated <- RunUMAP(single.integrated, features = VariableFeatures(single.integrated))
  single.integrated <- RunTSNE(single.integrated, features = VariableFeatures(single.integrated))
  single.integrated <- FindNeighbors(single.integrated, dims = 1:20)
  single.integrated <- FindClusters(single.integrated, resolution = res) # resolution 
  single.integrated[['cellcluster']] <- Idents(single.integrated)[rownames(single.integrated@meta.data)]
  
  # UMAP/TSNE plot
  P.umap.all <- UMAPPlot(single.integrated, label = TRUE) 
  P.umap.split <- UMAPPlot(single.integrated, split.by = "sample", label = TRUE)
  P.umap.sample <- UMAPPlot(single.integrated, group.by = "sample", label = TRUE)
  P.tsne.all <- TSNEPlot(single.integrated, label = TRUE) 
  P.tsne.split <- TSNEPlot(single.integrated, split.by = "sample", label = TRUE) 
  P.tsne.sample <- TSNEPlot(single.integrated, group.by = "sample") 
  SavePlot(od = outputdir, filename = "umap_all", data = P.umap.all)
  SavePlot(od = outputdir, filename = "umap_splt", data = P.umap.split)
  SavePlot(od = outputdir, filename = "umap_sample", data = P.umap.sample)
  SavePlot(od = outputdir, filename = "tsne_all", data = P.tsne.all)
  SavePlot(od = outputdir, filename = "tsne_splt", data = P.tsne.split)
  SavePlot(od = outputdir, filename = "tsne_sample", data = P.tsne.sample)
  
  pca.projection <- single.integrated@reductions$pca@cell.embeddings
  cells <- rownames(pca.projection)
  pca.out <- data.frame(Cells = cells, pca.projection)
  write.table(pca.out, file = file.path(outputdir, "pca_components.xls"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
  
  tsne.projection <- single.integrated@reductions$tsne@cell.embeddings
  cells <- rownames(tsne.projection)
  tsne.out <- data.frame(Cells = cells, tsne.projection)
  write.table(tsne.out, file = file.path(outputdir, "tsne_components.xls"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
  
  umap.projection <- single.integrated@reductions$umap@cell.embeddings
  cells <- rownames(umap.projection)
  umap.out <- data.frame(Cells = cells, umap.projection)
  write.table(umap.out, file = file.path(outputdir, "umap_components.xls"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
  print("reductions done...")
  
  # output other file
  outputdir <- od
  cluster.out <- data.frame(rownames(single.integrated@meta.data), single.integrated@meta.data$seurat_clusters)
  names(cluster.out) <- c("Barcode", "Cluster")
  write.table(cluster.out, file = file.path(outputdir, "All_ncells_clusters.xls"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
  
  # clusters
  # res <- 0.2
  col.nam <- paste("integrated_snn_res", res, sep = ".")
  stat.table <- table(single.integrated@meta.data[, col.nam], single.integrated@meta.data$sample)
  sample.cluster <- matrix(stat.table, ncol = 2)
  sample.cluster.out <- data.frame(rownames(stat.table), sample.cluster)
  names(sample.cluster.out) <- c("cluster", colnames(stat.table))
  write.table(sample.cluster.out, file = file.path(outputdir, "All_ncells_clusters.stat.xls"), sep = "\t", row.names = F, quote = FALSE)
  
  P.bar <- ggplot(data = single.integrated@meta.data, aes(x = seurat_clusters, fill = sample)) + geom_bar(position = "dodge") + 
    scale_fill_manual(values = brewer.pal(8, "Set2")[1:length(unique(single.integrated@meta.data$sample))]) + 
    theme_bw() + theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=0.5,colour="black"))
  SavePlot(od = outputdir, filename = "sample_barplot", data = P.bar)
  
  save(single.integrated, file = file.path(outputdir, "single_analysis.Rda"))
  saveRDS(single.integrated, file = file.path(outputdir, "single_seruat.Rds"))
}





