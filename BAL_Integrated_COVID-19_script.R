setwd("Set Data Path here")

"Load all libraries for Analysis"

library(stringr)
library(Seurat)
library(ggplot2)
library(SingleR)
library(tidyverse)
library(gridExtra)
library(cowplot)
library(DoubletFinder)
library(plotly)
library(ggExtra)
library(kableExtra)
library(knitr)

"Set the Data path "
data_path <- "Input Data Path here"

mtgenes= as.character(unlist(read.table("MT_genes.txt"),
                             use.names = F))

dirs <- list.dirs(data_path, recursive = FALSE)
dirs <- dirs[grepl("GEX_count_Mmul10$", dirs)]
dirs

path <- paste0(dirs,"/outs/raw_feature_bc_matrix.h5")

createobj <- function(obj,file){
  s <- gsub("_GEX_count_Mmul10","",str_split_fixed(obj,"BAL/",2))[,2]
  obj <- CreateSeuratObject(counts=Read10X_h5(file),project=s)
}


list_obj <- mapply(dirs,FUN=createobj,path)

BAL_merge <- merge(list_obj[[1]], y=list_obj[2:length(list_obj)], 
                   add.cell.ids=c("RLF10_Baseline_Treated","Sample5215_Baseline_Untreated",
                                  "RQv9_Baseline_Untreated","Sample5215_Day4_Untreated",
                                  "RLF10_Day4_Treated","RQv9_Day4_Untreated",
                                  "RVf12_Baseline_Treated","RHz12_Baseline_Untreated",
                                  "RVf12_Day4_Treated","RHz12_Day4_Untreated",
                                  "RHz12_Necropsy_Untreated"))

head(BAL_merge@meta.data, 5)
tail(BAL_merge@meta.data, 5)
BAL_merge$Sample <- str_match(row.names(BAL_merge@meta.data),
                              "[[:alnum:]]+_[[:alnum:]]+_[[:alnum:]]+")


table(BAL_merge$Sample)
head(BAL_merge@meta.data, 13)
BAL_merge@meta.data <- BAL_merge@meta.data %>% separate(Sample,c("Sample","Day","Type"),"_",
                                                        remove = FALSE)

BAL_merge[["percent.hbb"]]  = PercentageFeatureSet(BAL_merge, pattern = "^HBB")
BAL_merge[["percent.rps"]]  = PercentageFeatureSet(BAL_merge, pattern = "^RPS")
BAL_merge[["percent.rpl"]]  = PercentageFeatureSet(BAL_merge, pattern = "^RPL")
BAL_merge[["percent.mt"]]  = PercentageFeatureSet(BAL_merge, features = mtgenes)
BAL_merge$log10GenesPerUMI <- log10(BAL_merge$nFeature_RNA) / log10(BAL_merge$nCount_RNA)
metadata_merged = BAL_merge@meta.data

filtered_merged_seurat <- subset(x = BAL_merge, 
                                 subset= (nFeature_RNA >= 500) & (nFeature_RNA <=3500) &
                                   (nCount_RNA >= 250) & (log10GenesPerUMI >= 0.8)
                                 & (percent.hbb < 10) & 
                                   (percent.mt < 10) & (percent.rps < 10) & (percent.rpl < 10))

filtered_merged_seurat$Sample_name = paste0(filtered_merged_seurat@meta.data$Sample,
                                            "_", filtered_merged_seurat@meta.data$Day,"_",
                                            filtered_merged_seurat$Type)

Idents(filtered_merged_seurat) <- filtered_merged_seurat$Sample_name

#Testing ggscatterplot with density plot




pdf("BALs_Merged_VlnPlot.pdf", width =40, height =5)
VlnPlot(object = filtered_merged_seurat, features = c("nFeature_RNA", "nCount_RNA",
                                                      "log10GenesPerUMI","percent.hbb",
                                                      "percent.mt", "percent.rps",
                                                      "percent.rpl","precent.HTO"), ncol = 7)
dev.off()

##----------Gene Level Filtering------------##

"Output a logical vector for every gene on whether
          the more than zero counts per cell Extract counts"

counts <- GetAssayData(object = filtered_merged_seurat, slot = "counts")

"Output a logical vector for 
          every gene on whether the more than zero counts per cell"

nonzero <- counts > 0

"Sums all TRUE values and 
          returns TRUE if more than 10 TRUE values per gene"

keep_genes <- Matrix::rowSums(nonzero) >= 10

"Only keeping those genes 
          expressed in more than 10 cells"

filtered_counts <- counts[keep_genes, ]

"Reassign to filtered Seurat object"

filtered_merged_seurat <- CreateSeuratObject(filtered_counts, 
                                             meta.data = filtered_merged_seurat@meta.data)

metadata_clean_merged = filtered_merged_seurat@meta.data

##----------End of this Block-------------##


##----------VISUALIZATION-----------------##

" Visualize the number UMIs/transcripts per cell (Unfiltered)"

a = metadata_merged %>%
  ggplot(aes(color=orig.ident, x=nFeature_RNA, fill= orig.ident)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 150)+
  ggtitle("Unfiltered")

"Visualize the number UMIs/transcripts per cell (Filtered)"

b = metadata_clean_merged %>%
  ggplot(aes(color=orig.ident, x=nFeature_RNA, fill= orig.ident)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)+
  ggtitle("Filtered")
"Visualize the distribution of genes detected per cell via histogram (Unfiltered)"
c = metadata_merged %>%
  ggplot(aes(color=orig.ident, x=nCount_RNA, fill= orig.ident)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  scale_x_log10() +
  geom_vline(xintercept = 200)+
  ggtitle("Unfiltered")

"Visualize the distribution of genes detected per cell via histogram (filtered)"
d = metadata_clean_merged %>%
  ggplot(aes(color=orig.ident, x=nCount_RNA, fill= orig.ident)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  scale_x_log10() +
  geom_vline(xintercept = 700)+
  ggtitle("Filtered")

"Visualize the distribution of genes detected per cell via boxplot (Unfiltered)"
e = metadata_merged %>%
  ggplot(aes(x=orig.ident, y=log10(nCount_RNA), fill=orig.ident)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells vs NGenes (Unfiltered)")

"Visualize the distribution of genes detected per cell via boxplot (filtered)"
f = metadata_clean_merged %>%
  ggplot(aes(x=orig.ident, y=log10(nCount_RNA), fill=orig.ident)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells vs NGenes (Filtered)")


"Visualize the overall complexity of the gene expression
by visualizing the genes detected per UMI (Unfiltered)"
g = metadata_merged %>%
  ggplot(aes(x=log10GenesPerUMI, color = orig.ident, fill=orig.ident)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)+
  ggtitle("Unfiltered")

"Visualize the overall complexity of the
gene expression by visualizing the genes detected per UMI (filtered)"

h = metadata_clean_merged %>%
  ggplot(aes(x=log10GenesPerUMI, color = orig.ident, fill=orig.ident)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)+
  ggtitle("Filtered")


"Visualize the percent hbb (unfiltered)"
i = metadata_merged %>%
  ggplot(aes(x=percent.hbb, color = orig.ident, fill=orig.ident)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.1)+
  ggtitle("Unfiltered")

"Visualize the percent hbb (filtered)"
j = metadata_clean_merged %>%
  ggplot(aes(x=percent.hbb, color = orig.ident, fill=orig.ident)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.1)+
  ggtitle("Filtered")

"Visualize the percent mt (unfiltered)"
k = metadata_merged %>%
  ggplot(aes(x=percent.mt, color = orig.ident, fill=orig.ident)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.1)+
  ggtitle("Unfiltered")

"Visualize the percent mt (unfiltered)"
l = metadata_clean_merged %>%
  ggplot(aes(x=percent.mt, color = orig.ident, fill=orig.ident)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.1)+
  ggtitle("Filtered")

"Visualize the percent rps (unfiltered)"
m = metadata_merged %>%
  ggplot(aes(x=percent.rps, color = orig.ident, fill=orig.ident)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.1)+
  ggtitle("Unfiltered")

" Visualize the percent rps (filtered)"
n = metadata_clean_merged %>%
  ggplot(aes(x=percent.rps, color = orig.ident, fill=orig.ident)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.1)+
  ggtitle("Filtered")

"Visualize the percent rpl (unfiltered)"
o = metadata_merged %>%
  ggplot(aes(x=percent.rpl, color = orig.ident, fill=orig.ident)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.1)+
  ggtitle("Unfiltered")

"Visualize the percent rpl (filtered)"

p = metadata_clean_merged %>%
  ggplot(aes(x=percent.rpl, color = orig.ident, fill=orig.ident)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.1)+
  ggtitle("Filtered")

pdf(paste0("/Merged_BALs","_","InitialQC_Unfiltered_vs_Filtered.pdf"),
    width = 20, height = 30)
plot_grid(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p, nrow = 8 )
dev.off()

##----------End of this Block-------------##


##----------Perform Doublet detection and removal-------------##

"Split By Sample"

BAL.list_0 <- SplitObject(filtered_merged_seurat, split.by = "Sample_name")

for (i in names(BAL.list_0)) {
#normalization
  BAL.list_0[[i]] <- NormalizeData(BAL.list_0[[i]])
  BAL.list_0[[i]] <- FindVariableFeatures(BAL.list_0[[i]], selection.method = "vst",
                                              nfeatures = 3500)
  BAL.list_0[[i]] <- ScaleData(BAL.list_0[[i]])
  BAL.list_0[[i]] <- RunPCA(BAL.list_0[[i]])
  BAL.list_0[[i]] <- RunUMAP(BAL.list_0[[i]], dims = 1:30)

  sweep.res.list_Bal <- paramSweep_v3(BAL.list_0[[i]], PCs = 1:30, sct = F)
  #gt.calls = sweep_RLF10@meta.data[rownames(sweep_RLF10[[1]]), "GT"]
  sweep.stats_bal <- summarizeSweep(sweep.res.list_Bal, GT = F)
  bcmvn_bal <- find.pK(sweep.stats_bal)

  BAL.list_0[[i]] = FindNeighbors(BAL.list_0[[i]], reduction = "pca")
  BAL.list_0[[i]] = FindClusters(BAL.list_0[[i]], resolution = 0.6)

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
  annotations = BAL.list_0[[i]]$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)          
  nExp_poi <- round(0.075*length(colnames(BAL.list_0[[i]])))  
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
  BAL.list_0[[i]] <- doubletFinder_v3(BAL.list_0[[i]], PCs = 1:30, pN = 0.25, pK = 0.09,
                               nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  
  BAL.list_0[[i]]$pANN = BAL.list_0[[i]]@meta.data[[15]]

  BAL.list_0[[i]] <- doubletFinder_v3(BAL.list_0[[i]], PCs = 1:30, pN = 0.25, 
                                          pK = 0.09, nExp = nExp_poi.adj, 
                               reuse.pANN = "pANN",
                               sct = FALSE)
 
##Save all the Doublet-singlet Plots
  Idents(BAL.list_0[[i]]) = BAL.list_0[[i]]@meta.data[[16]]
  pdf(paste0(BAL.list_0[[i]]$Sample_name,"_", "Doublets_dimplot.pdf"))
  DimPlot(BAL.list_0[[i]])
  dev.off()
  
}

for (i in names(BAL.list_0)){
  BAL.list_0[[i]] = SubsetData(BAL.list_0[[i]], cells = rownames(BAL.list_0[[i]]@meta.data)[ which(BAL.list_0[[i]]@meta.data[[16]] == "Singlet") ])
  
}


"Perform SCTransform on 
          all samples in the list"

for (i in names(BAL.list_0)) {
  BAL.list_0[[i]] <- SCTransform(BAL.list_0[[i]], verbose = FALSE, 
                                  vars.to.regress = c("percent.hbb", 
                                                      "percent.mt",
                                                      "percent.rps", 
                                                      "percent.rpl"))
}

BAL_doubletremoved_sctransformed = BAL.list_0

"Select Integration features"

BAL.features <- SelectIntegrationFeatures(object.list = BAL_doubletremoved_sctransformed, nfeatures = 3500)

"Prepare 'SCT' based Integration"

BAL_doubletremoved_sctransformed <- PrepSCTIntegration(object.list = BAL_doubletremoved_sctransformed, anchor.features = BAL.features)

"Find Anchors within datasets(samples) 
      to integrate using 'SCT'"

BAL.anchors <- FindIntegrationAnchors(object.list = BAL_doubletremoved_sctransformed, normalization.method = "SCT", 
                                      anchor.features = BAL.features)

BAL.integrated <- IntegrateData(anchorset = BAL.anchors, normalization.method = "SCT")

"Perform PCA and UMAP"

BAL.integrated <- RunPCA(object = BAL.integrated, verbose = FALSE)
BAL.integrated <- RunUMAP(object = BAL.integrated, dims = 1:30)
BAL.integrated = FindNeighbors(BAL.integrated, reduction = "pca", dims = 1:13)
BAL.integrated <- FindClusters(object = BAL.integrated, graph.name = "integrated_snn" ,resolution = 0.6)
DimPlot(BAL.integrated, split.by = "seurat_clusters", ncol = 5)
DimPlot(BAL.integrated, split.by = "Type")


"Dimplot UMAP plotted by Groups- Sample and Day"

plots <- DimPlot(BAL.integrated, group.by = c("Sample", "Day"))
pdf("/Integrated_BALs_Split_by_Sample&Day_UMAP.pdf", width = 15, height = 10)
plots & theme(legend.position = "top") & 
  guides(color = guide_legend(nrow = 4, byrow = TRUE,override.aes = list(size = 2.5)))
dev.off()

"Dimplot of UMAP Split by Day"

pdf("/Integrated_BALs_Split_by_Sample&Day_UMAP.pdf", width = 15, height = 10)
DimPlot(BAL.integrated, reduction = "umap", split.by = "Day")
dev.off()

"Dimplot with Sample Name"

Idents(BAL.integrated) = BAL.integrated$Sample_name
pdf("/Integrated_BALs_Split_by_SampleName_UMAP.pdf", width = 15, height = 10)
DimPlot(BAL.integrated, reduction = "umap")
dev.off()

##----------End of this Block-------------##


##--------------Panel of Gene Markers--------------##

"Assign RNA as Default Assay"

DefaultAssay(BAL.integrated) <- "RNA"

"Normalize RNA data for 
          visualization purposes"

BAL.integrated <- NormalizeData(BAL.integrated, verbose = FALSE)

"Get a list of Featues (Genes) 
        to plot"

features = c("MS4A1","FOXP3","GATA3", "PTPRC","RORC","TBX21","BCL6", "CCR7","CD28","CD4", "CD69",
             "CD3D","CD3E", "FASLG","IL2RA", "CD8A", "GZMB","FCER1G","ITGAM","CD163",
             "CCR2","CX3CR1", "ITGAX", "FCGR3", "KLRG1", "AXL", "SIGLEC6","TNFRSF17",
             "CD101","CSF1R", "CXCR5","CXCR6","IRF4", "ITGA1","ITGAE","PDCD1",
             "PRDM1","SPI1","STAT3","STAT4","STAT5A", "TCF7", "TOX", "ISG15",
             "IFI6","CXCL8","CXCL3","CCL2")

"Plot FeaturePlot (UMAP) 
      for Immune Cell Markers"

p <- FeaturePlot(BAL.integrated, features=features, ncol=4)
pdf(paste0("/Integrated_BALs","_FeaturePlot_MarkerGenes_umap.pdf"), 
    height = 30, width = 15)
p
dev.off()

"Assign 'Sample' Column as Idents & 
          Plot ViolinPlots for Immune Cell Markers by Sample"

Idents(BAL.integrated) = BAL.integrated$Sample
p <- VlnPlot(BAL.integrated, features=features, pt.size = 0)
pdf(paste0("/Integrated_BALs","_VLNplot_MarkerGenes_bySample.pdf"), 
    height = 20, width = 30)
p
dev.off()

"Plot RidgePlots for 
      Immune cell markers group by Sample"

p = RidgePlot(BAL.integrated, features = features, group.by = "Sample")
pdf(paste0("/Integrated_BALs","_RidgePlots_MarkerGenes_bySample.pdf"),
    height = 20, width = 30)
p
dev.off()

"Plot RidgePlots for 
      Immune cell markers group by Day"

p = RidgePlot(BAL.integrated, features = features, group.by = "Day")
pdf(paste0("/Integrated_BALs","_RidgePlots_MarkerGenes_byDay.pdf"),
    height = 20, width = 30)
p
dev.off()

"Plot DimPlot for 
      Immune cell markers split by Day"

Idents(BAL.integrated) = BAL.integrated$Day
p = DimPlot(BAL.integrated,split.by  = "Day")
png(paste0("/Integrated_BALs","_DimPlots_MarkerGenes_byDay.pdf"),
    height = 10, width = 15)
p
dev.off()

"Plot Dimplot fo Immune cell markers Split
      by Type_Day"
BAL.integrated$Type_Day = paste0(BAL.integrated$Type, "_", BAL.integrated$Day)
Idents(BAL.integrated) = BAL.integrated$Type_Day
p = DimPlot(BAL.integrated,split.by  = "Type_Day", ncol = 4)
pdf(paste0("/Integrated_BALs","_DimPlots_MarkerGenes_byType_Day.pdf"),
    height = 10, width = 15)
p
dev.off()

##----------End of this Block-------------##

saveRDS(BAL.integrated, "/BAL_Integrated_Seurat_Obj.rds")

##-----------Run SingleR for BP Encode------------##

"Download BP encode
          database"

BPencode.se = readRDS("Analysis_AB/BPEncode.rds")
#hpce.se = readRDS("Analysis_AB/HumanEncode.rds")

BAL_Integrated.sce <- as.SingleCellExperiment(BAL.integrated)

"Convert the seurat object to 
          singlecellexperiment"

singleR_BAL_intergrated_bp = SingleR(test = BAL_Integrated.sce, ref = BPencode.se,
                                     labels = BPencode.se$label.main, fine.tune = T,
                                     prune = T, BPPARAM = MulticoreParam(4))

"Map co-ordinates of UMAP onto 
      singleR object"

BAL.integrated$singleRclusters_BP_pruned = singleR_BAL_intergrated_bp$pruned.labels

"Plot DimPlot for singleR clusters-BP
      split by singleRclusters"

Idents(BAL.integrated) = BAL.integrated$singleRclusters_BP_pruned
p = DimPlot(BAL.integrated,split.by  = "singleRclusters_BP_pruned", ncol =4)
pdf(paste0("/Integrated_BALs","_DimPlots_pruned-singleR.pdf"),
    height = 10, width = 15)
p
dev.off()

"Plot Dimplot for 
      SingleR cluster-BP"

p = DimPlot(BAL.integrated, ncol =4)
pdf(paste0("/Integrated_BALs","_DimPlots_pruned-singleR.pdf"),
    height = 10, width = 15)
p
dev.off()

"Plot Violin plots split 
      by singleR clusters-BP"

p <- VlnPlot(BAL.integrated, features=features, ncol = 5, pt.size = 0)
pdf(paste0("/Integrated_BALs","_VLNplot_MarkerGenes_bySingleRClusters-BP.pdf"),
    height = 25, width = 30)
p
dev.off()

"Plot DotPlots across 
      singleR clusters-BP encode"

p <- DotPlot(BAL.integrated, features=features, cols = c("Blue", "Red"),
             col.min = -1, col.max = 1, dot.scale = 10) + RotatedAxis()
pdf(paste0("/Integrated_BALs","_DotPlot_MarkerGenes_bySingleRClusters-BP.pdf"), 
    height = 20, width = 30)
p
dev.off()

"Scale data for sake of plotting heatmap"

data_Integrated_BAL <- ScaleData(object = BAL.integrated, features = rownames(BAL.integrated))

"Plot Heatmap split 
      by singleR clusters-BP"

p <- DoHeatmap(subset(BAL.integrated,downsample =50), features = features, 
               size = 3,  slot = "scale.data", angle = 45, raster = T, draw.lines = T,
               group.bar.height = 0.01, disp.min = -1, disp.max = 1, combine = T)+scale_fill_gradientn(
                 colors = rev(RColorBrewer::brewer.pal(n = 8,name = "RdBu")) ) + guides(color=FALSE)
pdf(paste0("/Integrated_BALs","_Heatmap_bySingleRClusters-BP.pdf"),
    height = 10, width = 20)
p
dev.off()


"Plot Heatmap split 
      by Type"

Idents(BAL.integrated) = "Type"
p <- DoHeatmap(subset(BAL.integrated,downsample =50), features = features, 
               size = 3,  slot = "scale.data", angle = 45, raster = T,
               group.bar.height = 0.01, disp.min = -1, disp.max = 1, combine = T)+scale_fill_gradientn(
                 colors = rev(RColorBrewer::brewer.pal(n = 8,name = "RdBu")) ) + guides(color=FALSE)
pdf(paste0("/Integrated_BALs","_Heatmap_byType.pdf"),
    height = 10, width = 20)
p
dev.off()

"Plot Heatmap split 
      by Type_Day"

Idents(BAL.integrated) = "Type_Day"
p <- DoHeatmap(subset(BAL.integrated,downsample =50), features = features, 
               size = 3,  slot = "scale.data", angle = 45, raster = T,
               group.bar.height = 0.01, disp.min = -1, disp.max = 1, combine = T)+scale_fill_gradientn(
                 colors = rev(RColorBrewer::brewer.pal(n = 8,name = "RdBu")) ) + guides(color=FALSE)
pdf(paste0("/Integrated_BALs","_Heatmap_byType_Day.pdf"),
    height = 10, width = 20)
p
dev.off()


##----------End of this Block-------------##


##---------Get Cell Counts and Proportions-----------##

"Save the file to an RDS object"
saveRDS(BAL.integrated, file = "/BAL_Integrated_SingleR_SeuratObj.rds")


#Get proportions of cells per Single R cluster

x = knitr::kable(prop.table(table(BAL.integrated$Celltype)),"html")
kable_styling(x, font_size = 50, bootstrap_options = "striped", full_width = F) %>% save_kable("BAL_Integrated/BAL_Integrated_Cell_Cluster_Proportions.png")

e#Get Count of cells per Single R cluster
x = knitr::kable(table(BAL.integrated$Celltype),"html")
kable_styling(x, font_size = 45, bootstrap_options = "striped", full_width = F) %>% save_kable("BAL_Integrated/BAL_Integrated_Cell_Cluster_counts.png")



##-----Make 3D umap------##


# Re-run UMAPs that you have accurate calculations for all UMAP(s)
BAL.integrated <- RunUMAP(BAL.integrated,
                            dims = 1:13,
                            n.components = 3L)

# Extract tSNE information from Seurat Object
umap_1 <- BAL.integrated[["umap"]]@cell.embeddings[,1]
umap_2 <- BAL.integrated[["umap"]]@cell.embeddings[,2]
umap_3 <- BAL.integrated[["umap"]]@cell.embeddings[,3]

# Visualize what headings are called so that you can extract them to form a dataframe
Embeddings(object = BAL.integrated, reduction = "umap")

# Prepare a dataframe for cell plotting
plot.data <- FetchData(object = BAL.integrated, vars = c("UMAP_1", "UMAP_2", "UMAP_3", "singleRclusters_BP_pruned"))

# Make a column of row name identities (these will be your cell/barcode names)
plot.data$label <- paste(rownames(plot.data))

# Plot your data, in this example my Seurat object had 21 clusters (0-20)
plot_ly(data = plot.data, 
        x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, 
        color = ~singleRclusters_BP_pruned, 
        colors = c("lightseagreen",
                   "gray50",
                   "darkgreen",
                   "red4",
                   "red",
                   "turquoise4",
                   "black",
                   "yellow4",
                   "royalblue1",
                   "lightcyan3",
                   "peachpuff3",
                   "khaki3",
                   "gray20",
                   "orange2",
                   "royalblue4",
                   "yellow3",
                   "gray80",
                   "darkorchid1",
                   "lawngreen",
                   "darkmagenta"),
        type = "scatter3d", 
        mode = "markers", 
        marker = list(size = 2, width=2), # controls size of points
        text=~label, #This is that extra column we made earlier for which we will use for cell ID
        hoverinfo="text")


"Manual Cell selection"
BAL.integrated@meta.data$Celltype = "NA"
Idents(BAL.integrated) =BAL.integrated$singleRclusters_BP_pruned
plot = DimPlot(BAL.integrated)
select.cells = CellSelector(plot = plot)

BAL.integrated@meta.data[select.cells,]$Celltype = "C1"

select.cells = CellSelector(plot = plot)

BAL.integrated@meta.data[select.cells,]$Celltype = "C2"

select.cells = CellSelector(plot = plot)

BAL.integrated@meta.data[select.cells,]$Celltype = "C3"

select.cells = CellSelector(plot = plot)

BAL.integrated@meta.data[select.cells,]$Celltype = "C4"

select.cells = CellSelector(plot = plot)

BAL.integrated@meta.data[select.cells,]$Celltype = "Macrophages/Monocytes"

select.cells = CellSelector(plot = plot)

BAL.integrated@meta.data[select.cells,]$Celltype = "C6"

select.cells = CellSelector(plot = plot)

BAL.integrated@meta.data[select.cells,]$Celltype = "C7"
#Plot UMAP with new cluster types

Idents(BAL.integrated) = BAL.integrated$Celltype
DimPlot(BAL.integrated)

#Assign cell types 

BAL.integrated$Celltype[BAL.integrated$Celltype == "C1"] <- "Macrophages"
BAL.integrated$Celltype[BAL.integrated$Celltype == "C2"] <- "DC" # Check
BAL.integrated$Celltype[BAL.integrated$Celltype == "C3"] <- "Unclassified"
BAL.integrated$Celltype[BAL.integrated$Celltype == "C4"] <- "HSCs"
BAL.integrated$Celltype[BAL.integrated$Celltype == "C5"] <- "T-cells/NKcells"
BAL.integrated$Celltype[BAL.integrated$Celltype == "C6"] <- "Epithelial Cells" # Check
BAL.integrated$Celltype[BAL.integrated$Celltype == "NA"] <- "Unclassified"
Idents(BAL.integrated) = BAL.integrated$Celltype



Idents(BAL.integrated) = "Celltype"
pdf(paste0("BAL_Integrated_UMAP_splitBy_TypeDay.pdf"), 
    height = 10, width = 15)
DimPlot(BAL.integrated, split.by = "Type_Day", ncol = 2,
        cols = c("steelblue2","gray","springgreen3","tan3","royalblue1",
                 "yellow3",
                 "red3","mediumturquoise"))
dev.off()

saveRDS(BAL.integrated, file = "/BAL_Integrated_New_Celltypes_object.rds")

"Read SingleRSeurat-object to a file"


BAL.integrated = readRDS(file = "BAL_Integrated/BAL_Integrated_New_Celltypes_object.rds")

Idents(BAL.integrated) = "Type_Day"

BAL.integrated = subset(BAL.integrated, subset =  (Type_Day != "Untreated_Necropsy"))

Idents(BAL.integrated) = "Celltype"
png("Panel_A_UMAP.png", height = 700, width = 1000)
DimPlot(BAL.integrated,cols = c("steelblue2","gray","springgreen3","tan3","royalblue1",
                                "yellow3",
                                "red3","mediumturquoise"), dims = c(1,2))&
  ggplot2::theme(axis.title = element_blank(), 
                 axis.text = element_blank(),
                 axis.ticks = element_blank(),
                 strip.text = element_blank(), 
                 legend.position = "right", legend.direction="vertical",legend.key.size = unit(0.5, "in"),
                 axis.line   = element_blank(), legend.text=element_text(size=25) ) & labs(title = "")
dev.off()

#Plot ISGs and Inflammatory Genes

ISGs = c("IFI6","IFI27","ISG15","MX1","ISG20","MX2","OAS2","IFIT3",
         "IRF7","OAS1")

CCK = c("CCL4L1", "CXCL10", "CXCL3","CXCL8")

IFGs = c("IL6", "TNF", 
             "IL10", "IL1B", "IFNB1")

Idents(BAL.integrated) = "Celltype"
BAL.integrated_scaled= ScaleData(BAL.integrated)
DefaultAssay(BAL.integrated_scaled) = "RNA"
Idents(BAL.integrated_scaled) = "Celltype"

p <- DoHeatmap(subset(BAL.integrated_scaled, downsample = 50), features = ISGs, 
               size = 3,  slot = "scale.data", angle = 45, raster = T,
               group.bar.height = 0.01, disp.min = -1, disp.max = 1, combine = T)+scale_fill_gradientn(
                 colors = rev(RColorBrewer::brewer.pal(n = 8,name = "RdBu"))) + guides(color=FALSE)
pdf(paste("BAL_Integrated_DimHeatmap_",st, ".pdf", sep = ""), height = 10, width = 10)
p
dev.off()

#Plotting DotPlots-InflammatoryGenes
BAL.integrated = subset(BAL.integrated, subset = (Type_Day != "Untreated_Necropsy"))
Subsetted_Macrophages = subset(BAL.integrated, subset = (Celltype == "Macrophages"))
Subsetted_Macrophages$CellType_Typeday = paste0(Subsetted_Macrophages$Celltype,"_", Subsetted_Macrophages$Type_Day)
Subsetted_Macrophages$Day_Type = paste0(Subsetted_Macrophages$Day,"_", Subsetted_Macrophages$Type)

Idents(Subsetted_Macrophages) = "Day_Type"
pdf(paste0("Macrophages_splitBy_Day_Type.pdf"))
Baseline_Treated <- WhichCells(Subsetted_Macrophages, idents = c("Baseline_Treated"))
Baseline_Untreated <- WhichCells(Subsetted_Macrophages, idents = c( "Baseline_Untreated"))
Day4_Treated <- WhichCells(Subsetted_Macrophages, idents = c( "Day4_Treated"))
Day4_Untreated <- WhichCells(Subsetted_Macrophages, idents = c( "Day4_Untreated"))


BAL.integrated$Day_Type = paste0(BAL.integrated$Day,"_", BAL.integrated$Type)
Idents(BAL.integrated) = "Day_Type"

p = DimPlot(BAL.integrated, split.by ="Day_Type",  cols = c("royalblue", 
                                                        "indianred2", 
                                                        "forestgreen","goldenrod2"), ncol = 2) &
  ggplot2::theme(axis.title = element_blank(), 
                 axis.text = element_blank(),
                 axis.ticks = element_blank(),
                 strip.text = element_blank(), 
                 legend.position = "right", legend.direction="vertical",legend.key.size = unit(0.5, "in"),
                 axis.line   = element_blank(), legend.text=element_text(size=25) ) & labs(title = "")

png(paste0("Integrated_splitBy_DayType.png"),
    width = 900, height = 600)
p
dev.off()



P = DimPlot(Subsetted_Macrophages, split.by ="Day_Type",  cols = c("royalblue", 
                                                            "indianred2", 
                                                            "forestgreen","goldenrod2"), ncol = 2) &
ggplot2::theme(axis.title = element_blank(), 
                 axis.text = element_blank(),
                 axis.ticks = element_blank(),
                 strip.text = element_blank(), 
                 legend.position = "bottom", legend.direction="horizontal",legend.key.size = unit(0.35, "in"),
                 axis.line   = element_blank() ) & labs(title = "")

png(paste0("Subsetted_Macrophages_splitBy_DayType.png"),
    width = 700, height = 600)
P
dev.off()

#Dotplots for IFGs
Idents(BAL.integrated) = "Day_Type"
p = DotPlot(BAL.integrated, features=IFGs, cols = c("navyblue","red"),
            col.min = -1,group.by  = "Day_Type",col.max = 1,  dot.scale = 10, scale.by = "size") + RotatedAxis()


pdf(paste0("Macrophages-IFGs_DotPlots_byTypeDay.pdf"), height = 5, width = 10)
p
dev.off()

#Dotplots for ISGs
Idents(BAL.integrated) = "Day_Type"
p = DotPlot(BAL.integrated, features=ISGs, cols = c("navyblue","red"),
            col.min = -1,group.by  = "Day_Type",col.max = 1,  dot.scale = 10, scale.by = "size") + RotatedAxis()


pdf(paste0("Macrophages-ISGs_DotPlots_byTypeDay.pdf"), height = 7, width = 10)
p
dev.off()

#Dotplots for CCKs
Idents(Subsetted_Macrophages) = "Day_Type"
p = DotPlot(Subsetted_Macrophages, features=CCK, cols = c("navyblue","red"),
            group.by = "Day_Type",
            col.min = -1, col.max = 1,  dot.scale = 5, scale.by = "size") + RotatedAxis()

pdf(paste0("Macrophages-CCKs_DotPlots_byTypeDay.pdf"), height = 5, width = 10)
p
dev.off()

#Feature Plots for ISGs

Idents(BAL.integrated) = "Celltype"

p1 = lapply( ISGs , function(x) { FeaturePlot( BAL.integrated, 
                                               features= x ,
                                               pt.size = 1,
                                               split.by = "Day_Type",
                                               cols = c("grey90","darkgreen") , combine = T,
                                               max.cutoff = "q60") &
    ggplot2::theme(axis.title = element_blank() , 
                   axis.text = element_blank(),
                   axis.ticks = element_blank(),
                   strip.text = element_blank() , 
                   legend.position = "left" , legend.direction="vertical",
                   legend.key.size = unit(1, "in"),legend.text=element_text(size=50, face = "bold"),
                   axis.line   = element_blank() ) & labs(title = "") })

png(paste0("ALLCells-ISGs_FeaturePlots_byTypeDay.png"),
    height = 7000, width = 4600)
patchwork::wrap_plots(p1, nrow = 10, ncol = 1) 
dev.off()

#Feature Plots for IFGs

Idents(BAL.integrated) = "Celltype"

p2 = lapply( IFGs , function(x) { FeaturePlot( BAL.integrated, 
                                              features= x ,
                                              pt.size = 1,
                                              split.by = "Day_Type",
                                              cols = c("grey90","brown1") , combine = T,
                                              max.cutoff = "q60") &
    ggplot2::theme(axis.title = element_blank() , 
                   axis.text = element_blank(),
                   axis.ticks = element_blank(),
                   strip.text = element_blank() , 
                   legend.position = "left" , legend.direction="vertical",legend.key.size = unit(0.35, "in"),
                   axis.line   = element_blank() ) & labs(title = "") })
png(paste0("ALLCells-IFGs_FeaturePlots_byTypeDay.png"),
    height = 2200, width = 1900)
patchwork::wrap_plots(p2, nrow = 5, ncol = 1, widths = 30, heights = 10) 
dev.off()

#Feature Plots for CCKs

Idents(BAL.integrated) = "Celltype"

p3 = lapply( CCK , function(x) { FeaturePlot( BAL.integrated, 
                                               features= x ,
                                               pt.size = 1,
                                               split.by = "Day_Type",
                                               cols = c("grey90","darkorchid1") , combine = T,max.cutoff = "q60") &
    ggplot2::theme(axis.title = element_blank() , 
                   axis.text = element_blank(),
                   axis.ticks = element_blank(),
                   strip.text = element_blank() , 
                   legend.position = "left" , legend.direction="vertical",legend.key.size = unit(0.35, "in"),
                   axis.line   = element_blank() ) & labs(title = "")})

png(paste0("ALLCells-CCKs_FeaturePlots_byTypeDay.png"),
    height = 1500, width = 1400)
patchwork::wrap_plots(p3, nrow = 4, ncol = 1, widths = 30, heights = 10) 
dev.off()

#Violin plots for IFGs

p = VlnPlot(Subsetted_Macrophages, features = IFGs, pt.size = 0, split.by = "Day_Type", ncol =4,
            group.by = "Day_Type", cols = c("#85D1FF","#FFC1AD","#00578B","#FF8C69")) &
  ggplot2::theme( text = element_text(size=20))

pdf(paste0("Macrophages-IFGs_VLnplots_byTypeDay.pdf"), 
    height = 10, width = 20)
p
dev.off()

#Violin plots for ISGs

p = VlnPlot(Subsetted_Macrophages, features = ISGs, pt.size = 0, split.by = "Day_Type", ncol =3,
            group.by = "Day_Type", cols = c("#85D1FF","#FFC1AD","#00578B","#FF8C69"))

pdf(paste0("Macrophages-ISGs_VLnplots_byTypeDay.pdf"), 
    height = 15, width = 20)
p
dev.off()

#Violin plots for CCKs

p = VlnPlot(Subsetted_Macrophages, features = CCK, pt.size = 0, split.by = "Day_Type", ncol =4,
            group.by = "Day_Type", cols = c("#85D1FF","#FFC1AD","#00578B","#FF8C69"))

pdf(paste0("Macrophages-CCKs_VLnplots_byTypeDay.pdf"), 
    height = 3.5, width = 15)
p
dev.off()

##----Get Session info----##
sessionInfo()

