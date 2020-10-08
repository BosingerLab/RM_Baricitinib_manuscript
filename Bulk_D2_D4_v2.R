library(DESeq2)
library(ggplot2)
library(patchwork)
library(ggrepel)
library(ComplexHeatmap)
library(circlize)
library("tidyverse")
library("reshape2")


setwd("~/Analysis")
path <- "~/Data"

sampleTable <- read.table("samplesheet_Bulk_baseline_D2_D4.txt", header = TRUE, row.names = 1, check.names = FALSE)

dds <- DESeqDataSetFromHTSeqCount(sampleTable, path, design = ~ Subject + Group)

# Keep protein coding genes only
prot_coding <- as.character(unlist(read.table("Data/protein_coding.txt")))
keep <- row.names(dds) %in% prot_coding
dds <- dds[keep,]

dds <- DESeq(dds, parallel = TRUE)
rnames <- resultsNames(dds)

rld <- rlog(dds, blind = FALSE, fitType = "parametric")
write.table(assay(rld), file = "BAL_rlog.txt", sep = "\t", quote = FALSE)

norm_counts <- counts(dds,normalized=TRUE)
write.table(norm_counts,"BAL_norm_counts.txt",quote=FALSE,sep="\t")

counts <- counts(dds,normalized=FALSE)
write.table(counts,"BAL_counts.txt",quote=FALSE,sep="\t")

############################################################################################################
# Get rlog files
############################################################################################################

st <- sampleTable
st$ID <- paste(st$Condition,st$AnimalID,st$Timepoint,sep=".")
head(st)

identical(colnames(rld),row.names(st))
colnames(rld)<- as.character(st$ID)

rld_filt <- rld[rowMeans(assay(rld)) > 0,]

write.table(assay(rld_filt)[,grepl("Treated.*baseline|Treated.*2dpi",colnames(rld_filt))], "rlog_Treated_2dpi_baseline_filt.txt",quote=FALSE,sep="\t")
write.table(assay(rld_filt)[,grepl("Treated.*baseline|Treated.*4dpi",colnames(rld_filt))], "rlog_Treated_4dpi_baseline_filt.txt",quote=FALSE,sep="\t")
write.table(assay(rld_filt)[,grepl("Treated.*2dpi|Treated.*4dpi",colnames(rld_filt))], "rlog_Treated_4dpi_2dpi_filt.txt",quote=FALSE,sep="\t")

write.table(assay(rld_filt)[,grepl("Untreated.*baseline|Untreated.*2dpi",colnames(rld_filt))], "rlog_Untreated_2dpi_baseline_filt.txt",quote=FALSE,sep="\t")
write.table(assay(rld_filt)[,grepl("Untreated.*baseline|Untreated.*4dpi",colnames(rld_filt))], "rlog_Untreated_4dpi_baseline_filt.txt",quote=FALSE,sep="\t")
write.table(assay(rld_filt)[,grepl("Untreated.*2dpi|Untreated.*4dpi",colnames(rld_filt))], "rlog_Untreated_4dpi_2dpi_filt.txt",quote=FALSE,sep="\t")

write.table(assay(rld_filt)[,grepl("4dpi",colnames(rld_filt))], "rlog_4dpi_Treated_Untreated_filt.txt",quote=FALSE,sep="\t")
write.table(assay(rld_filt)[,grepl("2dpi",colnames(rld_filt))], "rlog_2dpi_Treated_Untreated_filt.txt",quote=FALSE,sep="\t")

############################################################################################################
# DE analysis
############################################################################################################

#-----------------------------------------------------------------------------------------------------------
# Treated
#-----------------------------------------------------------------------------------------------------------
# 2dpi
res_trt_2dpi <- results(dds, contrast = c("Group","Treated_2dpi","Treated_baseline"), alpha = 0.05, lfcThreshold = log2(1.5))
summary(res_trt_2dpi)
res_trt_2dpi <- res_trt_2dpi[order(res_trt_2dpi$pvalue),]
write.table(res_trt_2dpi, file = "BAL_treated_2dpi_baseline_DEG.txt", sep = "\t", quote = FALSE)
res_trt_2dpi_sig <- na.omit(res_trt_2dpi)
res_trt_2dpi_sig <- res_trt_2dpi_sig[res_trt_2dpi_sig$padj < 0.05 & abs(res_trt_2dpi_sig$log2FoldChange) > log2(1.5) & res_trt_2dpi_sig$lfcSE < 1,]
write.table(res_trt_2dpi_sig, file = "BAL_treated_2dpi_baseline_DEG_sig.txt", sep = "\t", quote = FALSE)


# 4dpi
res_trt_4dpi <- results(dds, contrast = c("Group","Treated_4dpi","Treated_baseline"), alpha = 0.05, lfcThreshold = log2(1.5))
summary(res_trt_4dpi)
res_trt_4dpi <- res_trt_4dpi[order(res_trt_4dpi$pvalue),]
write.table(res_trt_4dpi, file = "BAL_treated_4dpi_baseline_DEG.txt", sep = "\t", quote = FALSE)
res_trt_4dpi_sig <- na.omit(res_trt_4dpi)
res_trt_4dpi_sig <- res_trt_4dpi_sig[res_trt_4dpi_sig$padj < 0.05 & abs(res_trt_4dpi_sig$log2FoldChange) > log2(1.5) & res_trt_4dpi_sig$lfcSE < 1,]
write.table(res_trt_4dpi_sig, file = "BAL_treated_4dpi_baseline_DEG_sig.txt", sep = "\t", quote = FALSE)


# 4dpi_2dpi
res_trt_4dpi_2dpi <- results(dds, contrast = c("Group","Treated_4dpi","Treated_2dpi"), alpha = 0.05, lfcThreshold = log2(1.5))
summary(res_trt_4dpi_2dpi)
res_trt_4dpi_2dpi <- res_trt_4dpi_2dpi[order(res_trt_4dpi_2dpi$pvalue),]
write.table(res_trt_4dpi_2dpi, file = "BAL_treated_4dpi_2dpi_baseline_DEG.txt", sep = "\t", quote = FALSE)
res_trt_4dpi_2dpi_sig <- na.omit(res_trt_4dpi_2dpi)
res_trt_4dpi_2dpi_sig <- res_trt_4dpi_2dpi_sig[res_trt_4dpi_2dpi_sig$padj < 0.05 & abs(res_trt_4dpi_2dpi_sig$log2FoldChange) > log2(1.5) & res_trt_4dpi_2dpi_sig$lfcSE < 1,]
write.table(res_trt_4dpi_2dpi_sig, file = "BAL_treated_4dpi_2dpi_DEG_sig.txt", sep = "\t", quote = FALSE)

#-----------------------------------------------------------------------------------------------------------
# Untreated
#-----------------------------------------------------------------------------------------------------------

# 2dpi
res_untrt_2dpi <- results(dds, contrast = c("Group","Untreated_2dpi","Untreated_baseline"), alpha = 0.05, lfcThreshold = log2(1.5))
summary(res_untrt_2dpi)
res_untrt_2dpi <- res_untrt_2dpi[order(res_untrt_2dpi$pvalue),]
write.table(res_untrt_2dpi, file = "BAL_untreated_2dpi_baseline_DEG.txt", sep = "\t", quote = FALSE)
res_untrt_2dpi_sig <- na.omit(res_untrt_2dpi)
res_untrt_2dpi_sig <- res_untrt_2dpi_sig[res_untrt_2dpi_sig$padj < 0.05 & abs(res_untrt_2dpi_sig$log2FoldChange) > log2(1.5) & res_untrt_2dpi_sig$lfcSE < 1,]
write.table(res_untrt_2dpi_sig, file = "BAL_untreated_2dpi_baseline_DEG_sig.txt", sep = "\t", quote = FALSE)

# 4dpi
res_untrt_4dpi <- results(dds, contrast = c("Group","Untreated_4dpi","Untreated_baseline"), alpha = 0.05, lfcThreshold = log2(1.5))
summary(res_untrt_4dpi)
res_untrt_4dpi <- res_untrt_4dpi[order(res_untrt_4dpi$pvalue),]
write.table(res_untrt_4dpi, file = "BAL_untreated_4dpi_baseline_DEG.txt", sep = "\t", quote = FALSE)
res_untrt_4dpi_sig <- na.omit(res_untrt_4dpi)
res_untrt_4dpi_sig <- res_untrt_4dpi_sig[res_untrt_4dpi_sig$padj < 0.05 & abs(res_untrt_4dpi_sig$log2FoldChange) > log2(1.5) & res_untrt_4dpi_sig$lfcSE < 1,]
write.table(res_untrt_4dpi_sig, file = "BAL_untreated_4dpi_baseline_DEG_sig.txt", sep = "\t", quote = FALSE)

# 4dpi vs 2dpi
res_untrt_4dpi_2dpi <- results(dds, contrast = c("Group","Untreated_4dpi","Untreated_2dpi"), alpha = 0.05, lfcThreshold = log2(1.5))
summary(res_untrt_4dpi_2dpi)
res_untrt_4dpi_2dpi <- res_untrt_4dpi_2dpi[order(res_untrt_4dpi_2dpi$pvalue),]
write.table(res_untrt_4dpi_2dpi, file = "BAL_untreated_4dpi_2dpi_baseline_DEG.txt", sep = "\t", quote = FALSE)
res_untrt_4dpi_2dpi_sig <- na.omit(res_untrt_4dpi_2dpi)
res_untrt_4dpi_2dpi_sig <- res_untrt_4dpi_2dpi_sig[res_untrt_4dpi_2dpi_sig$padj < 0.05 & abs(res_untrt_4dpi_2dpi_sig$log2FoldChange) > log2(1.5) & res_untrt_4dpi_2dpi_sig$lfcSE < 1,]
write.table(res_untrt_4dpi_2dpi_sig, file = "BAL_untreated_4dpi_2dpi_DEG_sig.txt", sep = "\t", quote = FALSE)

############################################################################################################
# Heatmap - normalize by median of all baseline samples
############################################################################################################

baseline_medians <- rowMedians(assay(rld[,grepl("baseline",colnames(rld))]))
hdata <- assay(rld) - baseline_medians

col.order.all <- c("Untreated.6_112.baseline","Untreated.5_215.baseline","Untreated.RQv9.baseline","Untreated.RHz12.baseline",
                   "Untreated.6_112.2dpi","Untreated.5_215.2dpi","Untreated.RQv9.2dpi","Untreated.RHz12.2dpi",
                   "Untreated.6_112.4dpi","Untreated.5_215.4dpi","Untreated.RQv9.4dpi","Untreated.RHz12.4dpi",
                   "Treated.RAt11.baseline","Treated.RLf10.baseline","Treated.7_141.baseline","Treated.RVf12.baseline",
                   "Treated.RAt11.2dpi","Treated.RLf10.2dpi","Treated.7_141.2dpi","Treated.RVf12.2dpi",
                   "Treated.RAt11.4dpi","Treated.RLf10.4dpi","Treated.7_141.4dpi","Treated.RVf12.4dpi")

path_leading <- "~/Analysis/GSEA/shortlist/Untreated_4dpi_2dpi.Gsea.1597974173866/"  # Fig 3
path_leading <- "~/Analysis/GSEA/shortlist/Untreated_4dpi_baseline.Gsea.1597974009109/"  # Fig S5d

# Neut degranulation
genes <- read.table(paste0(path_leading,"REACTOME_NEUTROPHIL_DEGRANULATION.tsv"), header=TRUE,sep="\t")
genes <- as.character(genes[genes$CORE.ENRICHMENT =="Yes",]$SYMBOL[1:35])
hd <- hdata[genes,col.order.all]

col_fun = colorRamp2(c(-3, 0, 3), c("blue", "white", "red"))

Heatmap(hd, 
        column_order = col.order.all, 
        show_column_dend = FALSE,
        cluster_columns = FALSE,
        show_row_dend = FALSE,
        cluster_rows = FALSE,
        col = col_fun,
        column_names_gp = gpar(fontsize = 10), 
        row_names_gp = gpar(fontsize = 10), 
        show_column_names = FALSE,
        column_names_rot = 90,
        rect_gp = gpar(col = "black", lwd = 1),
        column_split= factor(c(rep("Untreated\nbaseline",4),rep("Untreated\n2dpi",4),rep("Untreated\n4dpi",4),
                               rep("Treated\nbaseline",4),rep("Treated\n2dpi",4),rep("Treated\n4dpi",4)), 
                             levels=c("Untreated\nbaseline","Untreated\n2dpi","Untreated\n4dpi",
                                      "Treated\nbaseline","Treated\n2dpi","Treated\n4dpi")),
        column_gap = unit(2,"mm"),
        border = TRUE,
        heatmap_width = unit(5,"in"),
        heatmap_height = unit(6,"in")
) 

# TNFA
genes <- read.table(paste0(path_leading,"HALLMARK_TNFA_SIGNALING_VIA_NFKB.tsv"), header=TRUE,sep="\t")
genes <- as.character(genes[genes$CORE.ENRICHMENT =="Yes",]$SYMBOL[1:35])
hd <- hdata[genes,col.order.all]

Heatmap(hd, 
        column_order = col.order.all, 
        show_column_dend = FALSE,
        cluster_columns = FALSE,
        show_row_dend = FALSE,
        cluster_rows = FALSE,
        col = col_fun,
        column_names_gp = gpar(fontsize = 10), 
        row_names_gp = gpar(fontsize = 10), 
        show_column_names = FALSE,
        column_names_rot = 90,
        rect_gp = gpar(col = "black", lwd = 1),
        column_split= factor(c(rep("Untreated\nbaseline",4),rep("Untreated\n2dpi",4),rep("Untreated\n4dpi",4),
                               rep("Treated\nbaseline",4),rep("Treated\n2dpi",4),rep("Treated\n4dpi",4)), 
                             levels=c("Untreated\nbaseline","Untreated\n2dpi","Untreated\n4dpi",
                                      "Treated\nbaseline","Treated\n2dpi","Treated\n4dpi")),
        column_gap = unit(2,"mm"),
        border = TRUE,
        heatmap_width = unit(5,"in"),
        heatmap_height = unit(6,"in")
) 


# ISG
genes <- read.table(paste0(path_leading,"ISG.tsv"), header=TRUE,sep="\t")
genes <- as.character(genes[genes$CORE.ENRICHMENT =="Yes",]$SYMBOL)
genes
hd <- hdata[genes,col.order.all]

Heatmap(hd, 
        column_order = col.order.all, 
        show_column_dend = FALSE,
        cluster_columns = FALSE,
        show_row_dend = FALSE,
        cluster_rows = FALSE,
        col = col_fun,
        column_names_gp = gpar(fontsize = 10), 
        row_names_gp = gpar(fontsize = 10), 
        show_column_names = FALSE,
        column_names_rot = 90,
        rect_gp = gpar(col = "black", lwd = 1),
        column_split= factor(c(rep("Untreated\nbaseline",4),rep("Untreated\n2dpi",4),rep("Untreated\n4dpi",4),
                               rep("Treated\nbaseline",4),rep("Treated\n2dpi",4),rep("Treated\n4dpi",4)), 
                             levels=c("Untreated\nbaseline","Untreated\n2dpi","Untreated\n4dpi",
                                      "Treated\nbaseline","Treated\n2dpi","Treated\n4dpi")),
        column_gap = unit(2,"mm"),
        border = TRUE,
        heatmap_width = unit(5,"in"),
        heatmap_height = unit(3,"in")
) 


#IL6 JAK STAT

genes <- read.table(paste0(path_leading,"HALLMARK_IL6_JAK_STAT3_SIGNALING.tsv"), header=TRUE,sep="\t")
genes <- as.character(genes[genes$CORE.ENRICHMENT =="Yes",]$SYMBOL)

hd <- hdata[genes,col.order.all]

Heatmap(hd, 
        column_order = col.order.all, 
        show_column_dend = FALSE,
        cluster_columns = FALSE,
        show_row_dend = FALSE,
        cluster_rows = FALSE,
        col = col_fun,
        column_names_gp = gpar(fontsize = 10), 
        row_names_gp = gpar(fontsize = 10), 
        show_column_names = FALSE,
        column_names_rot = 90,
        rect_gp = gpar(col = "black", lwd = 1),
        column_split= factor(c(rep("Untreated\nbaseline",4),rep("Untreated\n2dpi",4),rep("Untreated\n4dpi",4),
                               rep("Treated\nbaseline",4),rep("Treated\n2dpi",4),rep("Treated\n4dpi",4)), 
                             levels=c("Untreated\nbaseline","Untreated\n2dpi","Untreated\n4dpi",
                                      "Treated\nbaseline","Treated\n2dpi","Treated\n4dpi")),
        column_gap = unit(2,"mm"),
        border = TRUE,
        heatmap_width = unit(5,"in"),
        heatmap_height = unit(4.11,"in")
) 


############################################################################################################
# NES barplot
############################################################################################################

nes_path <- "~/Analysis/GSEA/NES.txt";

d<-read.table(nes_path, header = FALSE)
head(d)
names(d) <- c("Geneset","Condition","NES","NOM p-val","FDR q-val")
d$Condition <- factor(d$Condition,levels=c("Treated","Untreated"))
d$`NOM p-val`<- round(d$`NOM p-val`, digits = 3)

x <- d[d$Condition=="Untreated" & d$Condition=="Untreated",]
geneset_order <- as.character(x[order(x$NES),]$Geneset)  
d$Geneset <- factor(d$Geneset, levels=geneset_order)
head(d)

d$`NOM p-val`[d$`NOM p-val`==0] = "<0.001"

ggplot(d,aes(Geneset,NES, fill=Condition)) + 
  geom_bar(stat = "identity", width = 0.8, position = "dodge") + coord_flip() +
  scale_fill_manual(values=c("#00578b","#ff8c69")) + 
  geom_text_repel(aes(label=`NOM p-val`),position=position_dodge(width=0.8), 
            direction = "x",
            hjust = -0.5,
            segment.size = 0) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 9, angle = 0, hjust = 1, color="black"),
        axis.text.y = element_text(size = 9, color="black"),
        axis.title = element_text(size = 14, color="black"),
        strip.text = element_text(size = 12, color="black"),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        line = element_line(colour = "black"),
        strip.background =element_rect(fill="white")) +
  ylim(c(-2.5,2.5))


############################################################################################################
# Line plots
############################################################################################################

st <- sampleTable
st$ID <- paste(st$Condition,st$AnimalID,st$Timepoint,sep=".")

identical(colnames(norm_counts),row.names(st))
colnames(norm_counts) <- st$ID

ra <- read.table("~/baseline_D2_D4_v2/GSEA/shortlist/D4_Untreated_Treated.Gsea.1597974108091/RHEUMATOID_ARTHRITIS.tsv",sep="\t", header=TRUE)

genes <- as.character(ra[ra$CORE.ENRICHMENT == "Yes",]$SYMBOL)
length(genes)

d <- t(norm_counts[genes,])
d <- melt(d)
names(d) <- c("Sample","Gene","Normalized expression")
head(d)
d <- d %>% separate("Sample",into =c("Condition","AnimalID","Timepoint"),sep="\\.", remove = FALSE)
d <- d[d$Timepoint == "4dpi" | d$Timepoint == "2dpi",]
d <- droplevels(d)
d$log10Expression <- log10(d$`Normalized expression` + 1)
d$Timepoint <- factor(d$Timepoint, levels=c("2dpi","4dpi"))
d$Condition <- factor(d$Condition, levels = c("Untreated","Treated"))

gd <- d %>% 
  group_by(Condition,Timepoint,Gene) %>% 
  summarise(Mean = mean(`Normalized expression`))
head(gd)

gd <- gd[gd$Timepoint=="4dpi",]

ggplot(gd,aes(Condition,log10(Mean))) + 
  geom_point(aes(color=Condition, fill=Condition),shape=21, size = 4) + 
  geom_line(aes(group=Gene)) + theme_bw() +
  geom_text_repel(data = subset(gd,Condition=="Treated"), aes(label=Gene),
                  direction= "y",
                  hjust        = 0, 
                  segment.size = 0.2,
                  nudge_x=0.15, force = 0.5,
                  segment.colour = "gray50"
                  ) +
  scale_fill_manual(values=c("#ff8c69","#00578b")) + 
  scale_color_manual(values=c("#ff8c69","#00578b"))  +
  ylab("log10 Average Normalized expression") +
  theme(axis.text.x = element_text(size = 12, angle = 00, color="black"),
        axis.text.y = element_text(size = 12, color="black"),
        axis.title = element_text(size = 14, color="black"),
        strip.text = element_text(size = 14, color="black"),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        line = element_line(colour = "black"),
        strip.background =element_rect(fill="white"),
        legend.position = "None")  


gd <- d %>% 
  group_by(Condition,Timepoint,Gene) %>% 
  summarise(Mean = mean(`Normalized expression`))
head(gd)

gd <- gd[gd$Timepoint=="2dpi",]

ggplot(gd,aes(Condition,log10(Mean))) + 
  geom_point(aes(color=Condition, fill=Condition),shape=21, size = 4) + 
  geom_line(aes(group=Gene)) + theme_bw() +
  geom_text_repel(data = subset(gd,Condition=="Treated"), aes(label=Gene),
                  direction= "y",
                  hjust        = 0, 
                  segment.size = 0.2,
                  nudge_x=0.15, force = 0.5,
                  segment.colour = "gray50"
  ) +
  scale_fill_manual(values=c("#ff8c69","#00578b")) + 
  scale_color_manual(values=c("#ff8c69","#00578b"))  +
  ylab("log10 Average Normalized expression") +
  theme(axis.text.x = element_text(size = 12, angle = 00, color="black"),
        axis.text.y = element_text(size = 12, color="black"),
        axis.title = element_text(size = 14, color="black"),
        strip.text = element_text(size = 14, color="black"),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        line = element_line(colour = "black"),
        strip.background =element_rect(fill="white"),
        legend.position = "None")  
  

