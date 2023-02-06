
library(dplyr)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(Seurat)
library(data.table)
library(Biobase)
library(BisqueRNA)

# Human data deconvolution
# Load human & macaque bulk RNA-seq data
load("apd_macaque_wgcna.rdata")

# Get human DLPFC and sgACC single-nucleus RNA-seq data
hbcc_sn2018 = readRDS("hbcc_sn2018.rds")

# Reformat data, so that sample name can be parsed from cell name
Idents(hbcc_sn2018) = "broad.class" #default is narrow class "ID"
hbcc_newnames = names(hbcc_sn2018$orig.ident)
hbcc_newnames = sub("DLPFC_","DLPFC-",hbcc_newnames)
hbcc_newnames = sub("sgACC_","sgACC-",hbcc_newnames)
hbcc_sn2018 = RenameCells(hbcc_sn2018, new.names = hbcc_newnames)

# Define single-cell expression set
hbcc_sc.eset = SeuratToExpressionSet(hbcc_sn2018, delimiter="-", position=1, version="v3")

# Define bulk expression set
pfc_bulk.eset = ExpressionSet(assayData = as.matrix(cmc_hbcc_dlpfc_f))

# Run Bisque deconvolution for human data (needs more RAM)
# 5 samples overlap and will be used to improve deconvolution
res_decon_pfc_hbcc = ReferenceBasedDecomposition(pfc_bulk.eset, hbcc_sc.eset, markers=NULL, use.overlap=TRUE)

# Filter 5 overlapping samples from metadata
cmc_metadata_dlpfc_ff = cmc_metadata_dlpfc_f[match(colnames(res_decon_hbcc$bulk.props), cmc_metadata_dlpfc_f$BrNum_Region),]

# Summarize deconvolution result
celltype_decon = t(res_decon_hbcc$bulk.props)
colnames(celltype_decon) = c("ExN", "InN", "Astro", "Oligo", "OPC", "Micro", "Endo")
cmc_metadata_dlpfc_ff = cbind(cmc_metadata_dlpfc_ff, celltype_decon)

# Add diagnosis info
df_dlpfc_celltypeDx = cmc_metadata_dlpfc_ff %>% filter(rnaSeq_dissection.Brain_Region=="DLPFC") %>% select(BrNum_Region, Dx, colnames(celltype_decon))
df_dlpfc_celltypeDx = reshape2::melt(df_dlpfc_celltypeDx)
colnames(df_dlpfc_celltypeDx)[3:4] = c("CellType", "Coef")

# Plot SCZ vs. Controls cell type estimates
ggplot(data = df_dlpfc_celltypeDx,  aes(x=Dx, y=Coef, fill=Dx))  + geom_boxplot(outlier.shape = NA) + scale_fill_brewer(palette = "Set2") + geom_jitter(position = position_jitter(.1)) + facet_wrap(~CellType, ncol = 7) +
  theme_bw() + theme(axis.text.x = element_text(size=12, angle = 90, hjust = 1, vjust = .5)) +
  stat_compare_means(ref.group = "Control", method = "wilcox.test", aes(label=..p.signif..))

df_dlpfc_celltypeDxAPD_class = cmc_metadata_dlpfc_ff  %>% select(BrNum_Region, Dx, DxAPD_class, colnames(celltype_decon))
df_dlpfc_celltypeDxAPD_class = reshape2::melt(df_dlpfc_celltypeDxAPD_class)
colnames(df_dlpfc_celltypeDxAPD_class)[4:5] = c("CellType", "Coef")

# Plot SCZ toxicology subgroups cell type estimates
ggplot(data = df_dlpfc_celltypeDxAPD_class,  aes(x=DxAPD_class, y=Coef, fill=Dx))  + geom_boxplot(outlier.shape = NA) + scale_fill_brewer(palette = "Set2") + geom_jitter(position = position_jitter(.1)) + facet_wrap(~CellType, ncol = 7, scales = "free") +
  theme_bw() + theme(axis.text.x = element_text(size=12, angle = 90, hjust = 1, vjust = .5))

# Restrict to 3 cell types for main Fig. 3a
idx_dlpfc_3types = df_dlpfc_celltypeDxAPD_class$CellType %in% c("ExN","InN","Astro")

p_fig3a = ggplot(data = droplevels(df_dlpfc_celltypeDxAPD_class[idx_dlpfc_3types,]),  aes(x=DxAPD_class, y=Coef, fill=Dx))  + geom_boxplot(outlier.shape = NA) + scale_fill_brewer(palette = "Set2") + geom_jitter(position = position_jitter(.1)) + facet_wrap(~CellType, ncol = 7, scales = "free") +
  theme_bw() + theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1, vjust = 1)) + xlab("") +ylab("Estimated cell type proportion")


# Test significance (non-parametric)
a=data.frame(celltype_decon,group=cmc_metadata_dlpfc_ff$DxAPD_class, Dx=cmc_metadata_dlpfc_ff$Dx)
kruskal.test(data=a, ExN~group) # KW p=0.9082
kruskal.test(data=droplevels(a[a$Dx!="Control",]), ExN~group)  # 0.7976

kruskal.test(data=a, InN~group) # KW p=0.03183
kruskal.test(data=droplevels(a[a$Dx!="Control",]), InN~group)  # 0.8753

kruskal.test(data=a, Astro~group) # KW p=0.06894
kruskal.test(data=droplevels(a[a$Dx!="Control",]), Astro~group)  # 0.5648

wilcox.test(data=a, ExN~Dx) # MWU p=0.6957
wilcox.test(data=a, InN~Dx) # MWU p=0.001587
wilcox.test(data=a, Astro~Dx) # MWU p=0.01762

# Macaque data deconvolution
# http://www.evolution.psychencode.org/
# Data from Zhu et al. 2018

load("Sestan.adultMonkeyNuclei.Psychencode.Rdata")
sc_pfc_mmul = CreateSeuratObject(counts = umi2, meta.data = meta2, project = "Monkey_DLPFC")

# Seurat pipeline for macaque data
sc_pfc_mmul = NormalizeData(sc_pfc_mmul, normalization.method = "LogNormalize", scale.factor = 10000)  %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000)  %>% ScaleData() %>%  RunPCA(npcs = 50, ndims.print = 1:5, nfeatures.print = 5)
ElbowPlot(sc_pfc_mmul, ndims = 50)

# UMAP plot of macaque data
set.seed(123)
sc_pfc_mmul = RunUMAP(sc_pfc_mmul, dims = 1:20)
DimPlot(sc_pfc_mmul, reduction = "umap", pt.size = 0.1, group.by = "subtype", label = T)

# Reformat data
# Merge vacular cells (Endo+Peri)
# Filter two clusters that did not map with the rest (Astro2, ExN9)
sc_pfc_mmul$broadclass = sub("[0-9]*$","",sc_pfc_mmul$subtype)
sc_pfc_mmul$broadclass = sub("Endo|Peri","EndoPeri",sc_pfc_mmul$broadclass)
sc_pfc_mmul$broadclass = factor(sc_pfc_mmul$broadclass, levels = c("ExN","InN","Astro","Oligo","OPC","EndoPeri"))
sc_pfc_mmul_f = subset(sc_pfc_mmul, subset = subtype %in% c("Astro2","ExN9"), invert = TRUE)

# Define bulk expression set
mmul_bulk.eset = Biobase::ExpressionSet(assayData = as.matrix(cmc_macaque_countdata))

# Define single-cell expression set
Idents(sc_pfc_mmul_f) = "broadclass" 
mmul_sc.eset = SeuratToExpressionSet(sc_pfc_mmul_f, delimiter="_", position=1, version="v3")
rownames(mmul_sc.eset) = sub(".*-ENSMMUG","ENSMMUG",rownames(mmul_sc.eset))

## Run Bisque deconvolution for Monkey (needs more RAM)
res_decon_mmul = ReferenceBasedDecomposition(bulk.eset = mmul_bulk.eset, sc.eset = mmul_sc.eset, markers=NULL, use.overlap=F)

# Summarize deconvolution result
mmul_celltype_decon = t(res_decon_mmul$bulk.props)
mmul_apd_celltype = data.frame(APD = mmul_metadata_simple$group, mmul_celltype_decon)

# Add APD group info
mmul_apd_celltype$APD = plyr::mapvalues(mmul_metadata_simple$group, from = c("placebo","clozapine","haloperidol_low","haloperidol_high"), to = c("placebo","CLZ","HAL.lo","HAL.hi"))

df_mmul_celltypeAPD = reshape2::melt(mmul_apd_celltype)
colnames(df_mmul_celltypeAPD)[2:3] = c("CellType", "Coef")

# Plot APD groups cell type estimates
ggplot(data = df_mmul_celltypeAPD,  aes(x=APD, y=Coef, fill=APD))  + geom_boxplot(outlier.shape = NA) + geom_jitter(position = position_jitter(.1)) + facet_wrap(~CellType, ncol = 7, scales = "free") +
  theme_bw() + theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1, vjust = 1)) + scale_fill_manual(values = RColorBrewer::brewer.pal(8, "Set2")[c(1,4,3,3)])

# Restrict to 3 cell types for main Fig. 3b
idx_mmul_3types = df_mmul_celltypeAPD$CellType %in% c("ExN","InN","Astro")

p_fig3b = ggplot(data = droplevels(df_mmul_celltypeAPD[idx_mmul_3types,]),  aes(x=APD, y=Coef, fill=APD))  + geom_boxplot(outlier.shape = NA) + geom_jitter(position = position_jitter(.1)) + facet_wrap(~CellType, ncol = 3, scales = "free") +
  theme_bw() + theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1, vjust = 1)) + scale_fill_manual(values = RColorBrewer::brewer.pal(8, "Set2")[c(1,4,3,3)]) + xlab("") +ylab("Estimated cell type proportion")

pdf("fig3_apd_decon.pdf", w=7,h=10, useDingbats = F)
plot_grid(p_fig3a, p_fig3b, labels = c("a","b"), label_size = 18, ncol=1, rel_heights = c(1.2,1))
dev.off()

# Test significance (non-parametric)
a = data.frame(APD = mmul_metadata_simple$group, mmul_celltype_decon)
a$APD = plyr::mapvalues(a$APD, from = c("placebo","clozapine","haloperidol_low","haloperidol_high"), to = c("placebo","CLZ","HAL.lo","HAL.hi"))

kruskal.test(data=a, ExN~APD) # KW p=0.0433
# Fewer excitatory neurons in high-dose haloperidol 
kruskal.test(data=a, InN~APD) # KW p=0.5471
kruskal.test(data=a, Astro~APD) # KW p=0.3355

# Supplemnts for cell type deconvolution

# Seurat pipeline for HBCC data
hbcc_sn2018 = NormalizeData(hbcc_sn2018, normalization.method = "LogNormalize", scale.factor = 10000)  %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000)  %>% ScaleData() %>%  RunPCA(npcs = 50, ndims.print = 1:5, nfeatures.print = 5)
ElbowPlot(sc_pfc_mmul, ndims = 50)

# UMAP for HBCC
set.seed(123)
hbcc_sn2018 = RunUMAP(hbcc_sn2018, dims = 1:20)
hbcc_sn2018$cluster_id = sub("Oligo Pre","OPC",hbcc_sn2018$ID)

p_umap_hbcc = DimPlot(hbcc_sn2018, reduction = "umap", pt.size = 0.1, group.by = "cluster_id", label = T, raster = 10^4, label.size = 3) + NoLegend() + ggtitle("Human single-nucleus RNA-seq") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 15))

p_umap_mmul = DimPlot(sc_pfc_mmul, reduction = "umap", pt.size = 0.1, group.by = "subtype", label = T, raster = 10^4, label.size = 3) + NoLegend() + ggtitle("Macaque single-nucleus RNA-seq") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 15))

# Plot human cell type estimates (SCZ vs. controls for all) 
p_sup_ctype_hbcc = ggplot(data = df_dlpfc_celltypeDx,  aes(x=Dx, y=Coef, fill=Dx))  + geom_boxplot(outlier.shape = NA) + scale_fill_brewer(palette = "Set2") + geom_jitter(position = position_jitter(.1)) + facet_wrap(~CellType, ncol = 7, scales = "free") +
  xlab("") +ylab("Estimated cell type proportion") + ggtitle("Human bulk RNA-seq deconvolution") + 
  theme_bw() + theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1, vjust = 1), plot.title = element_text(hjust = 0.5, face = "bold", size = 15))

# Plot macaque cell type estimates for APD groups (4 remaining cell types not shown in main Figure 3b) 
idx_mmul_4types = !df_mmul_celltypeAPD$CellType %in% c("ExN","InN","Astro")

p_sup_ctype_mmul = ggplot(data = droplevels(df_mmul_celltypeAPD[idx_mmul_4types,]),  aes(x=APD, y=Coef, fill=APD))  + geom_boxplot(outlier.shape = NA) + geom_jitter(position = position_jitter(.1)) + facet_wrap(~CellType, ncol = 3, scales = "free") +
  xlab("") +ylab("Estimated cell type proportion") + ggtitle("Macaque bulk RNA-seq deconvolution") + 
  theme_bw() + theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1, vjust = 1), plot.title = element_text(hjust = 0.5, face = "bold", size = 15)) + scale_fill_manual(values = RColorBrewer::brewer.pal(8, "Set2")[c(1,4,3,3)])

p_umaps = plot_grid(p_umap_hbcc, p_umap_mmul,align = "h", rel_widths = c(1.2,1), labels = c("a","b"), label_size = 18)
p_sup_ctype = plot_grid(p_sup_ctype_hbcc, p_sup_ctype_mmul, align = "h", rel_widths = c(1.4,1), labels = c("c","d"), label_size = 18)

pdf("sup_fig_umaps_ctypes.pdf", w=15,h=12, useDingbats = F)
plot_grid(p_umaps, p_sup_ctype, align = "v", label_size = 18, ncol = 1)
dev.off()


