
library(dplyr)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(edgeR)

# Human DLPFC data from CMC_HBCC
cmc_hbcc_dlpfc = read.delim("CMC_HBCC_DLPFC_Counts.tsv", row.names = 1)
cmc_metadata = read.csv("CMC_Human_rnaSeq_mergedMetadata.csv")

cmc_metadata_dlpfc = filter(.data = cmc_metadata, Institution == "NIMH-HBCC" & rnaSeq_dissection.Brain_Region=="DLPFC")
cmc_metadata_dlpfc = cmc_metadata_dlpfc[match(colnames(cmc_hbcc_dlpfc),cmc_metadata_dlpfc$Sample_RNA_ID),]

# Get gene models from BioMart and remove duplicate gene symbols
genes_ens86 = read.csv("ens_genes86.csv")
all_names_dlpfc = genes_ens86$Gene.name[match(sub("\\..*","",rownames(cmc_hbcc_dlpfc)), genes_ens86$Gene.stable.ID)]
dup_names_dlpfc = unique(all_names_dlpfc[duplicated(all_names_dlpfc)])
dup_ensid_dlpfc = genes_ens86$Gene.stable.ID[genes_ens86$Gene.name %in% dup_names_dlpfc]
cmc_hbcc_dlpfc_nodup = cmc_hbcc_dlpfc[!sub("\\..*","",rownames(cmc_hbcc_dlpfc)) %in% dup_ensid_dlpfc,]
rownames(cmc_hbcc_dlpfc_nodup) = genes_ens86$Gene.name[match(sub("\\..*","",rownames(cmc_hbcc_dlpfc_nodup)), genes_ens86$Gene.stable.ID)]

# Add general APD toxicology info
BrNum_hbcc_cmc = read.csv("CMC_HBCC_CapstoneID_with_BrNum.csv")
hbcc_apd_tox = read.csv("toxData.csv") # all major drug classes
hbcc_apd_tox = droplevels(hbcc_apd_tox[hbcc_apd_tox$Antipsychotics..ZANTIPSYCH. != "NULL",])

cmc_metadata_dlpfc$APD_tox = hbcc_apd_tox$Antipsychotics..ZANTIPSYCH.[match(cmc_metadata_dlpfc$Brain_ID, hbcc_apd_tox$BrainNumber)]
cmc_metadata_dlpfc$APD_tox[is.na(cmc_metadata_dlpfc$APD_tox)] = ""

# Filter data & metadata
idx_szctrl = cmc_metadata_dlpfc$Dx %in% c("Control","SCZ")
idx_pfc_over17 = cmc_metadata_dlpfc$ageOfDeath >= 17
idx_pfc_trep = duplicated(cmc_metadata_dlpfc$Brain_ID)
idx_noAPD # remove samples without toxicology and 2 cases negative in blood but positive in brain
idx_outlier # remove 2 samples defined in outlier analysis below
idx_pfc_f = idx_pfc_over17 & idx_szctrl & (!is.na(cmc_metadata_dlpfc$pH)) & (!idx_pfc_trep &!idx_noAPD &!idx_outlier)

cmc_hbcc_dlpfc_f = cmc_hbcc_dlpfc_nodup[,idx_pfc_f]
colnames(cmc_hbcc_dlpfc_f) = paste0(cmc_metadata_dlpfc$Brain_ID[idx_pfc_f], "_DLPFC")

cmc_metadata_dlpfc_f = cmc_metadata_dlpfc[idx_pfc_f,]
cmc_metadata_dlpfc_f$BrNum_Region = paste0(cmc_metadata_dlpfc_f$Brain_ID, "_", cmc_metadata_dlpfc_f$rnaSeq_dissection.Brain_Region)

# Add APD compounds
apd_types = read.csv("apd_types.csv") # typical vs. atypical for each compound
apd_compounds = read.csv("apd_compounds.csv") # APD comound for each BrNum

apd_types = apd_types[apd_types$APD %in% unique(apd_compounds$compound),]
apd_compounds$APD_class = plyr::mapvalues(x = as.character(apd_compounds$compound), from = as.character(apd_types$APD), to = as.character(apd_types$class))
apd_classes = apd_compounds[,c(1,5)]
apd_classes = apd_classes[!duplicated(apd_classes),]
apd_mixed = as.character(apd_classes$brainnumber[duplicated(apd_classes$brainnumber)])
apd_classes$APD_class[apd_classes$brainnumber %in% apd_mixed] = "mixed"
apd_classes = apd_classes[!duplicated(apd_classes),]

cmc_metadata_dlpfc_f$APD_class = apd_classes$APD_class[match(cmc_metadata_dlpfc_f$Brain_ID, apd_classes$brainnumber)]
cmc_metadata_dlpfc_f$APD_class[is.na(cmc_metadata_dlpfc_f$APD_class)] = ""

cmc_metadata_dlpfc_f$DxAPD = paste0(cmc_metadata_dlpfc_f$Dx,".",sub("itive|ative",".",cmc_metadata_dlpfc_f$APD_tox))
cmc_metadata_dlpfc_f$DxAPD_class = paste0(cmc_metadata_dlpfc_f$DxAPD, cmc_metadata_dlpfc_f$APD_class)

# Simplify metadata
metadata_szctrl_pfc = cmc_metadata_dlpfc_f
metadata_szctrl_pfc = data.frame(row.names = paste0("Br",metadata_szctrl_pfc$Brain_ID),
                                 Dx = factor(metadata_szctrl_pfc$Dx, levels = c("Control","SCZ")),
                                 Sex = factor(metadata_szctrl_pfc$Reported_Gender,levels = c("Male","Female")),
                                 Race = factor(metadata_szctrl_pfc$Ethnicity),
                                 batch_rna = factor(metadata_szctrl_pfc$rnaSeq_isolation.RNA_Isolation_Batch),
                                 batch_lib = factor(metadata_szctrl_pfc$rnaSeq_report.Library_Batch),
                                 batch_seq = factor(metadata_szctrl_pfc$rnaSeq_report.Flowcell_Batch),
                                 Age = metadata_szctrl_pfc$ageOfDeath,
                                 PMI = metadata_szctrl_pfc$PMI_.in_hours.,
                                 pH = metadata_szctrl_pfc$pH,
                                 RIN = metadata_szctrl_pfc$rnaSeq_isolation.RIN,
                                 total_tissue = metadata_szctrl_pfc$rnaSeq_dissection.Tissue_Amount_.grams.,
                                 total_RNA = metadata_szctrl_pfc$rnaSeq_isolation.Total_RNA_Yield,
                                 total_reads = metadata_szctrl_pfc$rnaSeq_report.Total_Reads,
                                 total_genes = metadata_szctrl_pfc$rnaSeq_report.Genes_Detected,
                                 rate_mapped = metadata_szctrl_pfc$rnaSeq_report.Percent_Aligned,
                                 rate_intron = metadata_szctrl_pfc$rnaSeq_report.Intronic_Rate,
                                 rate_intergenic = metadata_szctrl_pfc$rnaSeq_report.Intergenic_Rate,
                                 rate_rRNA = metadata_szctrl_pfc$rnaSeq_report.rRNA_Rate,
                                 rate_efficiency = metadata_szctrl_pfc$rnaSeq_report.Expression_Profiling_Efficiency,
                                 DxAPD = paste0(metadata_szctrl_pfc$Dx,".",sub("itive|ative",".",metadata_szctrl_pfc$APD_tox)))

# Define useful functions

mypca = function(object, col = "black", labels = colnames(object), main = NULL, cex=1) 
{
  pca = prcomp(t(object))
  percentVar = pca$sdev^2/sum(pca$sdev^2)
  plot(pca$x[, 1], pca$x[, 2], col = as.character(col), main = main, pch="",
       xlab = paste0("PC1: ", round(percentVar[1] * 100, digits = 1), "% variance"),
       ylab = paste0("PC2: ", round(percentVar[2] * 100, digits = 1), "% variance"))
  text(pca$x[, 1], pca$x[, 2], labels=labels, col = as.character(col), cex = cex)
  pca$percentVar=percentVar
  return(pca)
}


get_voom=function(counts, genes, design=NULL){
  require(edgeR)
  d=DGEList(counts = counts, genes = genes)
  d=calcNormFactors(d)
  voom(d,plot=T, design = design)
}

getFit=function(voom, design){
  require(edgeR)
  lfit=lmFit(voom, design = design)
  eBayes(lfit)
}


get_vresid = function(fit, voom){
  require(edgeR)
  resid = residuals.MArrayLM(object = fit, y = voom)
  vresid = voom
  vresid$E = resid
  return(vresid)
}

# Run PCA on logCPMs
v_pfc_f = get_voom(counts = cmc_hbcc_dlpfc_f, genes = genes_ens86[match(rownames(cmc_hbcc_dlpfc_f), genes_ens86$Gene.name),])

pca_pfc_f = mypca(object = v_pfc_f$E)

# Explore relationship of PCs and covariates
pval_pfc_f=matrix(ncol=20, nrow=20)
rownames(pval_pfc_f)=colnames(metadata_szctrl_pfc)
colnames(pval_pfc_f)=paste0("PC",1:20)

for (PC in 1:20) {for (varname in colnames(metadata_szctrl_pfc)){
  a=summary(lm(pca_pfc_f$x[,PC]~metadata_szctrl_pfc[,varname]))
  pval_pfc_f[varname,PC]=pf(a$fstatistic[1],a$fstatistic[2],a$fstatistic[3],lower.tail=FALSE)[1]
}}
round(-log10(pval_pfc_f))


##outlier test
#
#apply(abs(scale(pca_pfc_f$x[,1:20])), 2, which.max)
#apply(abs(scale(pca_pfc_f$x[,1:20])), 2, max)
#
#pc_pfc_2sd = abs(scale(pca_pfc_f$x[,1:20]))>2
#pc_pfc_3sd = abs(scale(pca_pfc_f$x[,1:20]))>3
#pc_pfc_5sd = abs(scale(pca_pfc_f$x[,1:20]))>5
#
#hist(rowSums(pc_pfc_2sd))
#hist(rowSums(pc_pfc_3sd))
#hist(rowSums(pc_pfc_5sd))
#
#rownames(pc_pfc_2sd)[rowSums(pc_pfc_2sd)>=7]
#rownames(pc_pfc_3sd)[rowSums(pc_pfc_3sd)>=3]
#rownames(pc_pfc_5sd)[rowSums(pc_pfc_5sd)>=2]
#
## Remove 2 outliers with abs(Z)>5 in 2 or more of the top 20 PCs (highest scores; only ones >5)
## Go back to filter and define idx_outlier to filter these 2 samples

# Define residuals after covariate regression and test DGE in SCZ via limma-voom
des_szctrl_pfc_cov1 = model.matrix(data=metadata_szctrl_pfc,~batch_lib+Sex+Race+Age+RIN+rate_efficiency+rate_intergenic+PMI+pH)
fit_szctrl_pfc_cov1 = getFit(voom = v_pfc_f, design = des_szctrl_pfc_cov1)
resid_szctrl_pfc_cov1 = get_vresid(fit = fit_szctrl_pfc_cov1, voom = v_pfc_f)

des_szctrl_pfc_Dx = model.matrix(data = metadata_szctrl_pfc, ~Dx)
fit_szctrl_pfc_resid1_Dx = getFit(voom = resid_szctrl_pfc_cov1, design = des_szctrl_pfc_Dx)

des_szctrl_pfc_DxAPD = model.matrix(data = metadata_szctrl_pfc, ~DxAPD)
fit_szctrl_pfc_resid1_DxAPD = getFit(voom = resid_szctrl_pfc_cov1, design = des_szctrl_pfc_DxAPD)

tt_HBCC.SCZ = topTable(fit_szctrl_pfc_resid1_Dx, coef = 2, sort.by = "none", n=Inf)

# Get all DE genes (nominal P<.05)
genes_sz_pfc_resid1 = rownames(resid_szctrl_pfc_cov1)[fit_szctrl_pfc_resid1_Dx$p.value[,"DxSCZ"]<.05]

# Define SCZ signature (aggregated DE score) based on raw and residualized gene expression values
sig_sz_pfc_resid1 = scale(t(resid_szctrl_pfc_cov1$E[genes_sz_pfc_resid1,]) %*% fit_szctrl_pfc_resid1_Dx$coefficients[genes_sz_pfc_resid1,"DxSCZ"])
sig_sz_pfc_raw = scale(x = t(v_pfc_f$E[genes_sz_pfc_resid1,]) %*% fit_szctrl_pfc_resid1_Dx$coefficients[genes_sz_pfc_resid1,"DxSCZ"])


# Plots for Fig. 1

df_apdtox = data.frame(SZ_signature = sig_sz_pfc_resid1, group = metadata_szctrl_pfc$DxAPD)

df_apdtox_class = data.frame(Dx = cmc_metadata_dlpfc_f$Dx, SZ_signature = sig_sz_pfc_resid1, group = cmc_metadata_dlpfc_f$DxAPD_class)
#my_comparisons = list(c("Control.Neg.", "SCZ.Neg."),c("Control.Neg.", "SCZ.Pos.atypical"),c("Control.Neg.", "SCZ.Pos.mixed"),c("Control.Neg.", "SCZ.Pos.typical"))
df_apdtox_class$group=factor(df_apdtox_class$group)

p_fig1a = ggplot(df_apdtox_class, aes(x=group, y=SZ_signature, fill=Dx)) + geom_boxplot(outlier.shape = NA) + geom_jitter(position = position_jitter(.1)) + scale_y_continuous() + scale_fill_brewer(palette = "Set2") + theme_classic() + theme(axis.text.x = element_text(size=12, angle = 90, hjust = 1, vjust = .5)) + xlab("") + ylab("SCZ expression signature (Z-score)") +
#  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", p.adjust.method = "BH", aes(label=..p.signif..)) +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1))

# test significance (parametric)
summary(aov(data=df_apdtox_class, SZ_signature~group)) #ANOVA p=1.29e-12 ***
TukeyHSD(aov(data=df_apdtox_class, SZ_signature~group)) 
# Shapiro-Wilks test for normality: p=0.0004608 (different from normal -> use non-parametric)
shapiro.test(sig_sz_pfc_resid1)
# test significance (non-parametric)
kruskal.test(data=df_apdtox_class, SZ_signature~group) # KW p=3.829e-10 ****
kruskal.test(data=droplevels(df_apdtox_class[df_apdtox_class$Dx!="Control",]), SZ_signature~group)  # within SCZ KW p=0.01382 *
FSA::dunnTest(data=df_apdtox_class, SZ_signature~group, method="bh")

# Get individual APDs
apd_compounds2 = apd_compounds[,-2]
apd_compounds2$compound[grep("Risperidone|risperidone")]
apd_compounds2$compound = gsub("8-Hydroxy |Cis-| Metabolite| - Total","", apd_compounds2$compound)
apd_compounds2$compound = gsub(".*risperidone","Risperidone", apd_compounds2$compound)
apd_compounds2$compound = gsub(".*clozapine","Clozapine", apd_compounds2$compound)

# label cases with multiple APDs or individual APDs with n<3 as 'other'
apd_compounds2 = apd_compounds2[!duplicated(apd_compounds2),]
apd_compounds_coll = apd_compounds2 %>% group_by(brainnumber)  %>% summarise(tox = paste(compound, collapse = ","))
apd_compounds_coll$tox_simplif = sub(".*,.*","other",apd_compounds_coll$tox)
apd_compounds_coll$tox_simplif = sub("Aripiprazole|Paliperidone|Ziprasidone|Quetiapine","other",apd_compounds_coll$tox_simplif)
apd_compounds_coll$tox_simplif = sub("Chlorpromazine|Perphenazine","other",apd_compounds_coll$tox_simplif)

cmc_metadata_dlpfc_f$APD_compound = apd_compounds_coll$tox_simplif[match(cmc_metadata_dlpfc_f$Brain_ID, apd_compounds_coll$brainnumber)]
cmc_metadata_dlpfc_f$APD_compound[is.na(cmc_metadata_dlpfc_f$APD_compound)] = ""
cmc_metadata_dlpfc_f$DxAPD_compound = paste0(cmc_metadata_dlpfc_f$DxAPD, cmc_metadata_dlpfc_f$APD_compound)

cmc_metadata_dlpfc_f$DxAPD_compound = factor(cmc_metadata_dlpfc_f$DxAPD_compound, levels = c("Control.Neg.", "SCZ.Neg.","SCZ.Pos.Clozapine","SCZ.Pos.Olanzapine","SCZ.Pos.Risperidone","SCZ.Pos.Haloperidol","SCZ.Pos.Fluphenazine","SCZ.Pos.other"))

# Plots for APD compound
df_apdtox_compound = data.frame(SZ_signature = sig_sz_pfc_resid1, group = cmc_metadata_dlpfc_f$DxAPD_compound, Dx.APD=cmc_metadata_dlpfc_f$DxAPD_class)
df_apdtox_compound=droplevels(df_apdtox_compound[df_apdtox_compound$group!="SCZ.Pos.other",])

p_fig1b = ggplot(df_apdtox_compound, aes(x=group, y=SZ_signature, fill=Dx.APD)) + geom_boxplot(outlier.shape = NA) + geom_jitter(position = position_jitter(.1)) + scale_y_continuous() + theme_classic() + theme(axis.text.x = element_text(size=12, angle = 90, hjust = 1, vjust = .5)) + xlab("") + ylab("SCZ expression signature (Z-score)") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1)) + scale_fill_manual(values = RColorBrewer::brewer.pal(8, "Set2")[c(1,2,4,3)])

# Combine plots for Fig. 1
pdf("fig1_scz_signature.pdf", w=10,h=5, useDingbats = F)
plot_grid(p_fig1a, p_fig1b,align = "h", rel_widths = c(1,1.37), labels = c("a","b"), label_size = 18)
dev.off()

# Use Gandal data to define SCZ significant genes and betas for SCZ signature (aggregated DE score)
gandal_dge = read.csv("gandal_dge.csv")
genes_sz_gandal = gandal_dge$ensembl_gene_id[gandal_dge$SCZ.p.value<.05 & gandal_dge$ensembl_gene_id %in% fit_szctrl_pfc_cov1$genes$Gene.stable.ID]
sig_sz_gandal_resid1 = scale(t(resid_szctrl_pfc_cov1$E[match(genes_sz_gandal, fit_szctrl_pfc_cov1$genes$Gene.stable.ID),]) %*% gandal_dge$SCZ.log2FC[match(genes_sz_gandal, gandal_dge$ensembl_gene_id)])

# add Gandal signature to data frame for each individual
df_apdtox_class = data.frame(Dx = cmc_metadata_dlpfc_f$Dx, SZ_signature = sig_sz_pfc_resid1, SZ_signature_gandal = sig_sz_gandal_resid1, group = cmc_metadata_dlpfc_f$DxAPD_class)
df_apdtox_class$group=factor(df_apdtox_class$group)

df_apdtox_compound = data.frame(SZ_signature = sig_sz_pfc_resid1, SZ_signature_gandal = sig_sz_gandal_resid1, group = cmc_metadata_dlpfc_f$DxAPD_compound, Dx.APD=cmc_metadata_dlpfc_f$DxAPD_class)
df_apdtox_compound=droplevels(df_apdtox_compound[df_apdtox_compound$group!="SCZ.Pos.other",])

# Make plots with SCZ signature derived from Gandal (Supplementary Fig. S1)
p_sup_fig1a = ggplot(df_apdtox_class, aes(x=group, y=SZ_signature_gandal, fill=Dx)) + geom_boxplot(outlier.shape = NA) + geom_jitter(position = position_jitter(.1)) + scale_y_continuous() + scale_fill_brewer(palette = "Set2") + theme_classic() + theme(axis.text.x = element_text(size=12, angle = 90, hjust = 1, vjust = .5)) + xlab("") + ylab("SCZ expression signature (Z-score)\nwith PsychENCODE SCZ betas") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1))

p_sup_fig1b =ggplot(df_apdtox_compound, aes(x=group, y=SZ_signature_gandal, fill=Dx.APD)) + geom_boxplot(outlier.shape = NA) + geom_jitter(position = position_jitter(.1)) + scale_y_continuous() + theme_classic() + theme(axis.text.x = element_text(size=12, angle = 90, hjust = 1, vjust = .5)) + xlab("") + ylab("SCZ expression signature (Z-score)\nwith PsychENCODE SCZ betas") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1)) + scale_fill_manual(values = RColorBrewer::brewer.pal(8, "Set2")[c(1,2,4,3)])

pdf("sup_fig_scz_signature_gandal.pdf", w=10,h=5.25, useDingbats = F)
plot_grid(p_sup_fig1a, p_sup_fig1b,align = "h", rel_widths = c(1,1.37), labels = c("a","b"), label_size = 18)
dev.off()


# UCI data DGE and signature validation

# Import Salmon-derived transcript-level estimated counts
dlpfc_uci30_genes_tx=read.table("salmon_30samples_nreads.txt", sep = "|", skip = 1)
colnames(dlpfc_uci30_genes_tx)=c("transcript_id", "gene_id", "Havana_gene_id", "Havana_transcript_id",
                                 "transcript_name", "gene_name", "transcript_length", "transcript_biotype", "rest")
dlpfc_uci30_nreads = read.delim("salmon_30samples_nreads.txt", row.names = 1)
dlpfc_uci30_names=colnames(dlpfc_uci30_nreads)
rownames(dlpfc_uci30_nreads)=dlpfc_uci30_genes_tx$transcript_id

# summarize Salmon transcript-level data at the gene level
dlpfc_uci30_gene_counts = dlpfc_uci30_nreads %>% group_by(dlpfc_uci30_genes_tx$gene_id) %>%summarise_all(sum)
dlpfc_uci30_gene_counts=data.frame(dlpfc_uci30_gene_counts[,-1], row.names = sub("\\..*","",dlpfc_uci30_gene_counts$`dlpfc_uci30_genes_tx$gene_id`))

# Import metadata (including APD toxicology)
metadata_dlpfc_uci30=read.csv("metadata_uci30_apd.csv", stringsAsFactors = T)
table(match(colnames(dlpfc_uci30_gene_counts), paste0("Sample_",metadata_dlpfc_uci30$HSB))==1:30)

# Filter UCI data
gene_data_uci30=data.frame(gene_id=rownames(dlpfc_uci30_gene_counts), gene_name=dlpfc_uci30_genes_tx$gene_name[match(rownames(dlpfc_uci30_gene_counts),sub("\\..*","",dlpfc_uci30_genes_tx$gene_id))])
cpm_uci30=10^6*t(t(dlpfc_uci30_gene_counts)/rowSums(t(dlpfc_uci30_gene_counts)))
idx_cpm1_uci30=rowSums(cpm_uci30>1)>=3 # expression filter (expressed at CPM>1 in at least 3 samples)

idx_uci30_lowRIN = metadata_dlpfc_uci30$RIN.cblm<3 # remove low RIN sample
metadata_dlpfc_uci30_f = droplevels(metadata_dlpfc_uci30[!idx_uci30_lowRIN,])
metadata_dlpfc_uci30_f = metadata_dlpfc_uci30_f[,c("Dx","Gender","Race","Cohort","Age","RIN.cblm","pH","PMI","APD_tox","APD_compound")]
metadata_dlpfc_uci30_f$Dx_APD_tox = paste0(metadata_dlpfc_uci30_f$Dx, ".", metadata_dlpfc_uci30_f$APD_tox)
metadata_dlpfc_uci30_f$Dx_APD_compound = paste0(metadata_dlpfc_uci30_f$Dx, ".", metadata_dlpfc_uci30_f$APD_tox, ".", metadata_dlpfc_uci30_f$APD_compound)
metadata_dlpfc_uci30_f = data.frame(metadata_dlpfc_uci30_f, stringsAsFactors = T)

# Limma-voom for UCI
v_uci30_f =get_voom(counts = dlpfc_uci30_gene_counts[idx_cpm1_uci30,!idx_uci30_lowRIN], genes = gene_data_uci30[idx_cpm1_uci30,])
des30_Dx=model.matrix(~Dx, data = metadata_dlpfc_uci30_f)

# PCA & exploration of covariates
pca_uci30_f = mypca(object = v_uci30_f$E)

pval_uci30_f=matrix(ncol=20, nrow=10)
rownames(pval_uci30_f)=colnames(metadata_dlpfc_uci30_f)[1:10]
colnames(pval_uci30_f)=paste0("PC",1:20)

for (PC in 1:20) {for (varname in colnames(metadata_dlpfc_uci30_f)[1:10]){
  a=summary(lm(pca_uci30_f$x[,PC]~metadata_dlpfc_uci30_f[,varname]))
  pval_uci30_f[varname,PC]=pf(a$fstatistic[1],a$fstatistic[2],a$fstatistic[3],lower.tail=FALSE)[1]
}}
round(-log10(pval_uci30_f))

# UCI model fit, calculate residuals and DGE for SCZ
des_uci30_Dx=model.matrix(~Dx, data = metadata_dlpfc_uci30_f)
des_uci30_cov1=model.matrix(~Gender+RIN.cblm+pH+PMI, data = metadata_dlpfc_uci30_f)

fit_uci30_cov1 = getFit(voom = v_uci30_f, design = des_uci30_cov1)
resid_uci30_cov1 = get_vresid(fit = fit_uci30_cov1, voom = v_uci30_f)

fit_uci30_resid1_Dx = getFit(voom = resid_uci30_cov1, design = des_uci30_Dx)

# explore DE results
table(decideTests(fit_uci30_resid1_Dx[,"DxSCZ"]))
table(decideTests(fit_uci30_resid1_Dx[,"DxSCZ"], adjust.method = "none"))
topTable(fit_uci30_resid1_Dx)

# compare with Gandal and HBCC DGE in SCZ
intergenes_uci30_cmc = intersect(fit_uci30_resid1_Dx$genes$gene_id, fit_szctrl_pfc_cov1$genes$Gene.stable.ID)
intergenes_gandal_cmc = intersect(gandal_dge$ensembl_gene_id, fit_szctrl_pfc_cov1$genes$Gene.stable.ID)
intergenes_gandal_uci30 = intersect(gandal_dge$ensembl_gene_id[!is.na(gandal_dge$SCZ.log2FC)], fit_uci30_resid1_Dx$genes$gene_id)

cor(fit_szctrl_pfc_resid1_Dx$coefficients[match(intergenes_uci30_cmc,fit_szctrl_pfc_resid1_Dx$genes$Gene.stable.ID),"DxSCZ"], fit_uci30_resid1_Dx$coefficients[intergenes_uci30_cmc,"DxSCZ"])
# 0.4208854

cor(fit_szctrl_pfc_resid1_Dx$coefficients[match(intergenes_gandal_cmc,fit_szctrl_pfc_resid1_Dx$genes$Gene.stable.ID),"DxSCZ"], gandal_dge$SCZ.log2FC[match(intergenes_gandal_cmc,gandal_dge$ensembl_gene_id)])
# 0.7214159

cor(fit_uci30_resid1_Dx$coefficients[intergenes_gandal_uci30,"DxSCZ"], gandal_dge$SCZ.log2FC[match(intergenes_gandal_uci30,gandal_dge$ensembl_gene_id)])
# 0.4292515

# get genes with p<.05 for UCI, calculate SCZ (aggregated DE score) signature on UCI data
genes_sz_uci30_resid1 = fit_uci30_resid1_Dx$genes$gene_id[fit_uci30_resid1_Dx$p.value[,"DxSCZ"]<.05]
sig_uci30_resid1 = scale(t(resid_uci30_cov1$E[genes_sz_uci30_resid1,]) %*% fit_uci30_resid1_Dx$coefficients[genes_sz_uci30_resid1,"DxSCZ"])

# Make plots for UCI SCZ signature (Supplementary Fig. S2)
df_apdtox_uci = data.frame(SZ_signature = sig_uci30_resid1, Dx = metadata_dlpfc_uci30_f$Dx,
                           Dx.APD = metadata_dlpfc_uci30_f$Dx_APD_tox, Dx.APD.compound = metadata_dlpfc_uci30_f$Dx_APD_compound, stringsAsFactors = T)

p_uci1 = ggplot(df_apdtox_uci, aes(x=Dx.APD, y=SZ_signature, fill=Dx)) + geom_boxplot(outlier.shape = NA) + geom_jitter(position = position_jitter(.1)) + scale_y_continuous() + scale_fill_brewer(palette = "Set2") + theme_classic() + theme(axis.text.x = element_text(size=12, angle = 90, hjust = 1, vjust = .5)) + xlab("") + ylab("SCZ expression signature (Z-score)\nin UCI validation dataset") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1))

p_uci2 = ggplot(df_apdtox_uci, aes(x=Dx.APD.compound, y=SZ_signature, fill=Dx.APD)) + geom_boxplot(outlier.shape = NA) + geom_jitter(position = position_jitter(.1)) + scale_y_continuous() + theme_classic() + theme(axis.text.x = element_text(size=12, angle = 90, hjust = 1, vjust = .5)) + xlab("") + ylab("SCZ expression signature (Z-score)\nin UCI validation dataset") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1)) + scale_fill_manual(values = RColorBrewer::brewer.pal(8, "Set2")[c(1,2,4)])

# test significance for UCI (parametric)
summary(aov(data=df_apdtox_uci, SZ_signature~Dx.APD)) #ANOVA p=7.12e-06 ***
TukeyHSD(aov(data=df_apdtox_uci, SZ_signature~Dx.APD)) 
# Shapiro-Wilks test for normality: p=0.0004608 (different from normal -> use non-parametric)
shapiro.test(sig_uci30_resid1)
# test significance for UCI (non-parametric)
kruskal.test(data=df_apdtox_uci, SZ_signature~Dx.APD) # KW p=7.78e-05 ****
FSA::dunnTest(data=df_apdtox_uci, SZ_signature~Dx.APD, method="bh")

# Make plots with SCZ sfor UCI (Supplementary Fig. S2)
pdf("sup_fig_uci_sz_signature.pdf", w=7,h=5, useDingbats = F)
plot_grid(p_uci1, p_uci2,align = "h", rel_widths = c(1,1.2), labels = c("a","b"), label_size = 18)
dev.off()

# DGE for subgroups (Supplementary Table S1)
des_szctrl_pfc_DxAPD_class = model.matrix(data = data.frame(metadata_szctrl_pfc, grp=cmc_metadata_dlpfc_f$DxAPD_class), ~grp)
fit_szctrl_pfc_resid1_DxAPD_class = getFit(voom = resid_szctrl_pfc_cov1, design = des_szctrl_pfc_DxAPD_class)

colnames(fit_szctrl_pfc_resid1_DxAPD_class)

tt_HBCC.SCZ.Neg = topTable(fit_szctrl_pfc_resid1_DxAPD_class, coef = 2, sort.by = "none", n=Inf)
tt_HBCC.SCZ.Pos.atypical = topTable(fit_szctrl_pfc_resid1_DxAPD_class, coef = 3, sort.by = "none", n=Inf)
tt_HBCC.SCZ.Pos.mixed = topTable(fit_szctrl_pfc_resid1_DxAPD_class, coef = 4, sort.by = "none", n=Inf)
tt_HBCC.SCZ.Pos.typical = topTable(fit_szctrl_pfc_resid1_DxAPD_class, coef = 5, sort.by = "none", n=Inf)

# check number of DE genes (highest for atypical APDs)
table(decideTests(fit_szctrl_pfc_resid1_DxAPD_class, adjust.method = "BH")[,2])
table(decideTests(fit_szctrl_pfc_resid1_DxAPD_class, adjust.method = "BH")[,3])
table(decideTests(fit_szctrl_pfc_resid1_DxAPD_class, adjust.method = "BH")[,4])
table(decideTests(fit_szctrl_pfc_resid1_DxAPD_class, adjust.method = "BH")[,5])

l_subgrp_p = list(SCZ.Neg = tt_HBCC.SCZ.Neg[,c(1,2,7,9:11)] %>% filter(P.Value<.05),
                  SCZ.Pos.atypical = tt_HBCC.SCZ.Pos.atypical[,c(1,2,7,9:11)] %>% filter(P.Value<.05),
                  SCZ.Pos.mixed = tt_HBCC.SCZ.Pos.mixed[,c(1,2,7,9:11)] %>% filter(P.Value<.05),
                  SCZ.Pos.typical = tt_HBCC.SCZ.Pos.typical[,c(1,2,7,9:11)] %>% filter(P.Value<.05))

openxlsx::write.xlsx(l_subgrp_p, file="tab_dge_scz_subgrp.xlsx")

save.image("apd_scz_toxicology.rdata")