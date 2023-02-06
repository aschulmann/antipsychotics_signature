
library(dplyr)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(ggrepel)
library(ComplexHeatmap)
library(edgeR)
library(WGCNA)
library(GEOquery)
library(affy)

# Macaque DGE

# Load human data
load("apd_scz_toxicology.rdata")
# Get macaque data
cmc_macaque_countdata=read.delim("CMC_RhesusMacaque_DLPFC_geneExpressionRaw.txt", row.names = 1)
cmc_macaque_metadata=read.csv("CMC_RhesusMacaque_Clinical_DLPFCmRNA-metaData.csv")

# Filter & reformat macaque data & metadata
idx_mmul_outlier # 1 outlier sample ro be filtered (see below)
cmc_macaque_metadata = droplevels(cmc_macaque_metadata[!idx_mmul_outlier,])

cmc_macaque_metadata$DLPFC_RNA_isolation_Batch = factor(paste0("B",cmc_macaque_metadata$DLPFC_RNA_isolation_Batch))
cmc_macaque_metadata$DLPFC_RNA_Sequencing_Library_Batch = factor(paste0("B",cmc_macaque_metadata$DLPFC_RNA_Sequencing_Library_Batch))
cmc_macaque_metadata$DLPFC_RNA_Sequencing_Flowcell_Batch = factor(paste0("B",cmc_macaque_metadata$DLPFC_RNA_Sequencing_Flowcell_Batch))
cmc_macaque_metadata$DLPFC_RNA_Sequencing_Ribozero_Batch = factor(paste0("B",cmc_macaque_metadata$DLPFC_RNA_Sequencing_Ribozero_Batch))
cmc_macaque_metadata$sex = factor(cmc_macaque_metadata$sex, levels = c("male", "female"))

cmc_macaque_metadata$group=as.character(cmc_macaque_metadata$treatmentType)
cmc_macaque_metadata$group[which(cmc_macaque_metadata$dose..mg.kg.d.==0.14)]="haloperidol_low"
cmc_macaque_metadata$group[which(cmc_macaque_metadata$dose..mg.kg.d.==4)]="haloperidol_high"
cmc_macaque_metadata$group=factor(cmc_macaque_metadata$group, levels = unique(cmc_macaque_metadata$group))
cmc_macaque_metadata=cmc_macaque_metadata[match(colnames(cmc_macaque_countdata)[!idx_mmul_outlier],cmc_macaque_metadata$DLPFC_RNA_Sequencing_Sample_ID),]

mmul_metadata_simple = data.frame(group = cmc_macaque_metadata$group,
                                  CLZ_dose = ifelse(cmc_macaque_metadata$treatmentType == "clozapine", cmc_macaque_metadata$dose..mg.kg.d., 0),
                                  HAL_dose = ifelse(cmc_macaque_metadata$treatmentType == "haloperidol", cmc_macaque_metadata$dose..mg.kg.d., 0),
                                  Sex = cmc_macaque_metadata$sex,
                                  batch_rna = cmc_macaque_metadata$DLPFC_RNA_isolation_Batch,
                                  batch_ribo0 = cmc_macaque_metadata$DLPFC_RNA_Sequencing_Ribozero_Batch,
                                  batch_lib = cmc_macaque_metadata$DLPFC_RNA_Sequencing_Library_Batch,
                                  batch_seq = cmc_macaque_metadata$DLPFC_RNA_Sequencing_Flowcell_Batch,
                                  Age = cmc_macaque_metadata$Age.of.Death,
                                  RIN = cmc_macaque_metadata$DLPFC_RNA_isolation_RIN,
                                  total_RNA = cmc_macaque_metadata$DLPFC_RNA_isolation_TotalYield_ug,
                                  total_reads = cmc_macaque_metadata$DLPFC_RNA_Sequencing_Total_Reads,
                                  total_genes = cmc_macaque_metadata$DLPFC_RNA_Sequencing_Genes_Detected,
                                  rate_mapped = cmc_macaque_metadata$DLPFC_RNA_Sequencing_Percent_Aligned,
                                  rate_intron = cmc_macaque_metadata$DLPFC_RNA_Sequencing_Intronic_Rate,
                                  rate_intergenic = cmc_macaque_metadata$DLPFC_RNA_Sequencing_Intergenic_Rate,
                                  rate_rRNA = cmc_macaque_metadata$DLPFC_RNA_Sequencing_rRNA_Rate,
                                  rate_efficiency = cmc_macaque_metadata$DLPFC_RNA_Sequencing_Expression_Profiling_Efficiency)

cmc_macaque_cpm=10^6*t(t(cmc_macaque_countdata)/rowSums(t(cmc_macaque_countdata)))
idx_mmul_cpm1=rowSums(cmc_macaque_cpm[,!idx_mmul_outlier]>1)>=3 # expression filter (expressed at CPM>1 in at least 3 samples)

# Limma/voom for macaque
v_mmul=get_voom(counts = cmc_macaque_countdata[idx_mmul_cpm1,!idx_mmul_outlier], genes=NULL)

# Run PCA for macaque
pca_mmul = mypca(object = v_mmul$E)

## check outliers
#a=apply(abs(scale(pca_mmul$x[,1:20])), 2, which.max)
#unique(a[duplicated(a)])
#apply(abs(scale(pca_mmul$x[,1:20])), 2, max)

#pc_mmul_2sd = abs(scale(pca_mmul$x[,1:20]))>2
#pc_mmul_3sd = abs(scale(pca_mmul$x[,1:20]))>3

#hist(rowSums(pc_mmul_2sd))
#hist(rowSums(pc_mmul_3sd))
#hist(rowSums(pc_mmul_5sd))

#rownames(pc_mmul_3sd)[rowSums(pc_mmul_3sd)>=2]

## remove 1 outlier based on abs(Z)>3 in 2 or more of the top 20 PCs (also highest Z)
## Go back to filter and define idx_mmul_outlier to filter this sample

# Explore PCA & covariates
pval_pc_mmul=matrix(ncol=30, nrow=18)
rownames(pval_pc_mmul)=colnames(mmul_metadata_simple)
colnames(pval_pc_mmul)=paste0("PC",1:30)

for (PC in 1:30) {for (varname in colnames(mmul_metadata_simple)){
  a=summary(lm(pca_mmul$x[,PC]~mmul_metadata_simple[,varname]))
  pval_pc_mmul[varname,PC]=pf(a$fstatistic[1],a$fstatistic[2],a$fstatistic[3],lower.tail=FALSE)[1]
}}
round(-log10(pval_pc_mmul), digits = 1)

# Import gene models and human orthologs

mmul_anno_grch37=read.delim("mmul_anno_ens75.tsv")
mmul_anno_grch37=read.delim("ens75_mmul2hu.tsv")
table(rownames(cmc_macaque_countdata) %in% mmul_anno_grch37$Ensembl.Gene.ID)

gene_data_mmul_f=mmul_anno_grch37[match(rownames(cmc_macaque_countdata)[idx_mmul_cpm1],mmul_anno_grch37$Ensembl.Gene.ID),]
gene_data_mmul=mmul_anno_grch37[match(rownames(cmc_macaque_countdata),mmul_anno_grch37$Ensembl.Gene.ID),]

# Limma/voom model fit and DGE for macaques
v_mmul=get_voom(counts = cmc_macaque_countdata[idx_mmul_cpm1,!idx_mmul_outlier], genes=gene_data_mmul_f) # add genes to voom

des_mmul_cov1 = model.matrix(~Sex + batch_rna + RIN + rate_intergenic + rate_efficiency , data = mmul_metadata_simple)
fit_mmul_cov1 = getFit(v_mmul, design = des_mmul_cov1)
resid_mmul_cov1 = get_vresid(fit = fit_mmul_cov1, voom = v_mmul)

des_group = model.matrix(data=mmul_metadata_simple, ~group)
fit_mmul_resid1_group = getFit(voom = resid_mmul_cov1, design = des_group)

tt_mmul_grp_CLZ = topTable(fit_mmul_resid1_group, coef = 2, sort.by = "none", confint = T, n=Inf)
tt_mmul_grp_HAL.lo = topTable(fit_mmul_resid1_group, coef = 3, sort.by = "none", confint = T, n=Inf)
tt_mmul_grp_HAL.hi = topTable(fit_mmul_resid1_group, coef = 4, sort.by = "none", confint = T, n=Inf)

# GO term enrichment

# BioMart annotations for GRCh37 mapped to macaque
go_grch37 = read.delim("go_grch37.txt")
go_grch37 = go_grch37[go_grch37$GO.term.name!="",]
go_grch37 = split(x = go_grch37$Gene.stable.ID, f = go_grch37$GO.term.name)

# BioMart annotation for Ensembl v86 for human
go_hu_ens86 = read.delim("go_hu_ens86.txt")
go_hu_ens86 = go_hu_ens86[go_hu_ens86$GO.Term.Name!="",]
go_hu_ens86 = split(x = go_hu_ens86$Ensembl.Gene.ID, f = go_hu_ens86$GO.Term.Name)

# Define function testing Fisher's Exact test for over-representation
overrep.test=function(geneList,geneUniverse,geneNames,go2genes){
  et=list()
  for (i in 1:length(go2genes)){
    gl_hits=intersect(go2genes[[i]],geneList)
    bg_hits=intersect(go2genes[[i]],geneUniverse)
    hit_names = paste(geneNames[match(gl_hits, geneUniverse)], collapse = ",")
    et[[i]]=fisher.test(matrix(c(length(gl_hits),
                                 length(bg_hits),
                                 length(geneList),
                                 length(geneUniverse)),
                               nrow = 2,ncol = 2))
    et[[i]]=c(length(gl_hits), hit_names, et[[i]]$estimate,et[[i]]$p.value)
    names(et[[i]])=c( "n_hits", "gene_name","odds","pval")
  }
  et=data.frame(terms=names(go2genes), t(simplify2array(et)))
  et$pval=as.numeric(et$pval)
  et$odds=as.numeric(et$odds)
  et=et[order(et$pval),]
  et=droplevels(et[et$odds>1,])
  rownames(et)=NULL
  return(et)
}

human_HBCC_SCZ = overrep.test(geneList = tt_HBCC.SCZ$Gene.stable.ID[tt_HBCC.SCZ$P.Value<.05],
                              geneUniverse = tt_HBCC.SCZ$Gene.stable.ID, geneNames = tt_HBCC.SCZ$Gene.name,
                              go2genes = go_hu_ens86)
monkey_apd_go = list()
monkey_apd_go[["CLZ"]] = overrep.test(geneList = tt_mmul_grp_CLZ$Human.Ensembl.Gene.ID[tt_mmul_grp_CLZ$P.Value<.05],
                                      geneUniverse = tt_mmul_grp_CLZ$Human.Ensembl.Gene.ID, geneNames = tt_mmul_grp_CLZ$Associated.Gene.Name,
                                      go2genes = go_grch37)
monkey_apd_go[["HAL_lo"]] = overrep.test(geneList = tt_mmul_grp_HAL.lo$Human.Ensembl.Gene.ID[tt_mmul_grp_HAL.lo$P.Value<.05],
                                         geneUniverse = tt_mmul_grp_HAL.lo$Human.Ensembl.Gene.ID, geneNames = tt_mmul_grp_HAL.lo$Associated.Gene.Name,
                                         go2genes = go_grch37)
monkey_apd_go[["HAL_hi"]] = overrep.test(geneList = tt_mmul_grp_HAL.hi$Human.Ensembl.Gene.ID[tt_mmul_grp_HAL.hi$P.Value<.05],
                                         geneUniverse = tt_mmul_grp_HAL.hi$Human.Ensembl.Gene.ID, geneNames = tt_mmul_grp_HAL.hi$Associated.Gene.Name,
                                         go2genes = go_grch37)


human_HBCC_SCZ_dn = overrep.test(geneList = tt_HBCC.SCZ$Gene.stable.ID[head(order(tt_HBCC.SCZ$t),n=200)],
                                 geneUniverse = tt_HBCC.SCZ$Gene.stable.ID, geneNames = tt_HBCC.SCZ$Gene.name,
                                 go2genes = go_hu_ens86)
monkey_apd_go_dn = list()
monkey_apd_go_dn[["CLZ"]] = overrep.test(geneList = tt_mmul_grp_CLZ$Human.Ensembl.Gene.ID[head(order(tt_mmul_grp_CLZ$t),n=200)],
                                         geneUniverse = tt_mmul_grp_CLZ$Human.Ensembl.Gene.ID, geneNames = tt_mmul_grp_CLZ$Associated.Gene.Name,
                                         go2genes = go_grch37)
monkey_apd_go_dn[["HAL_lo"]] = overrep.test(geneList = tt_mmul_grp_HAL.lo$Human.Ensembl.Gene.ID[head(order(tt_mmul_grp_HAL.lo$t),n=200)],
                                            geneUniverse = tt_mmul_grp_HAL.lo$Human.Ensembl.Gene.ID, geneNames = tt_mmul_grp_HAL.lo$Associated.Gene.Name,
                                            go2genes = go_grch37)
monkey_apd_go_dn[["HAL_hi"]] = overrep.test(geneList = tt_mmul_grp_HAL.hi$Human.Ensembl.Gene.ID[head(order(tt_mmul_grp_HAL.hi$t),n=200)],
                                            geneUniverse = tt_mmul_grp_HAL.hi$Human.Ensembl.Gene.ID, geneNames = tt_mmul_grp_HAL.hi$Associated.Gene.Name,
                                            go2genes = go_grch37)

human_HBCC_SCZ_up = overrep.test(geneList = tt_HBCC.SCZ$Gene.stable.ID[tail(order(tt_HBCC.SCZ$t),n=200)],
                                 geneUniverse = tt_HBCC.SCZ$Gene.stable.ID, geneNames = tt_HBCC.SCZ$Gene.name,
                                 go2genes = go_hu_ens86)
monkey_apd_go_up = list()
monkey_apd_go_up[["CLZ"]] = overrep.test(geneList = tt_mmul_grp_CLZ$Human.Ensembl.Gene.ID[tail(order(tt_mmul_grp_CLZ$t),n=200)],
                                         geneUniverse = tt_mmul_grp_CLZ$Human.Ensembl.Gene.ID, geneNames = tt_mmul_grp_CLZ$Associated.Gene.Name,
                                         go2genes = go_grch37)
monkey_apd_go_up[["HAL_lo"]] = overrep.test(geneList = tt_mmul_grp_HAL.lo$Human.Ensembl.Gene.ID[tail(order(tt_mmul_grp_HAL.lo$t),n=200)],
                                            geneUniverse = tt_mmul_grp_HAL.lo$Human.Ensembl.Gene.ID, geneNames = tt_mmul_grp_HAL.lo$Associated.Gene.Name,
                                            go2genes = go_grch37)
monkey_apd_go_up[["HAL_hi"]] = overrep.test(geneList = tt_mmul_grp_HAL.hi$Human.Ensembl.Gene.ID[tail(order(tt_mmul_grp_HAL.hi$t),n=200)],
                                            geneUniverse = tt_mmul_grp_HAL.hi$Human.Ensembl.Gene.ID, geneNames = tt_mmul_grp_HAL.hi$Associated.Gene.Name,
                                            go2genes = go_grch37)

# Write out results for GO enrichment (Supplementary Table S3)
tab_dge_go = list(CLZ_DGE = filter(tt_mmul_grp_CLZ, P.Value<.05) %>% arrange(P.Value), CLZ_GO_up = monkey_apd_go_up$CLZ[1:50,], CLZ_GO_dn = monkey_apd_go_dn$CLZ[1:50,],
                  HAL.lo_DGE = filter(tt_mmul_grp_HAL.lo, P.Value<.05) %>% arrange(P.Value), HAL.lo_GO_up = monkey_apd_go_up$HAL_lo[1:50,], HAL.lo_GO_dn = monkey_apd_go_dn$HAL_lo[1:50,],
                  HAL.hi_DGE = filter(tt_mmul_grp_HAL.hi, P.Value<.05) %>% arrange(P.Value), HAL.hi_GO_up = monkey_apd_go_up$HAL_hi[1:50,], HAL.hi_GO_dn = monkey_apd_go_dn$HAL_hi[1:50,],
                  SCZ_DGE = filter(tt_HBCC.SCZ, P.Value<.05) %>% arrange(P.Value), SCZ_GO_up = human_HBCC_SCZ_up[1:50,], SCZ_GO_dn = human_HBCC_SCZ_dn[1:50,])

openxlsx::write.xlsx(tab_dge_go, file="tab_dge_go.xlsx")

# Plots for GO term enrichment (Supplementary Fig. 3)

GO_plot = function(df){
  df$terms = factor(df$terms, levels = rev(df$terms))
  df$n_hits = as.numeric(df$n_hits)
  require(ggplot2)
  ggplot(data = df, aes(y=terms, x=-log10(pval), size=n_hits, color=log(odds))) + geom_point() + theme_bw() + scale_color_gradientn(colours = RColorBrewer::brewer.pal(9, "YlOrRd"), limits = c(0,4.5)) + scale_size(range = c(min(df$n_hits)*0.5,max(df$n_hits)*0.5))
}

pgo_CLZ_up = GO_plot(df = head(monkey_apd_go_up$CLZ, n=20)) + ggtitle("CLZ (upgregulated)")
pgo_CLZ_dn = GO_plot(df = head(monkey_apd_go_dn$CLZ, n=20)) + ggtitle("CLZ (downgregulated)")

pgo_HAL.lo_up = GO_plot(df = head(monkey_apd_go_up$HAL_lo, n=20)) + ggtitle("HAL.lo (upgregulated)")
pgo_HAL.lo_dn = GO_plot(df = head(monkey_apd_go_dn$HAL_lo, n=20)) + ggtitle("HAL.lo (downgregulated)")

pgo_HAL.hi_up = GO_plot(df = head(monkey_apd_go_up$HAL_hi, n=20)) + ggtitle("HAL.hi (upgregulated)")
pgo_HAL.hi_dn = GO_plot(df = head(monkey_apd_go_dn$HAL_hi, n=20)) + ggtitle("HAL.hi (downgregulated)")

pgo_SCZ_up = GO_plot(df = head(human_HBCC_SCZ_up, n=20)) + ggtitle("SCZ (upgregulated)")
pgo_SCZ_dn = GO_plot(df = head(human_HBCC_SCZ_dn, n=20)) + ggtitle("SCZ (downgregulated)")

# Prepare macaque and human data for WGCNA
# Restrict limma/voom to 1-to-1 orthologs between Macaque and Human
# Get residuals with technical (but not demographic) variables regressed out

idx_1to1 = gene_data_mmul$Homology.Type == "ortholog_one2one"
table(idx_1to1)
v_mmul_1to1=get_voom(counts = cmc_macaque_countdata[idx_mmul_cpm1&idx_1to1,!idx_mmul_outlier], genes = gene_data_mmul[idx_mmul_cpm1&idx_1to1,])

des_mmul_cov2 = model.matrix(~ batch_rna + RIN + rate_intergenic + rate_efficiency , data = mmul_metadata_simple)
fit_mmul_cov2_1to1 = getFit(v_mmul_1to1, design = des_mmul_cov2)
resid_mmul_cov2_1to1 = get_vresid(fit = fit_mmul_cov2_1to1, voom = v_mmul_1to1)

des_szctrl_pfc_cov2 = model.matrix(data=metadata_szctrl_pfc,~batch_lib+RIN+rate_efficiency+rate_intergenic+PMI+pH)
fit_szctrl_pfc_cov2 = getFit(voom = v_pfc_f, design = des_szctrl_pfc_cov2)
resid_szctrl_pfc_cov2 = get_vresid(fit = fit_szctrl_pfc_cov2, voom = v_pfc_f)

## create reference table for human and macaque joint analysis
intergenes_mmul_human = intersect(fit_mmul_cov2_1to1$genes$Human.Ensembl.Gene.ID, fit_szctrl_pfc_cov2$genes$Gene.stable.ID)

idx_mmul_cov1.uniq = match(intergenes_mmul_human, fit_mmul_cov1$genes$Human.Ensembl.Gene.ID)
df_pfc_human_mmul = data.frame(fit_mmul_cov1$genes[idx_mmul_cov1.uniq,1:5],
                               CLZ = tt_mmul_grp_CLZ[idx_mmul_cov1.uniq,c(9,13,14)],
                               HAL.lo = tt_mmul_grp_HAL.lo[idx_mmul_cov1.uniq,c(9,13,14)],
                               HAL.hi = tt_mmul_grp_HAL.hi[idx_mmul_cov1.uniq,c(9,13,14)])
table(is.na(df_pfc_human_mmul)) # make sure there are no NAs

idx_human_cov1.uniq = match(df_pfc_human_mmul$Human.Ensembl.Gene.ID, fit_szctrl_pfc_cov1$genes$Gene.stable.ID)
df_pfc_human_mmul = data.frame(df_pfc_human_mmul, HBCC.SCZ = tt_HBCC.SCZ[idx_human_cov1.uniq,c(7,9,10)])


# Run WGCNA (needs more RAM)
library(WGCNA)
options(stringsAsFactors = FALSE)
enableWGCNAThreads()

idx_human_cov2.uniq = match(df_pfc_human_mmul$Human.Ensembl.Gene.ID, fit_szctrl_pfc_cov2$genes$Gene.stable.ID)
idx_mmul_cov2.uniq = match(df_pfc_human_mmul$Human.Ensembl.Gene.ID, fit_mmul_cov2_1to1$genes$Human.Ensembl.Gene.ID)

expr_human = t(resid_szctrl_pfc_cov2$E[idx_human_cov2.uniq,])
expr_mmul = t(resid_mmul_cov2_1to1$E[idx_mmul_cov2.uniq,])
colnames(expr_mmul)=resid_mmul_cov2_1to1$genes$Human.Ensembl.Gene.ID[idx_mmul_cov2.uniq]
colnames(expr_human)=resid_szctrl_pfc_cov2$genes$Gene.stable.ID[idx_human_cov2.uniq]
table(colnames(expr_mmul) == colnames(expr_human)) # check that gene names are the same

multiExpr = multiData(human = expr_human, monkey = expr_mmul)

# Explore soft power (code from WGCNA Vignette)
# Choose a set of soft-thresholding powers
powers = c(seq(1,10,by=1), seq(12,20, by=2));
# Initialize a list to hold the results of scale-free analysis
powerTables = vector(mode = "list", length = 2);
# Call the network topology analysis function for each set in turn
for (set in 1:2){
  powerTables[[set]] = list(data = pickSoftThreshold(multiExpr[[set]]$data, powerVector=powers,
                                                     verbose = 2, networkType = "signed")[[2]]);
}
par(mfrow=c(2,2))
plot(powerTables[[1]]$data$Power,powerTables[[1]]$data$SFT.R.sq)
plot(powerTables[[1]]$data$Power,powerTables[[1]]$data$mean.k.)
plot(powerTables[[2]]$data$Power,powerTables[[2]]$data$SFT.R.sq)
plot(powerTables[[2]]$data$Power,powerTables[[2]]$data$mean.k.)

# Use sft of 12
# One-step command for consensus WGCNA

net = blockwiseConsensusModules(
  multiExpr, power = 12, networkType = "signed",
  maxBlockSize = 15000, # larger than dataset to avoid chunks
  minModuleSize = 30, deepSplit = 2,
  pamRespectsDendro = FALSE, 
  mergeCutHeight = 0.25, numericLabels = TRUE,
  minKMEtoStay = 0.3,
  saveTOMs = FALSE, verbose = 5)

# Extract module eigengenes
consMEs = net$multiMEs
moduleLabels = net$colors
# Convert the numeric labels to color labels
moduleColors = labels2colors(moduleLabels)
consTree = net$dendrograms[[1]]

# Plot WGCNA dendrogram
plotDendroAndColors(consTree, moduleColors,
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Human & Monkey Consensus gene dendrogram and module colors")

# Create lookup table
table(moduleLabels) # 37 modules
df_wgcna = data.frame(module = net$colors, moduleCol = moduleColors, df_pfc_human_mmul)
df_wgcna = df_wgcna[df_wgcna$module!=0,]
df_wgcna$module = factor(df_wgcna$module, levels = 1:37)

cols_human_mmul = as.character(df_wgcna$moduleCol[match(levels(df_wgcna$module),df_wgcna$module)])
df_wgcna$moduleCol = factor(df_wgcna$moduleCol, levels = cols_human_mmul)

# Module eigengene linear regression for each macaque APD and SCZ
lm_MEs = matrix(ncol=4,nrow=37)
lmt_MEs = matrix(ncol=4,nrow=37)
pval_MEs = matrix(ncol=4,nrow=37)
colnames(lm_MEs) = colnames(pval_MEs) = colnames(lmt_MEs) = c("CLZ","HAL.lo","HAL.hi", "SCZ")
rownames(lm_MEs) = rownames(pval_MEs) = rownames(lmt_MEs) = paste0("M",1:37)

for(i in 1:37){
  a=as.data.frame(consMEs$monkey)[,paste0("data.ME",i)]
  a=summary(lm(data = data.frame(mmul_metadata_simple,ME=a), ME~group))
  b=as.data.frame(consMEs$human)[,paste0("data.ME",i)]
  b=summary(lm(data = data.frame(metadata_szctrl_pfc, ME=b), ME~Dx))
  lm_MEs[i,] = c(a$coefficients[2:4,1], b$coefficients[2,1])
  lmt_MEs[i,] = c(a$coefficients[2:4,3], b$coefficients[2,3])
  pval_MEs[i,] = c(a$coefficients[2:4,4], b$coefficients[2,4])
}

# Make significance labels (APD)
signif_MEs = symnum(pval_MEs, corr = FALSE, na = FALSE,
                     cutpoints = c(0,0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", ".", " "), legend = F)

# Get prioritized genes from PGC3 SCZ paper (Trubetskoy et al, 2022)
pgc3_prioritized = read.csv("PGC3_SCZ_S12_prioritized.csv")
pgc3_scz_ens = intersect(pgc3_prioritized$Ensembl.ID, df_wgcna$Human.Ensembl.Gene.ID)

fisher_pgc = matrix(ncol=4,nrow=37)
for(i in 1:37){
  a=df_wgcna[df_wgcna$module==i,]
  b=df_wgcna
  f=fisher.test(matrix(c(sum(a$Human.Ensembl.Gene.ID %in% pgc3_scz_ens),
                         sum(b$Human.Ensembl.Gene.ID %in% pgc3_scz_ens),
                         nrow(a),
                         nrow(b)),
                       nrow = 2,ncol = 2))
  gnames = intersect(a$Human.Ensembl.Gene.ID, pgc3_scz_ens)
  gnames = pgc3_prioritized$Symbol.ID[match(gnames,pgc3_prioritized$Ensembl.ID)]
  gnames = paste(gnames, collapse = ",")
  fisher_pgc[i,] = c(paste0("M",i), gnames, f$estimate, f$p.value)
}

colnames(fisher_pgc) = c("Module", "Gene hits", "Odds Ratio", "P")

# Make significance labels (PGC3)
signif_pgc = symnum(x = as.numeric(fisher_pgc[,4]), corr = FALSE, na = FALSE,
                     cutpoints = c(0,0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", ".", " "), legend = F)

# GO term enrichment for WGCNA (Supplementary Table S4)
l_go_wgcna = list()
for(i in 1:37){
  print(paste("Working on Module", i))
  l_go_wgcna[[i]] = overrep.test(geneList = df_wgcna$Human.Ensembl.Gene.ID[df_wgcna$module==i],
                                 geneUniverse = df_wgcna$Human.Ensembl.Gene.ID, geneNames = df_wgcna$Associated.Gene.Name,
                                 go2genes = go_grch37)
}

# GO enrichment with more than 2 hits
go_wgcna2 = sapply(l_go_wgcna, function(x){
  x=x[as.numeric(x[,2])>2,]
  return(x[1,1])
})

# GO enrichment with more than 3 hits
go_wgcna3 = sapply(l_go_wgcna, function(x){
  x=x[as.numeric(x[,2])>3,]
  return(x[1,1])
})


# Plot overview of consensus WGCNA results (Fig. 2a)
# Module eigengene stats, GWAS enrichment, module size, top GO terms

moduleSize = as.data.frame(table(moduleLabels))$Freq[-1]

ha_go = rowAnnotation("# genes\nin module" = anno_barplot(moduleSize), GO = anno_text(go_wgcna3, gp = gpar(fontsize = 10)), annotation_name_gp= gpar(fontsize = 9))
ha_pgc = Heatmap(as.matrix(data.frame(SCZ.GWAS=log10(as.numeric(fisher_pgc[,"P"])))), col=circlize::colorRamp2(seq(-4,0, length = 9), rev(RColorBrewer::brewer.pal(9, "Greens"))), 
                 column_names_rot = 45, rect_gp = gpar(col = 1, lwd = 1), heatmap_legend_param = list(color_bar="continuous", title="GWAS\nenrichment"),
                 cell_fun = function(j,i,x, y, width, height, fill) {
                   grid.text(as.matrix(data.frame(signif_pgc))[i,j], x, y, gp = gpar(fontsize = 10))})
ha_pgc

ha_cols = HeatmapAnnotation(text = anno_text(rownames(lmt_MEs), just = "left", gp = gpar(col=ifelse(rownames(lmt_MEs)%in%example_modules,"red","black")) ), which = "row")

set.seed(1000)

h_fig2a = Heatmap(lmt_MEs, col=circlize::colorRamp2(seq(-2,2, length = 7), rev(RColorBrewer::brewer.pal(7, "RdBu"))), rect_gp = gpar(col = 1, lwd = 1), column_names_gp = gpar(fontsize = 12),
                  heatmap_legend_param = list(color_bar="continuous", title="t-statistic"), column_title_gp = gpar(fontsize = 7), row_names_gp = gpar(fontsize = 12), row_names_side = "left",
                  cluster_rows = T, cluster_columns = F, row_km = 4, column_names_rot = 45, column_split = factor(c(rep("Monkey",3),"Human"), levels = c("Monkey","Human")),
                  left_annotation = ha_cols, show_row_names = F,
                  cell_fun = function(j, i, x, y, width, height, fill) {
                    grid.text(sprintf("%s", as.character(signif_MEs[i,j])), x, y, gp = gpar(fontsize = 10))}) + ha_pgc + ha_go 

p_fig2a = grid.grabExpr(draw(h_fig2a)) 

# Select genes with low P-values to highlight
idx_signif_CLZ = df_wgcna$HBCC.SCZ.P.Value<.05 & df_wgcna$CLZ.P.Value<.25
idx_signif_HAL.hi = df_wgcna$HBCC.SCZ.P.Value<.05 & df_wgcna$HAL.hi.P.Value<.25
idx_signif_HAL.lo = df_wgcna$HBCC.SCZ.P.Value<.05 &  df_wgcna$HAL.lo.P.Value<.25

# Create table for all genes
df_apd_wgcna = data.table::rbindlist(list(CLZ = data.frame(df_wgcna[,c(colnames(df_wgcna)[1:5],"HBCC.SCZ.t")], t_statistic=df_wgcna$CLZ.t),
                                          HAL.lo = data.frame(df_wgcna[,c(colnames(df_wgcna)[1:5],"HBCC.SCZ.t")], t_statistic=df_wgcna$HAL.lo.t),
                                          HAL.hi = data.frame(df_wgcna[,c(colnames(df_wgcna)[1:5],"HBCC.SCZ.t")], t_statistic=df_wgcna$HAL.hi.t)), idcol = "APD")
df_apd_wgcna$APD = factor(df_apd_wgcna$APD, levels = c("CLZ","HAL.lo","HAL.hi"))
df_apd_wgcna$module = factor(paste0("M",df_apd_wgcna$module), levels = paste0("M",1:37))

df_apd_wgcna_signif = data.table::rbindlist(list(CLZ = data.frame(df_wgcna[idx_signif_CLZ,c(colnames(df_wgcna)[1:5],"HBCC.SCZ.t")], t_statistic=df_wgcna$CLZ.t[idx_signif_CLZ]),
                                                 HAL.lo = data.frame(df_wgcna[idx_signif_HAL.lo,c(colnames(df_wgcna)[1:5],"HBCC.SCZ.t")], t_statistic=df_wgcna$HAL.lo.t[idx_signif_HAL.lo]),
                                                 HAL.hi = data.frame(df_wgcna[idx_signif_HAL.hi,c(colnames(df_wgcna)[1:5],"HBCC.SCZ.t")], t_statistic=df_wgcna$HAL.hi.t[idx_signif_HAL.hi])), idcol = "APD")
df_apd_wgcna_signif$APD = factor(df_apd_wgcna_signif$APD, levels = c("CLZ","HAL.lo","HAL.hi"))
df_apd_wgcna_signif$module = factor(paste0("M",df_apd_wgcna_signif$module), levels = paste0("M",1:37))

# Plot 5 example modules (Fig. 2b)
example_modules = paste0("M",c(16,9,11,25,31))
df_wgcna_example = df_apd_wgcna[df_apd_wgcna$module %in% example_modules,]
df_wgcna_example$module = factor(df_wgcna_example$module, levels = example_modules)
df_wgcna_example_signif = droplevels(df_apd_wgcna_signif[df_apd_wgcna_signif$module %in% example_modules,])
df_wgcna_example_signif$module = factor(df_wgcna_example_signif$module, levels = example_modules)

p_fig2b=ggplot(data = droplevels(df_wgcna_example), aes(x = t_statistic, y = HBCC.SCZ.t, label = Associated.Gene.Name)) + geom_point(color="gray") + theme_bw() + 
  geom_vline(xintercept=0, linetype="dotted") + geom_hline(yintercept=0, linetype="dotted") + facet_grid(vars(module), vars(APD)) + ylab("SCZ t-statistic (Human)") + xlab("APD t-statistic (Monkey)") +
  theme(legend.position="none") + geom_point(data = df_wgcna_example_signif, color="red") +
  geom_text_repel(data = df_wgcna_example_signif, max.overlaps = 20, size = 1.75, point.size=1)

pdf("fig2_wgcna_monkey_human.pdf", w=12.75,h=9.5, useDingbats = F)
plot_grid(p_fig2a, p_fig2b, rel_widths = c(1.05,1), labels = c("a","b"), label_size = 18)
dev.off()


# Correlate betas for macaque with mouse (Rizig et al. 2012)

# Obtain data from GEO 
library(GEOquery)
library(limma)
library(affy)

gset = getGEO("GSE6467", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx = grep("GPL339", attr(gset, "names")) else idx = 1
gset = gset[[idx]]

# Correct gene labels 
fvarLabels(gset) = make.names(fvarLabels(gset))

# Get expression values
ex = exprs(gset)

# Create metadata for mouse
metadata_mouse = data.frame(sample_name=gset$title, group=sub(" .*","",gset$title))
metadata_mouse$group = factor(metadata_mouse$group, levels = c("Untreated","Clozapine","Haloperidol"))

# DGE for mouse
des_mouse_group = model.matrix(~group,data=metadata_mouse)
fit_mouse_group = getFit(voom = gset, design = des_mouse_group)

tt_mouse_CLZ = topTable(fit_mouse_group, coef = 2, sort.by = "none", confint = T, n=Inf)
tt_mouse_HAL = topTable(fit_mouse_group, coef = 3, sort.by = "none", confint = T, n=Inf)

# Gene mapping between mouse and macaque
mouse_macaque_ens75 = read.delim("mouse_macaque_ens75.tsv")
mouse_affy_moe430a_ens75 = read.delim("mouse_affy_moe430a_ens75.tsv")

gene_data_mouse = fit_mouse_group$genes[,c(1,3,4)]
gene_data_mouse$mouse_ensembl_id = mouse_affy_moe430a_ens75$Ensembl.Gene.ID[match(gene_data_mouse$ID, mouse_affy_moe430a_ens75$Affy.moe430a.probeset)]
gene_data_mouse = data.frame(gene_data_mouse,
                             mouse_macaque_ens75[match(gene_data_mouse$mouse_ensembl_id, mouse_macaque_ens75$Ensembl.Gene.ID),])

gene_data_mouse_f = filter(gene_data_mouse, Macaque.Ensembl.Gene.ID %in% gene_data_mmul_f$Ensembl.Gene.ID, 
                           Homology.Type=="ortholog_one2one")

idx_mouse_genes_f = match(gene_data_mouse_f$ID, fit_mouse_group$genes$ID)

gene_data_mouse_f$Macaque.Ensembl.Gene.ID=factor(gene_data_mouse_f$Macaque.Ensembl.Gene.ID)

# Average multiple ptobes for the same gene 
av_mouse_f = avereps(ex[idx_mouse_genes_f,], ID = gene_data_mouse_f$Macaque.Ensembl.Gene.ID)

fit_av_mouse_f_group = getFit(voom = av_mouse_f, design = des_mouse_group)
tt_av_mouse_f_CLZ = topTable(fit_av_mouse_f_group, coef = 2, sort.by = "none", confint = T, n=Inf)
tt_av_mouse_f_HAL = topTable(fit_av_mouse_f_group, coef = 3, sort.by = "none", confint = T, n=Inf)

# Generate table with DGE results combined for mouse and macaque 
idx_mmul2mouse_f = match(rownames(tt_av_mouse_f_CLZ), fit_mmul_resid1_group$genes$Ensembl.Gene.ID)

df_mouse_mmul_f = data.frame(fit_mmul_resid1_group$genes[idx_mmul2mouse_f,1:2],
                             mouse_CLZ = tt_av_mouse_f_CLZ[,c("logFC","t","P.Value")], 
                             mouse_HAL = tt_av_mouse_f_HAL[,c("logFC","t","P.Value")],
                             mmul_CLZ = tt_mmul_grp_CLZ[idx_mmul2mouse_f,c("logFC","t","P.Value")],
                             mmul_HAL.lo = tt_mmul_grp_HAL.lo[idx_mmul2mouse_f,c("logFC","t","P.Value")],
                             mmul_HAL.hi = tt_mmul_grp_HAL.hi[idx_mmul2mouse_f,c("logFC","t","P.Value")])

# Test significance of correlation
summary(lm(df_mouse_mmul_f$mouse_HAL.t~df_mouse_mmul_f$mmul_HAL.hi.t))
summary(lm(df_mouse_mmul_f$mouse_HAL.t~df_mouse_mmul_f$mmul_HAL.lo.t))
summary(lm(df_mouse_mmul_f$mouse_CLZ.t~df_mouse_mmul_f$mmul_CLZ.t))

# Plot mouse and macaque comparison (Supplementary Fig. S4)
p_mouse_HAL.hi = ggplot(df_mouse_mmul_f, aes(x=mouse_HAL.t, y=mmul_HAL.hi.t)) + geom_point(alpha = 0.05) + geom_smooth(method = "lm") + theme_classic() +
  stat_cor() + xlab("Mouse HAL t-statistic") + ylab("Macaque HAL.hi t-statistic") + ggtitle("Macaque HAL.hi ~ Mouse HAL") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1), plot.title = element_text(hjust = 0.5, face = "bold", size = 15))

p_mouse_HAL.lo = ggplot(df_mouse_mmul_f, aes(x=mouse_HAL.t, y=mmul_HAL.lo.t)) + geom_point(alpha = 0.05) + geom_smooth(method = "lm") + theme_classic() +
  stat_cor() + xlab("Mouse HAL t-statistic") + ylab("Macaque HAL.lo t-statistic") + ggtitle("Macaque HAL.lo ~ Mouse HAL") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1), plot.title = element_text(hjust = 0.5, face = "bold", size = 15))

p_mouse_CLZ = ggplot(df_mouse_mmul_f, aes(x=mouse_CLZ.t, y=mmul_CLZ.t)) + geom_point(alpha = 0.05) + geom_smooth(method = "lm") + theme_classic() +
  stat_cor() + xlab("Mouse CLZ t-statistic") + ylab("Macaque CLZ t-statistic") + ggtitle("Macaque CLZ ~ Mouse CLZ") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1), plot.title = element_text(hjust = 0.5, face = "bold", size = 15))

pdf("sup_fig_mouse_macaque_cor.pdf", w=12,h=4, useDingbats = F)
plot_grid(p_mouse_CLZ, p_mouse_HAL.lo, p_mouse_HAL.hi, align = "h", labels = c(letters[1:3]), label_size = 18, ncol=3)
dev.off()

save.image("apd_macaque_wgcna.rdata")
