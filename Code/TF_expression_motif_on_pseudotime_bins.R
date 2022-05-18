library(Seurat)
library(Signac)
library(monocle3)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(reshape2)
library(cicero)
library(gplots)
library(uwot)
library(dplyr)
library(plyr)
blank_theme <- theme_minimal()+
  theme(
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold"),
    panel.background = element_rect(fill = "white", colour = "grey50")
  )

#################Gene expression/TF motif enrichment on pseudotime bins Figure 1i,j
#################Example shown for Branch point 1 in 1Figure 1i,analysis of Branch point 2 in Figure 1j is similar, but using ENP3, SC-EC and SC-beta-cells in D14, D21 
###########ChromVAR based motif enrichment score matrix using snATAC-seq data UMAP embeddings
diff_ds<-readRDS("/oasis/tscc/scratch/hazhu/share/upload/diff_atac_chromvar.rds")
###########Inferred gene expression matrix using snATAC-seq data UMAP embeddings
diff_rna<-readRDS("/oasis/tscc/scratch/hazhu/share/upload/diff_rna.rds")

DefaultAssay(diff_ds)<-"TF_motif"
FeaturePlot(diff_ds,features = "MA0069.1-PAX6",cols = c("grey85", "thistle1","Red4"))
Idents(diff_ds)<-"anno"
diff_ds_sub<-subset(diff_ds,subset=anno%in%c("diff_ENP1","diff_ENP3",
                                             "diff_ENPalpha","diff_ENP2","diff_alpha",
                                             "diff_EC","diff_beta","diff_delta"))
rm(diff_ds)
######take cells from D11 and D14
diff_ds_sub_2<-subset(diff_ds_sub,subset=time%in%c("D11","D14"))
dna_all_cds <- as.cell_data_set(diff_ds_sub_2)
dna_all_cds <- cluster_cells(cds = dna_all_cds, reduction_method = "UMAP")
dna_all_cds <- learn_graph(dna_all_cds, use_partition = F,close_loop = F,
                           learn_graph_control=list(ncenter=500,minimal_branch_len=10))

dna_all_cds <- order_cells(dna_all_cds, reduction_method = "UMAP")
######Extract gene expression information from the same set of cells after data integration
dna_all_cds.rna <- as.cell_data_set(diff_rna)
dna_all_cds.rna@assays@data$counts<-dna_all_cds.rna@assays@data$logcounts
rm(diff_rna)

################node1
###ENP-alpha
cds_subset<-dna_all_cds[,dna_all_cds@colData$anno%in%c("diff_ENP1","diff_ENPalpha","diff_alpha")]
######manually select alpha cell lineage
cds_subset <- choose_graph_segments(cds_subset)
cds_subset<-dna_all_cds[,colnames(cds_subset)]
pData(cds_subset)$Pseudotime <- pseudotime(cds_subset)
pData(cds_subset)$cell_subtype <- cut(pseudotime(cds_subset), 12)
bin_cds_subset <- aggregate_by_cell_bin(cds_subset, "cell_subtype")

d11.d14.enpalpha.chromvar<-bin_cds_subset@assays@data$counts

######extract RNA expression on the same lineage trajectory
dna_all_cds.rna.1<-dna_all_cds.rna[,colnames(cds_subset)]
pData(dna_all_cds.rna.1)$Pseudotime <- pData(cds_subset)$Pseudotime
pData(dna_all_cds.rna.1)$cell_subtype <- pData(cds_subset)$cell_subtype
bin_cds_subset.rna <- aggregate_by_cell_bin(dna_all_cds.rna.1, "cell_subtype")
d11.d14.enpalpha.rna<-bin_cds_subset.rna@assays@data$counts
rownames(d11.d14.enpalpha.rna)<-rownames(bin_cds_subset.rna)
#######trajectory cell composition
d11.d14.enpalpha.celltype<-cds_subset@colData[,c("anno","cell_subtype")]
prop.d11.d14.enpalpha.celltype<-tapply(d11.d14.enpalpha.celltype[,1], d11.d14.enpalpha.celltype[,2], function(x) prop.table(table(x)))
prop.d11.d14.enpalpha.celltype<-read.csv("OneDrive - UC San Diego/beta_cell_maturation/Single_cell/snATAC/pseudotime/new_time/prop.d11.d14.enpalpha.celltype.csv",row.names = 1)
prop.d11.d14.enpalpha.celltype<-cbind("pseudotime"=paste("pseudotime",c(1:12),sep = ""),prop.d11.d14.enpalpha.celltype)
prop.d11.d14.enpalpha.celltype.l<-melt(prop.d11.d14.enpalpha.celltype)
colnames(prop.d11.d14.enpalpha.celltype.l)<-c("pseudotime","anno","comp")
prop.d11.d14.enpalpha.celltype.l$pseudotime<-factor(prop.d11.d14.enpalpha.celltype.l$pseudotime,levels = paste("pseudotime",c(1:12),sep = ""))
col_dna_1<-c("ENP1"="thistle1",
             "ENPalpha"="mediumorchid2","Scalpha"="purple2")

bp<- ggplot(prop.d11.d14.enpalpha.celltype.l, aes(x=pseudotime, y=comp, fill=anno))+
  geom_bar(width = 1, stat = "identity")+scale_fill_manual(values =col_dna_1)
bp+blank_theme

###ENP-2
cds_subset<-dna_all_cds[,dna_all_cds@colData$anno%in%c("diff_ENP1","diff_ENP2")]
#######manually select non-alphal cell lineage
cds_subset <- choose_graph_segments(cds_subset)
cds_subset<-dna_all_cds[,colnames(cds_subset)]
pData(cds_subset)$Pseudotime <- pseudotime(cds_subset)
pData(cds_subset)$cell_subtype <- cut(pseudotime(cds_subset), 12)
bin_cds_subset <- aggregate_by_cell_bin(cds_subset, "cell_subtype")
d11.d14.enp2.chromvar<-bin_cds_subset@assays@data$counts

rownames(d11.d14.enp2.chromvar)<-rownames(bin_cds_subset)

######extract RNA expression on the same lineage trajectory
dna_all_cds.rna.2<-dna_all_cds.rna[,colnames(cds_subset)]
pData(dna_all_cds.rna.2)$Pseudotime <- pData(cds_subset)$Pseudotime
pData(dna_all_cds.rna.2)$cell_subtype <- pData(cds_subset)$cell_subtype
bin_cds_subset.rna <- aggregate_by_cell_bin(dna_all_cds.rna.2, "cell_subtype")
d11.d14.enp2.rna<-bin_cds_subset.rna@assays@data$counts
rownames(d11.d14.enp2.rna)<-rownames(bin_cds_subset.rna)

#####ENP2 cell type composition
d11.d14.enp2.celltype<-cds_subset@colData[,c("anno","cell_subtype")]
prop.d11.d14.enp2.celltype<-tapply(d11.d14.enp2.celltype[,1], d11.d14.enp2.celltype[,2], function(x) prop.table(table(x)))

prop.d11.d14.enp2.celltype<-read.csv("OneDrive - UC San Diego/beta_cell_maturation/Single_cell/snATAC/pseudotime/new_time/prop.d11.d14.enp2.celltype.csv",row.names = 1)
prop.d11.d14.enp2.celltype<-cbind("pseudotime"=paste("pseudotime",c(1:12),sep = ""),prop.d11.d14.enp2.celltype)
prop.d11.d14.enp2.celltype.l<-melt(prop.d11.d14.enp2.celltype)
colnames(prop.d11.d14.enp2.celltype.l)<-c("pseudotime","anno","comp")
prop.d11.d14.enp2.celltype.l$pseudotime<-factor(prop.d11.d14.enp2.celltype.l$pseudotime,levels = paste("pseudotime",c(1:12),sep = ""))
col_dna_2<-c("ENP1"="thistle1","ENP2"="tan1")

bp<- ggplot(prop.d11.d14.enp2.celltype.l, aes(x=pseudotime, y=comp, fill=anno))+
  geom_bar(width = 1, stat = "identity")+scale_fill_manual(values =col_dna_2)
bp+blank_theme

######combine two lineages

d11.d14.chromvar<-cbind(d11.d14.enp2.chromvar[,order(12:1)],d11.d14.enpalpha.chromvar)
colnames(d11.d14.chromvar)<-paste("pseudotime",c(1:24),sep = "")

d11.d14.chromvar.scale<-apply(d11.d14.chromvar,1,scale)
d11.d14.chromvar.scale<-t(d11.d14.chromvar.scale)
########plot key motifs
d11.d14.chromvar.scale.ordered.key<-d11.d14.chromvar.scale.ordered[c("MA0077.1-SOX9",
                                                                     "MA1515.1-KLF2",
                                                                     "MA1121.1-TEAD2",
                                                                     "MA0774.1-MEIS2",
                                                                     "MA0669.1-NEUROG2",
                                                                     "MA0069.1-PAX6",
                                                                     "MA0761.2-ETV1",
                                                                     "MA0798.2-RFX3",
                                                                     "MA1593.1-ZNF317",
                                                                     "MA0682.2-PITX1",
                                                                     "MA0132.2-PDX1",
                                                                     "MA0674.1-NKX6.1",
                                                                     "MA0068.2-PAX4",
                                                                     "MA0484.2-HNF4G",
                                                                     "MA0702.2-LMX1A",
                                                                     "MA0466.2-CEBPB",
                                                                     "MA0693.2-VDR"),]

hmp.var<-heatmap.2(as.matrix(d11.d14.chromvar.scale.ordered.key),
                   trace="none",
                   density="none",
                   Colv=NA,Rowv=NA,
                   key = F,
                   col=color.palette.1,
                   breaks=100,
                   lmat=rbind(c(0, 3, 4), c(2,1,0 )), 
                   lwid=c(0.5, 1.5, 1.2),
                   lhei=c(0,6))


###########corresponding TF RNA expression
d11.d14.rna<-cbind(d11.d14.enp2.rna[,order(12:1)],d11.d14.enpalpha.rna)
colnames(d11.d14.rna)<-paste("pseudotime",c(1:24),sep = "")

d11.d14.rna.key<-d11.d14.rna[c("SOX9",
                               "MEIS2","KLF2","TEAD2",
                               "MEIS2","NEUROG3",
                               "PAX6","ETV1",
                               "RFX3","RFX6","ZNF317",
                               "PITX1",
                               "PDX1","NKX6-1","PAX4",
                               "HNF4G","LMX1A",
                               "CEBPB","VDR"),]
d11.d14.rna.scale<-apply(d11.d14.rna.key,1,scale)
d11.d14.rna.scale<-t(d11.d14.rna.scale)
color.palette.2 = colorRampPalette(c("grey30","grey90","YellowGreen","Forestgreen"), space="Lab")

hmp.var<-heatmap.2(as.matrix(d11.d14.rna.scale),
                   trace="none",
                   density="none",
                   Colv=NA,Rowv=NA,
                   key = F,
                   col=color.palette.2,
                   breaks=100,
                   lmat=rbind(c(0, 3, 4), c(2,1,0 )), 
                   lwid=c(0.5, 1.5, 1.2),
                   lhei=c(0,6))

