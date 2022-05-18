library(Seurat)
library(Signac)
library(monocle3)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(reshape2)
library(cicero)
library(gplots)
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

################ temporal order of transcriptional programs, Figure 3
###use alpha cell lineage as an example
load("/oasis/tscc/scratch/hazhu/share/upload/diff_atac_chromvar.rds")
diff_ds_sub<-subset(diff_ds,subset=anno%in%c("diff_ENP1","diff_ENPalpha","diff_alpha"))
col_dna<-c("diff_ENP1"="thistle1",
           "diff_ENPalpha"="mediumorchid2","diff_alpha"="purple2")

dna_all_cds <- as.cell_data_set(diff_ds_sub)
dna_all_cds <- cluster_cells(cds = dna_all_cds, reduction_method = "UMAP")
dna_all_cds <- learn_graph(dna_all_cds, use_partition = F,close_loop = F,
                           learn_graph_control=list(ncenter=500,minimal_branch_len=10))

dna_all_cds <- order_cells(dna_all_cds, reduction_method = "UMAP")
######manually select alpha cell lineage
cds_subset<-choose_graph_segments(dna_all_cds)
cds_subset<-dna_all_cds[,colnames(cds_subset)]
p2<-plot_cells(dna_all_cds,
               color_cells_by = "anno",
               label_groups_by_cluster=FALSE,
               trajectory_graph_color = "grey0",
               trajectory_graph_segment_size = 1.5,
               group_label_size = 5,
               label_cell_groups = FALSE,
               label_leaves=FALSE,
               label_branch_points=FALSE,
               label_roots = FALSE)+ scale_color_manual(values =col_dna )+NoLegend()



alpha.pseudotime.mtx<-cds_subset@assays@data$counts

alpha.pseudotime.mtx<-alpha.pseudotime.mtx[rownames(alpha.pseudotime.mtx)%in%rownames(ccre.umap)[ccre.umap$module%in%c("ENP1","ENP_alpha","SC_alpha")],]
alpha.pseudotime.mtx<-as.data.frame(t(alpha.pseudotime.mtx))
alpha.pseudotime.mtx<-cbind(alpha.pseudotime.mtx,"pseudotime"=pseudotime(cds_subset))
alpha.pseudotime.mtx<-alpha.pseudotime.mtx[order(alpha.pseudotime.mtx$pseudotime),]

#####cCRE pseudotime
alpha.ccre.pseudotime<-matrix(0,0,2)
for (i in colnames(alpha.pseudotime.mtx[,-5961])) {
  temp.pseudotime<-alpha.pseudotime.mtx$pseudotime[alpha.pseudotime.mtx[,i]==1]
  pseodotime.sumit<-density(temp.pseudotime)$x[which.max(density(temp.pseudotime)$y)]
  temp<-c(i,pseodotime.sumit)
  alpha.ccre.pseudotime<-rbind(alpha.ccre.pseudotime,temp)
}
rownames(alpha.ccre.pseudotime)<-colnames(alpha.pseudotime.mtx[,-5961])
alpha.pseudotime.umap<-ccre.umap[rownames(alpha.ccre.pseudotime),]
alpha.pseudotime.umap$pseudotime<-as.numeric(alpha.ccre.pseudotime[,2])

############Plot cCRE pseudotime on cCRE UMAP, fupplementary figure 3j
ccre.umap<-read.csv("/oasis/tscc/scratch/hazhu/share/upload/sc.islet.diff.ccre.umap.csv",row.names = 1)
p.umap.pseudotime<-ggplot(ccre.umap,aes(x=UMAP1,y=UMAP2))+geom_point(size=0.01,colour = "grey80")+
  geom_point(data = alpha.pseudotime.umap,aes(x=UMAP1,y=UMAP2,color=pseudotime),size=0.01)+
  scale_colour_viridis_c(option = "inferno")+
  blank_theme

############construct RNA pseudotime
diff_rna<-readRDS("/oasis/tscc/scratch/hazhu/share/upload/diff_rna.rds")
rna_all_cds <- as.cell_data_set(diff_rna)
rna_all_cds <- cluster_cells(cds = rna_all_cds, reduction_method = "UMAP")
rna_subset<-rna_all_cds[,colnames(cds_subset)]
rm(rna_all_cds)
rm(diff_rna)

#######look for alpha genes and alpha TFs
alpha.peaks<-colnames(alpha.pseudotime.mtx[,1:(ncol(alpha.pseudotime.mtx)-1)])
alpha.genes<-sig_target_gene_predicton[sig_target_gene_predicton$peak_id%in%alpha.peaks,"gene"]
############include master TFs found in Figure 2d
master.tfs<-read.csv("/oasis/tscc/scratch/hazhu/share/upload/sc.islet.diff.module.specific.tf.sig.csv",row.names = 1)
alpha.master.TFs<-master.tfs[master.tfs$module%in%c("ENP1","ENP_alpha","SC_alpha"),"TF"]
alpha.genes<-unique(c(alpha.genes,alpha.master.TFs))

alpha.rna.pseudotime.mtx<-rna_subset@assays@data$logcounts
alpha.rna.pseudotime.mtx<-alpha.rna.pseudotime.mtx[rownames(alpha.rna.pseudotime.mtx)%in%alpha.genes,]
alpha.rna.pseudotime.mtx<-as.data.frame(t(alpha.rna.pseudotime.mtx))
alpha.rna.pseudotime.mtx<-cbind(alpha.rna.pseudotime.mtx,"pseudotime"=pseudotime(cds_subset))
alpha.rna.pseudotime.mtx<-alpha.rna.pseudotime.mtx[order(alpha.rna.pseudotime.mtx$pseudotime),]

alpha.rna.pseudotime<-matrix(0,0,2)
for (i in colnames(alpha.rna.pseudotime.mtx[,1:(ncol(alpha.rna.pseudotime.mtx)-1)])) {
  temp.fit<-smooth.spline(alpha.rna.pseudotime.mtx$pseudotime,alpha.rna.pseudotime.mtx[,colnames(alpha.rna.pseudotime.mtx)==i])
  temp.pseudotime<-temp.fit$x[which.max(temp.fit$y)]
  temp<-c(i,temp.pseudotime)
  alpha.rna.pseudotime<-rbind(alpha.rna.pseudotime,temp)
}
rownames(alpha.rna.pseudotime)<-colnames(alpha.rna.pseudotime.mtx[,1:(ncol(alpha.rna.pseudotime.mtx)-1)])
colnames(alpha.rna.pseudotime)<-c("names","pseudotime")
alpha.rna.pseudotime<-as.data.frame(alpha.rna.pseudotime)
alpha.rna.pseudotime$pseudotime<-as.numeric(alpha.rna.pseudotime$pseudotime)
alpha.rna.pseudotime$module<-NA
alpha.rna.pseudotime$origin<-rep("RNA",nrow(alpha.rna.pseudotime))
colnames(alpha.ccre.pseudotime)<-c("names","pseudotime")
alpha.ccre.pseudotime<-as.data.frame(alpha.ccre.pseudotime)
alpha.ccre.pseudotime$pseudotime<-as.numeric(alpha.ccre.pseudotime$pseudotime)
alpha.ccre.pseudotime$module<-alpha.pseudotime.umap$module
alpha.ccre.pseudotime$origin<-rep("cCRE",nrow(alpha.ccre.pseudotime))


alpha.all.pseudotime<-rbind(alpha.rna.pseudotime,alpha.ccre.pseudotime)
alpha.all.pseudotime$axis<-rep(1,nrow(alpha.all.pseudotime))
alpha.all.pseudotime<-alpha.all.pseudotime[order(alpha.all.pseudotime$pseudotime),]
alpha.all.pseudotime$rank<-c(1:nrow(alpha.all.pseudotime))

#######plot TF, TF target cCREs and TF target genes
TF.gene<-"ZNF414"
TF.bs<-C_TF[C_TF$gene==TF.gene & C_TF$stat>0.5,]
TF.bs<-TF.bs[TF.bs$peak_ID%in%rownames(alpha.all.pseudotime),]
alpha.all.pseudotime$scc<-rep(0,nrow(alpha.all.pseudotime))
for (i in TF.bs$peak_ID) {
  alpha.all.pseudotime$scc[rownames(alpha.all.pseudotime)==i]<-TF.bs[TF.bs$peak_ID==i,"stat"]
}
pseudotime.TF.bs<-alpha.all.pseudotime[!alpha.all.pseudotime$scc==0,]
tf.exp<-alpha.all.pseudotime[TF.gene,]
ccre.target<-sig_target_gene_predicton[sig_target_gene_predicton$peak_id%in%rownames(pseudotime.TF.bs),"gene"]
ccre.target<-unique(ccre.target)
ccre.target.pseudotime<-alpha.all.pseudotime[ccre.target,]
p.pseudotime.2<-ggplot(alpha.all.pseudotime,aes(x=pseudotime,y=origin))+geom_point(size=0.5,colour="grey85")+
  geom_point(data = ccre.target.pseudotime,aes(x=pseudotime,y=origin),size=0.5,colour="tan3")+
  geom_point(data =tf.exp, aes(x=pseudotime,y=origin),size=2, colour="forestgreen")+
  geom_point(data=pseudotime.TF.bs,aes(x=pseudotime,y=origin,color=scc),size=0.5)+
  scale_color_gradientn(colours = c("grey90","Tomato1","red4"))+
  blank_theme



