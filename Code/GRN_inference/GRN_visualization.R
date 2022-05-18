library(ggplot2)
library(patchwork)
library(reshape2)
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

#################GRN on cCRE UMAP.
#####Using sc.islet.diff GRN as example, similar analysis can be applied to the comparison between SC and primary beta cells, using similar input files start with sc.primary.islet**
###cCRE UMAP, Figure 2b
ccre.umap<-read.csv("/oasis/tscc/scratch/hazhu/share/upload/sc.islet.diff.ccre.umap.csv",row.names = 1)
col_cCRE<-c("ENP1"="thistle1","ENP_alpha"="khaki3","ENP3"="darkseagreen3",
            "SC_alpha"="purple2","SC_EC"="tan1","SC_beta"="royalblue2","SC_delta"="red2")
p.umap<-ggplot(ccre.umap,aes(x=UMAP1,y=UMAP2,color=module))+geom_point(size=0.01)+
  blank_theme+scale_color_manual(values =col_cCRE)
p.umap
####heatmap, Figure 2c use ENP1/2 module as example
sig_target_gene_predicton<-read.csv("/oasis/tscc/scratch/hazhu/share/upload/sc.islet.diff.cCRE.target.gene.pairs.csv",row.names = 1)
gene_cpm_mtx_ordered<-read.csv("/oasis/tscc/scratch/hazhu/share/upload/sc.islet.diff.gene.cpm.mtx.csv",row.names = 1)
#####plot cCRE accessibility
enp1.peaks<-as.matrix(ccre.umap[ccre.umap$module=="ENP1",c(3:14)])
hclu.enp1.peaks<-hclust(as.dist(1-abs(cor(t(enp1.peaks),method='spearman'))),method='ward')
enp1.peaks<-enp1.peaks[hclu.enp1.peaks$order,]
color.palette.1 = colorRampPalette(c("grey30","grey90", "CornflowerBlue","Royalblue4"), space="Lab")
hmp.var<-heatmap.2(enp1.peaks,trace="none",density="none",
                   Colv=NA,Rowv=NA,col=color.palette.1,breaks=100)

#####plot target gene expression
enp1.genes<-matrix(0,0,12)
for (i in 1:nrow(enp1.peaks)) {
  peak.id<-rownames(enp1.peaks)[i]
  gene.syb<-sig_target_gene_predicton[sig_target_gene_predicton$peak_id==peak.id,"gene"]
  temp<-as.data.frame(gene_cpm_mtx_ordered[gene.syb,4:15])
  enp1.genes<-rbind(enp1.genes,temp)
}
enp1.genes<-as.matrix(enp1.genes)
enp1.genes.scale<-apply(enp1.genes,1,scale)
enp1.genes.scale<-t(enp1.genes.scale)
colnames(enp1.genes.scale)<-colnames(enp1.genes)

color.palette.2 = colorRampPalette(c("grey30","grey90","YellowGreen","Forestgreen"), space="Lab")
hmp.var<-heatmap.2(enp1.genes.scale,trace="none",density="none",
                   Colv=NA,Rowv=NA,col=color.palette.2,breaks=100)

#######Identify key TF regulators of each lineage
C_TF.pos<-read.csv("/oasis/tscc/scratch/hazhu/share/upload/sc.islet.diff.TF.cCRE.pairs.csv",row.names = 1)
peak.TF.cor<-matrix(0,nrow(ccre.umap),0)
for (i in unique(C_TF.pos$gene)) {
  temp.C.TF<-C_TF.pos[C_TF.pos$gene==i,]
  peak.ids<-rownames(ccre.umap)[rownames(ccre.umap)%in%temp.C.TF$peak_ID]
  temp.C.TF<-temp.C.TF[match(peak.ids,temp.C.TF$peak_ID),]
  temp.cor<-rep(0,nrow(ccre.umap))
  temp.cor[rownames(ccre.umap)%in%temp.C.TF$peak_ID]<-temp.C.TF$stat
  peak.TF.cor<-cbind(peak.TF.cor,temp.cor)
}
rownames(peak.TF.cor)<-rownames(ccre.umap)
colnames(peak.TF.cor)<-unique(C_TF.pos$gene)
t.peak.TF.cor<-t(peak.TF.cor)
peak.TF.cor<-as.data.frame(peak.TF.cor)

#####use fisher exact test
peak.TF<-!peak.TF.cor==0
peak.TF<-as.data.frame(peak.TF)
peak.TF<-peak.TF[,!apply(peak.TF, 2,function(x) sum(as.logical(x)))==0]

########how to extract master TFs for each peak module
C_TF_crn<-C_TF[C_TF$peak_ID%in%rownames(ccre.umap),]
C_TF_pos<-C_TF_crn[C_TF_crn$stat>0.3,]
tfs<-unique(C_TF_pos$gene)

####fisher exact test
peak.TF$module<-ccre.umap$module
peak.TF$module<-factor(peak.TF$module,levels = unique(peak.TF$module))
module.specific.tf<-matrix(0,0,4)
for(i in unique(peak.TF$module)) for (j in 1:266) {
  test<-table(peak.TF[,j],peak.TF$module)
  test<-test[,colnames(test)==i]
  backgound<-table(peak.TF[,j])
  stat.tfbs<-fisher.test(cbind(backgound,test))
  temp<-c(i,colnames(peak.TF)[j],stat.tfbs$estimate,stat.tfbs$p.value)
  module.specific.tf<-rbind(module.specific.tf,temp)
}
module.specific.tf<-as.data.frame(module.specific.tf)
colnames(module.specific.tf)<-c("module","TF","odd.ratio","p.value")
rownames(module.specific.tf)<-paste(module.specific.tf$module,module.specific.tf$TF,sep = "_")
module.specific.tf$module<-factor(module.specific.tf$module,levels = unique(peak.TF$module))
module.specific.tf$FDR<-p.adjust(module.specific.tf$p.value,method = "BH")
module.specific.tf$logFDR<-(-log10(module.specific.tf$FDR))
module.specific.tf.sig<-module.specific.tf[module.specific.tf$odd.ratio>1&module.specific.tf$logFDR>1.3,]

#########plot key TF regulators in Figure 2d
module.specific.tf.key<-module.specific.tf[module.specific.tf$TF%in%c("NEUROG3",
                                                                      "RFX6",
                                                                      "PITX1","PAX6","ETV1","MEF2A",
                                                                      "RFX3","CDX2",
                                                                      "FEV","ASCL2",
                                                                      "NKX6-1","RXRA",
                                                                      "LMX1A","HNF4G",
                                                                      "ISL1","BACH2",
                                                                      "PDX1","SCRT1",
                                                                      "GRHL1"),]
module.specific.tf.key<-module.specific.tf.key[,c(1:3,6)]
module.specific.tf.key$TF<-factor(module.specific.tf.key$TF,levels = c("NEUROG3",
                                                                       "RFX6",
                                                                       "PITX1","PAX6","ETV1","MEF2A",
                                                                       "RFX3",
                                                                       "FEV","CDX2","ASCL2",
                                                                       "NKX6-1","RXRA",
                                                                       "LMX1A","HNF4G",
                                                                       "ISL1","BACH2",
                                                                       "PDX1","SCRT1",
                                                                       "GRHL1"))
module.specific.tf.key$module<-factor(module.specific.tf.key$module,levels = c("SC_delta","SC_beta","SC_EC",
                                                                               "SC_alpha","ENP3","ENP_alpha","ENP1"))
module.specific.tf.key$odd.ratio<-as.numeric(module.specific.tf.key$odd.ratio)
module.specific.tf.key$logFDR<-as.numeric(module.specific.tf.key$logFDR)
p<-ggplot(module.specific.tf.key,aes(x=TF,y=module))+
  geom_point(aes(fill=logFDR,size=odd.ratio),shape=21)+
  scale_size(range = c(1, 10))+
  scale_fill_gradientn(colors=c("white","gold1","gold2","goldenrod1","goldenrod2","goldenrod3","tan3","tan4","salmon4","navy"))+
  blank_theme
p

###Find cCRE module specific TF interaction, Figure 2i-k
TF.dependency.odd<-matrix(0,0,5)
colnames(TF.dependency.odd)<-c("test.TF","background.TF","ccre.module","odd.ratio","p.value")
peak.TF<-peak.TF[,colnames(peak.TF)%in%unique(module.specific.tf.sig$TF)]

for (m in unique(module.specific.tf.sig$module)) {
  peak.TF.module<-peak.TF[,colnames(peak.TF)%in%unique(module.specific.tf.sig$TF[module.specific.tf.sig$module==m])]
  peak.TF.module<-peak.TF.module[rownames(peak.TF.module)%in%rownames(ccre.umap[ccre.umap$module==m,]),]
  for (i in 1:ncol(peak.TF.module)) for (j in 1:ncol(peak.TF.module)){
    background<-table(peak.TF.module[,i])
    dependent.bs<-peak.TF.module[peak.TF.module[,j]==TRUE,]
    test<-table(dependent.bs[,i])
    if (dim(test)==1) {
      if (names(test)=="FALSE") {
        test[2]<-0
        names(test)<-c("FALSE","TRUE")
      } else {
        test[2]<-test[1]
        test[1]<-0
        names(test)<-c("FALSE","TRUE")
      }
    }
    stat.dependency<-fisher.test(cbind(backgound,test))
    temp<-c(colnames(peak.TF.module)[i],colnames(peak.TF.module)[j],m,stat.dependency$estimate,stat.dependency$p.value)
    TF.dependency.odd<-rbind(TF.dependency.odd,temp)
  }
}
TF.dependency.odd<-as.data.frame(TF.dependency.odd)
TF.dependency.odd<-TF.dependency.odd[!TF.dependency.odd$test.TF==TF.dependency.odd$background.TF,]
TF.dependency.odd$FDR<-p.adjust(TF.dependency.odd$p.value,method = "BH")
TF.dependency.odd$logFDR<-(-log10(TF.dependency.odd$FDR))
sig.TF.dependency.odd<-TF.dependency.odd[TF.dependency.odd$logFDR>5&TF.dependency.odd$odd.ratio>20,]

##plot Figure 2j
tf1<-"RXRA"
tf2<-"NR1H4"
ccre.module<-c("SC_EC","SC_beta")
tf.interaction.umap<-ccre.umap[,1:2]
tf.interaction.umap$tf.intersect<-rep("none",nrow(tf.interaction.umap))
tf.interaction.umap$tf.intersect[peak.TF[,tf1]]<-"TF.1"
tf.interaction.umap$tf.intersect[peak.TF[,tf2]]<-"TF.2"
tf.interaction.umap$tf.intersect[peak.TF[,tf1]&peak.TF[,tf2]]<-"co.bind"
tf.interaction.umap$tf.intersect<-factor(tf.interaction.umap$tf.intersect,levels = c("none","TF.1","TF.2","co.bind"))
col_intersect<-c("TF.1"="green","TF.2"="red","co.bind"="yellow","none"="grey90")
p.umap<-ggplot(tf.interaction.umap,aes(x=UMAP1,y=UMAP2,color=tf.intersect))+geom_point(size=0.1)+
  geom_point(data = ccre.umap[ccre.umap$module%in%ccre.module,],aes(x=UMAP1,y=UMAP2),color="grey65",size=1)+
  geom_point(data = tf.interaction.umap[tf.interaction.umap$tf.intersect=="TF.1",],aes(x=UMAP1,y=UMAP2,color=tf.intersect),size=0.5)+
  geom_point(data = tf.interaction.umap[tf.interaction.umap$tf.intersect=="TF.2",],aes(x=UMAP1,y=UMAP2,color=tf.intersect),size=0.5)+
  geom_point(data = tf.interaction.umap[tf.interaction.umap$tf.intersect=="co.bind",],aes(x=UMAP1,y=UMAP2,color=tf.intersect),size=0.5)+
  blank_theme+scale_color_manual(values =col_intersect)
p.umap


