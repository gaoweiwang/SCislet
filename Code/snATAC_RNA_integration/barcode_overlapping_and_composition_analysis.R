blank_theme <- theme_minimal()+
  theme(
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold"),
    panel.background = element_rect(fill = "white", colour = "grey50")
  )
################
##############barcode overlapping and composition analysis (Figure 1c, Figure S1g-h, Figure S6j, k, l)
################using Figure 1 as example, similar procedure is applied to Figure 6 using this matrix /oasis/tscc/scratch/hazhu/share/upload/primary.sc.islet.atac.rna.labels.csv
integrated.mtx<-read.csv("/oasis/tscc/scratch/hazhu/share/upload/sc.islet.diff.atac.rna.labels.csv",row.names = 1)
##############barcode overlapping (Figure 1c)
integrated.mtx$snATAC.seq.ident<-factor(integrated.mtx$snATAC.seq.ident,levels = c("PP-1","PP-2",
                                                                                   "ENP-1","ENP-alpha","ENP-2","ENP-3",
                                                                                   "SC-alpha","SC-EC","SC-beta","SC-delta"))
integrated.mtx$predicted.scRNA.seq.ident<-factor(integrated.mtx$predicted.scRNA.seq.ident,levels = c("PP-1","PP-2",
                                                                                   "ENP-1","ENP-alpha","ENP-2","ENP-3",
                                                                                   "SC-alpha","SC-EC","SC-beta","SC-delta"))

id.match.ratio<-prop.table(table(integrated.mtx$snATAC.seq.ident,integrated.mtx$predicted.scRNA.seq.ident),margin = 1)
heatmap.2(id.match.ratio,
          Rowv = NULL, Colv = NULL,
          trace="none",
          col=colorRampPalette(brewer.pal(8, "Blues"))(25),
          colsep=c(1:9),
          rowsep=c(1:9),
          sepcolor="white",
          sepwidth=c(0.05,0.05))
##############composition analysis Figure S1g
col_dna_1<-c("ENP-1"="thistle1","ENP-3"="darkseagreen3","PP-1"="wheat1","PP-2"="olivedrab3",
             "ENP-alpha"="khaki3","ENP-2"="darkgreen","SC-alpha"="purple2",
             "SC-EC"="tan1","SC-beta"="royalblue2","SC-delta"="red2")

atac.composition<-prop.table(table(integrated.mtx$stage,integrated.mtx$snATAC.seq.ident),margin = 1)
l.atac.composition<-melt(atac.composition)
colnames(l.atac.composition)<-c("stage","cell.type","composition")
l.atac.composition<-as.data.frame(l.atac.composition)
l.atac.composition$cell.type<-factor(l.atac.composition$cell.type,levels = c("PP-1","PP-2","ENP-1",
                                                                             "ENP-alpha","ENP-2","ENP-3",
                                                                             "SC-alpha","SC-EC","SC-beta","SC-delta"))
bp<- ggplot(l.atac.composition, aes(x=stage, y=composition, fill=cell.type))+
  geom_bar(width = 0.8, stat = "identity")+scale_fill_manual(values =col_dna_1)
bp+blank_theme




