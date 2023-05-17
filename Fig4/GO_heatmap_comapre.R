#hepatocyte
dir.create('~/CH/CH-cross/Fig4/hepatocyte')
dir.create('~/CH/CH-cross/Fig4/hepatocyte/fatbody_oenocyte/')
load('~/CH/CH-cross/Fig3/SAMap/SAMap/plot/barplot/species_gene_list.rda')

hu_gene.use=hu_gene.ls[['Hepatocyte']]
ma_gene.use=ma_gene.ls[['Hepatocyte']]
mo_gene.use=mo_gene.ls[['Hepatocyte']]
ze_gene.use=ze_gene.ls[['Hepatocyte']]

dm_gene.use1=dm_gene.ls[['Fatbody']]
dm_gene.use2=dm_gene.ls[['Oenocyte']]

#dm oneocyte and fatbody compare
#dm GO
# gene.use=setdiff(dm_gene.use1,dm_gene.use2)
# gene.use=setdiff(dm_gene.use2,dm_gene.use1)
gene.use=dm_gene.use1
gene_list <- bitr(gene.use,fromType = "SYMBOL",toType =  "ENTREZID",
                  OrgDb =org.Dm.eg.db)
go_list_bp<- enrichGO(gene = gene_list$ENTREZID,
                      OrgDb = org.Dm.eg.db,
                      ont = 'BP',
                      pvalueCutoff =0.01, qvalueCutoff = 0.01,
                      readable = TRUE)
#go_list_bp<- enrichKEGG(gene = gene_list$ENTREZID,organism = 'dme',keyType = 'kegg')
dm_go_1=go_list_bp@result

gene.use=dm_gene.use2
gene_list <- bitr(gene.use,fromType = "SYMBOL",toType =  "ENTREZID",
                  OrgDb =org.Dm.eg.db)
go_list_bp<- enrichGO(gene = gene_list$ENTREZID,
                      OrgDb = org.Dm.eg.db,
                      ont = 'BP',
                      pvalueCutoff =0.01, qvalueCutoff = 0.01,
                      readable = TRUE)
#go_list_bp<- enrichKEGG(gene = gene_list$ENTREZID,organism = 'dme',keyType = 'kegg')
dm_go_2=go_list_bp@result

dm_go_1=dm_go_1[dm_go_1$pvalue<=0.05,]
dm_go_2=dm_go_2[dm_go_2$pvalue<=0.05,]

dm_go_ov=intersect(dm_go_1$Description,dm_go_2$Description)
dm_go_ov=as.data.frame(dm_go_ov)

go_use=c('immune system process','carboxylic acid catabolic process','alpha-amino acid metabolic process','cellular amino acid biosynthetic process','Toll signaling pathway','fatty acid metabolic process','very long-chain fatty acid biosynthetic process','very long-chain fatty acid metabolic process','biosynthetic process of antibacterial peptides active against Gram-negative bacteria')

dm_go_use_1=dm_go_1[dm_go_1$Description%in%go_use,]
dm_go_use_2=dm_go_2[dm_go_2$Description%in%go_use,]
dm_go_use_1$score=-log10(dm_go_use_1$pvalue)
dm_go_use_2$score=-log10(dm_go_use_2$pvalue)

dm_go_use_1$ct='Fatbody'
dm_go_use_2$ct='Oenocyte'

dm_go_use=rbind(dm_go_use_1,dm_go_use_2)
dm_go_use=dm_go_use[,c('Description','score','ct')]

colnames(dm_go_use)[2]='value'
library(tidyr)
library(tibble)
mat=pivot_wider(dm_go_use, names_from = Description) %>%column_to_rownames("ct") %>%as.matrix()
mat=t(mat)
mat[is.na(mat)]=0

pal=colorRampPalette(c('#3f6a94','white','#ad486f'))(1000)

mycolor = c(brewer.pal(8,'Set3'))
mycolor = colorRampPalette(mycolor)(10)
show_col(mycolor)
mycolor=mycolor[c(1,2)]

ct=unique(colnames(mat))
color_df=data.frame('ct'= ct)
rownames(color_df)=ct

names(mycolor)=ct
ann_color=list(ct=mycolor)
mat=mat[go_use,]
p=pheatmap::pheatmap(mat,
                     color = paletteContinuous("whiteBlue"),
                     scale = 'column',
                     cluster_rows = F,
                     cluster_cols = F,
                     angle_col = 45,
                     fontsize = 10,
                     width = 3,
                     height = 3,show_rownames = T,annotation_colors = ann_color,annotation_col =color_df,border_color = 'NA')
p
