setwd('~/CH/CH-cross/Fig4/hepatocyte/')

load('~/CH/CH-cross/Fig3/SAMap/SAMap/plot/branch/ct_gene.ls.rds')

names(ew_gene.ls)
names(hu_gene.ls)


ea_gene.use=ew_gene.ls[['Erythrocruorin']]
hu_gene.use=hu_gene.ls[['Hepatocyte']]

library(data.table)
gene_pair=fread('~/CH/CH-cross/Fig3/SAMap/maps/huea/ea_to_hu.txt')
gene_pair=gene_pair[gene_pair$V3>=30,]

hu_ea_orth=gene_pair[gene_pair$V1%in%ea_gene.use,]$V2
hu_ea_orth=hu_ea_orth[!duplicated(hu_ea_orth)]

library(VennDiagram)
gene_list=list()
gene_list[['hu_Hepatocyte']]=hu_gene.use
gene_list[['ea_Erythrocruorin']]=hu_ea_orth
gene_ov=intersect(hu_gene.use,hu_ea_orth)
gene_ov[order(gene_ov)]

library(scales)
fill_colors=c('#0073C2FF','#EFC000FF')
show_col(fill_colors)
p=venn.diagram(gene_list,col='white',filename = NULL,fill=fill_colors,lwd=0.5,width = 2000,height = 2000,cex=0.5,cat.cex=0.5)
p
dir.create('~/CH/CH-cross/Fig4/hepatocyte/ea_hu/')
pdf('~/CH/CH-cross/Fig4/hepatocyte/ea_hu/venn_ew_hu.pdf',height = 10,width = 10)
grid.newpage()
grid.draw(p)
dev.off()


#featureplot
library(ArchR)
addArchRThreads(threads = 10)
setwd("~/CH/CH-ew/bed_total/")
proj2 <- readRDS('~/CH/CH-cross/Fig3/species/ew/ew_bigFigure_20230306.rds')


plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = 'main_ct', embedding = 'UMAP_peak_res1')


gene=getFeatures(proj2)
gene[1:5]
#gene=reshape2::colsplit(gene,':',names = c('c1','c2'))$c2
#foxa1
markerGenes<-c("evm.TU.Chr06.367","evm.TU.Chr11.200")
markerGenes <- c("evm.TU.Chr05.2330","evm.TU.Chr10.903")#not clear
markerGenes <- c("evm.TU.Chr03.779","evm.TU.Chr08.1311","evm.TU.Chr10.778" )

markerGenes=tmp
markerGenes=markerGenes[markerGenes%in%gene]
markerGenes

proj2=addImputeWeights(proj2,reducedDims = 'IterativeLSI_peak_res1')
p <- plotEmbedding(
  ArchRProj = proj2, 
  colorBy = "GeneScoreMatrix", 
  name = markerGenes, 
  embedding = 'UMAP_peak_res1',
  quantCut = c(0.01, 0.95),
  imputeWeights =getImputeWeights(proj2)
)
p

markersGS <- getMarkerFeatures(
  ArchRProj = proj2,
  useMatrix = "GeneScoreMatrix",
  groupBy = 'main_ct',
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
ew_markersGS=markersGS
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.05 & Log2FC >=0.5")
markerList@listData$Erythrocruorin$name[1:10]
ew_markerList=markerList

ea_ery=ew_markerList@listData$Erythrocruorin%>%as.data.frame()

g=intersect(gene_pair[gene_pair$V1%in%gene,]$V2%>%unique(),hu_ea_orth)
g[order(g)]
markersGS <- getMarkerFeatures(
  ArchRProj = proj2,
  useMatrix = "GeneScoreMatrix",
  groupBy = 'main_ct',
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
ew_markersGS=markersGS

#human
library(ArchR)
addArchRThreads(10)
setwd("~/CH/CH-cross/cluster/Renbing/bed/")
proj2 <- readRDS('~/CH/CH-cross/Fig3/species/intergrate_hg38/hg38_intergrate_bigFigure_0311.rds')

proj2=addImputeWeights(proj2,reducedDims = 'IterativeLSI_peak_res1')
markerGenes=c('TRIB1','SPRY2','CYP2C8','ABCA1','SLC2A2')
markerGenes=c('F2','F11')
p <- plotEmbedding(
  ArchRProj = proj2, 
  colorBy = "GeneScoreMatrix", 
  name = markerGenes, 
  embedding = 'UMAP_peak_res1',
  quantCut = c(0.01, 0.95),
  imputeWeights =getImputeWeights(proj2)
)
p

ea_gene.use1=ea_gene.use[ea_gene.use%in%ea_ery$name]
hu_gene.use1=hu_gene.use[hu_gene.use%in%hu_hepa$name]

gene_pair=fread('~/CH/CH-cross/Fig3/SAMap/maps/huea/ea_to_hu.txt')
gene_pair=gene_pair[gene_pair$V3>=30,]

hu_ea_orth1=gene_pair[gene_pair$V1%in%ea_gene.use1,]$V2
hu_ea_orth1=hu_ea_orth1[!duplicated(hu_ea_orth1)]

#change color
library(ArchR)
addArchRThreads(10)
setwd("~/CH/CH-cross/cluster/Renbing/bed/")
proj2 <- readRDS('~/CH/CH-cross/Fig3/species/intergrate_hg38/hg38_intergrate_bigFigure_0311.rds')
proj2=addImputeWeights(proj2,reducedDims = 'IterativeLSI_peak_res1')
ArchRProj=proj2

colorMat <- getMatrixFromProject(ArchRProj,useMatrix = 'GeneScoreMatrix')
colorMat=colorMat@assays@data$GeneScoreMatrix
rownames(colorMat)=getFeatures(ArchRProj)

colorMat=colorMat[c('F2','ABCA1','F11'),]

source('~/CH/CH-cross/Fig4/hepatocyte/ea_hu/util.R')


colorMat <- imputeMatrix(mat = as.matrix(colorMat), 
                         imputeWeights = getImputeWeights(ArchRProj), logFile =  createLogFile("plotEmbedding"))

df <- getEmbedding(ArchRProj, embedding = 'UMAP_peak_res1', returnDF = TRUE)
colorMat <- colorMat[, rownames(df), drop = FALSE]

colorMat
score1=colorMat[1,]
proj2$score1=score1
score2=colorMat[2,]
proj2$score2=score2


mydata <- proj2@embeddings$UMAP_peak_res1@listData$df %>% as.data.frame() %>% cbind(score=proj2$score2)
colnames(mydata) <- c('UMAP1','UMAP2','lineage')
mydata$plot <- 0
idx = which(proj2$score2>=2.5)
mydata$plot[idx] <- 1
table(mydata$plot)
# assign colors
mypalette <- data.frame(plot = c(0,1),tag = c('#DCDCDC','#E41A1C'))
mydata$color <- mypalette$tag[match(mydata$plot,mypalette$plot)]
table(mydata$color)
# reorder data.frame, in order to plot red points in the last step (important!)
mydata$plot <- as.factor(mydata$plot)
mydata <- mydata[order(mydata$plot),]
p <- ggplot(mydata, aes(x = UMAP1, y = UMAP2)) +
  geom_point(stat = "identity",
             aes(color = plot),
             color = mydata$color, size = 0.2, alpha = 0.5) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line.x = element_line(size = 0.5),
        axis.line.y = element_line(size = 0.5)) +
  theme(plot.title = element_text(size = 20, hjust = 0.5))  # main title aesthetic

p
