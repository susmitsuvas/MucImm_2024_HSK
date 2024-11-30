library(Seurat)
library(SeuratDisk)
library(tidyverse)
library(Nebulosa)
library(devtools)
library(pbapply)
library(ComplexHeatmap)
library(ggpubr)

session_info()

#working directory
setwd("/tmp/scRNAseq/myProject")

#read saved integrated data
pso <- readRDS("savedData/Integrated_Conditions.rds")

#generate Elbow Plot
png("plots/qcPlots/ElbowPlot.png",width=2000,height=2000,res=300)
elbowPlot <- ElbowPlot(pso,ndim = 100) +
theme_light() +
theme(panel.grid = element_blank(),
panel.border = element_blank())
print(elbowPlot)
dev.off()

#generate integrated plot by group
png("plots/qcPlots/Grouped_TSNE.png",res=300,width=4000,height=2000)
print(DimPlot(pso,reduction="tsne",group.by="Condition",split.by="Condition"))
dev.off()

png("plots/qcPlots/Grouped_UMAP.png",res=300,width=4000,height=2000)
print(DimPlot(pso,reduction="umap",group.by="Condition",split.by="Condition"))
dev.off()

png("plots/qcPlots/Grouped_TSNE_clustered.png",res=300,width=4000,height=2000)
print(DimPlot(pso,reduction="tsne",split.by="Condition",label = TRUE, repel = TRUE))
dev.off()

png("plots/qcPlots/Grouped_UMAP_clustered.png",res=300,width=4000,height=2000)
print(DimPlot(pso,reduction="umap",split.by="Condition",label = TRUE, repel = TRUE))
dev.off()

# selected dim = 20 and resolution = 0.8 as the optimal SO for processing
pso <- readRDS("savedData/Permutation_Data_d20_r0.8.rds")
png("plots/qcPlots/d20_r0.8_tsne.png",res=300,width=2000,height=2000)
print(DimPlot(pso,reduction="tsne",label=TRUE,repel=TRUE))
dev.off()


png("plots/qcPlots/d20_r0.8_umap.png",res=300,width=2000,height=2000)
print(DimPlot(pso,reduction="umap",label=TRUE,repel=TRUE))
dev.off()


#cell proportions analysis
CP <- table(pso$orig.ident, factor(Idents(pso)))
CP <- CP/rowSums(CP)
CP <- reshape2::melt(CP)
colnames(CP) <- c('Sample', 'CT', 'Proportion')

CP$CT <- as.character(CP$CT)

png(paste0('plots/qcPlots/CellProportions.png'), width = 3500, height = 1500, res = 300)
ggplot(CP, aes(Proportion, Sample, fill = CT)) +
geom_bar(stat = 'identity') +
theme_light() +
theme(legend.position = 'bottom',
panel.grid = element_blank(),
panel.border = element_blank(),
strip.text = element_text(color = 'black')) +
labs(fill = 'Cell Type')
dev.off()

####Perform scProportion Test
library(scProportionTest)
#read saved integrated data
pso <- readRDS("savedData/Permutation_Data_d20_r0.8.rds")


prop_test <- sc_utils(pso)
prop_test <- permutation_test(
	prop_test, cluster_identity = "seurat_clusters",
	sample_1 = "SubClinical", sample_2 = "Clinical",
	sample_identity = "orig.ident"
)

permutation_plot(prop_test)

# change cell ids for cluster 10
pso.c10 <- subset(pso, subset = seurat_clusters == 10)

#Add cellid as Meta data
pso.c10 <- AddMetaData(pso.c10,rownames(pso.c10@meta.data),"cell.id")
#remove cells that are seemingly not part of cluster 10
pso.c10@reductions$tsne@cell.embeddings[c(order(as.numeric(pso.c10@reductions$tsne@cell.embeddings[,1]))),]

r10 <- rownames(pso.c10@reductions$tsne@cell.embeddings[c(order(as.numeric(pso.c10@reductions$tsne@cell.embeddings[,1]))),])
pso.c10 <- subset(pso.c10, subset = cell.id %in% r10[-c(351,352)])



#generate feature Plots
genes_comb <- list(c("Krt14","Ifitm3","Gpha2","Slurp1","Krt12"),c("Gja1","Itgb4","Slurp1","Krt12","Col17a1"), c("Dsg1a","Dsp","Dsc2"), c("Cldn23","Muc4"),c("Cd74","H2-Ab1"),c("Hbb-bs","Hba-a2"),c("Ccl5","Cd3d","Cd3g","Cd4","Cd8a","Gzma","Gzmb","Icos","Ctla2a","Trbc1","Ifng","Ms4a4b"))

genes_all <- unique(unlist(genes_comb))


i=1
for(gene in genes_all){
	pngfile <- paste0("plots/featurePlots/",gene,"_tsne.png")
	png(pngfile,res=300,width=2000,height=2000)
	print(FeaturePlot(pso,features=gene,reduction="tsne",min.cutoff = "q9",label=TRUE,repel=TRUE,cols=c("grey","red")))
	dev.off()
	i = i + 1
}

i=1
for(gene in genes_all){
	pngfile <- paste0("plots/featurePlots/",gene,"_umap.png")
	png(pngfile,res=300,width=2000,height=2000)
	print(FeaturePlot(pso,features=gene,reduction="umap",min.cutoff = "q9",label=TRUE,repel=TRUE,cols=c("grey","red")))
	dev.off()
	i = i + 1
}

####Generate Cluster Markers Table
DefaultAssay(pso) <- "RNA"
sample.markers <- FindAllMarkers(pso, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(sample.markers, file = "topMarkers/UnAnnotated_Marker_Table_d20_r0.8.csv")


genes_violin <- c("Ifng","Tgfbi", "Dapl1", "Mgarp", "Ifi27l2a","Fabp5","Rpl32")


for(geneID in genes_violin){

	#geneID <- genes_violin[2]
    df <- data.frame(I = Idents(pso.c10), E = pso.c10@assays$RNA@data[geneID,], G = pso.c10$orig.ident)
    newFile <- paste0('plots/featurePlots/', geneID, '_vln.png')
    FC <- sapply(split(df$E, df$G), mean)
    if(FC[1] > FC[2]){
        colorPattern <- c('red', 'blue')
    } else {
        colorPattern <- c('blue', 'red')
    }
    png(newFile, width = 350, height = 750, res = 300)
    P <- ggplot(df, aes(G, E, fill = G, color = G)) +
    geom_violin() +
    theme_classic() +
    ylab('Normalized Expression') +
    xlab('Genotype') +
    theme(legend.position = 'None', axis.text.x=element_text(angle=90, vjust = 0.5, hjust = 1)) +
    scale_fill_manual(values = colorPattern) + 
    scale_color_manual(values = colorPattern) + 
    stat_compare_means(label = 'p.signif', label.x = 1.35, label.y = max(df$E)) +
    ylim(c(0, max(df$E)*1.1)) +
    labs(title = geneID) 
    print(P)
    dev.off()

}


for(geneID in genes_violin){

	#geneID <- genes_violin[1]
    df <- data.frame(I = Idents(pso.c10), E = pso.c10@assays$RNA@data[geneID,], G = pso.c10$orig.ident)
    newFile <- paste0('plots/featurePlots/', geneID, '_box.png')
    FC <- sapply(split(df$E, df$G), mean)
    if(FC[1] > FC[2]){
        colorPattern <- c('red', 'blue')
    } else {
        colorPattern <- c('blue', 'red')
    }
    png(newFile, width = 350, height = 750, res = 300)
    P <- ggplot(df, aes(G, E, fill = G, color = G)) +
    geom_boxplot(aes(lower=mean-sd,upper=mean+sd,middle=mean,ymax=mean+3*sd),stat="identity") +
    theme_classic() +
    ylab('Normalized Expression') +
    xlab('Genotype') +
    theme(legend.position = 'None', axis.text.x=element_text(angle=90, vjust = 0.5, hjust = 1)) +
    scale_fill_manual(values = colorPattern) + 
    scale_color_manual(values = colorPattern) + 
    stat_compare_means(label = 'p.signif', label.x = 1.35, label.y = max(df$E)) +
    ylim(c(0, max(df$E)*1.1)) +
    labs(title = geneID) 
    print(P)
    dev.off()

}


geneID <- genes_violin[7]
df <- data.frame(I = Idents(pso.c10), E = pso.c10@assays$RNA@data[geneID,], G = pso.c10$orig.ident)
P <- ggplot(df, aes(G, E, fill = G, color = G)) +
    geom_violin() +
    theme_classic() +
    ylab('Normalized Expression') +
    xlab('Genotype') +
    theme(legend.position = 'None', axis.text.x=element_text(angle=90, vjust = 0.5, hjust = 1)) +
	stat_compare_means()
P
