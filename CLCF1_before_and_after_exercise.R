#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages(c("BiocManager","dplyr","gplots"))
#BiocManager::install(c("GEOquery","limma"))
#BiocManager::install(c("limma"))

#########################################################################
# load libraries


library(Biobase)
library(GEOquery)
library(limma)
library(dplyr)
library(gplots)
library(reshape2)
library(ggplot2)
#########################################################################
setwd("/home/bio0/BIO6/DrYang_Muscle_Aging")

# Download GDS file
gds5218 <- getGEO('GDS5218', destdir=".")

#Or, open an existing GDS file 
#gds5218 <- getGEO(filename='GDS5218.soft.gz')
#########################################################################
# 
# Meta(gds5218)$channel_count
# Meta(gds5218)$description
# Meta(gds5218)$feature_count
# Meta(gds5218)$platform
# Meta(gds5218)$sample_count
# Meta(gds5218)$sample_organism
# Meta(gds5218)$sample_type
# Meta(gds5218)$title
# Meta(gds5218)$type
#
#########################################################################
# To generate a expression dataset
eset <- GDS2eSet(gds5218, do.log2=TRUE)
eset
e <- exprs(eset)
head(e[,1:5])
#write.csv(e,"All.Expression.csv")

# to see the sample phenotypes
#show(pData(phenoData(eset))[,1:3])
boxplot(e, col=c(rep("grey30",15),rep("grey",16),
                 rep("grey50",15),rep("lightgrey",16),
                 rep("grey30",12),rep("grey",12),
                 rep("grey50",12),rep("lightgrey",12)))


# PCA
pca.data <- e

grps <- c(rep("yBase",15),rep("y1st",16),
          rep("yB4Final",15),rep("yFinal",16),
          rep("oBase",12),rep("o1st",12),
          rep("oB4Final",12),rep("oFinal",12))

grpcol <- c(rep(rgb(0,0,1,0.6),15),rep(rgb(.2,.5,1,0.6),16),
            rep(rgb(1,1,1,0),15),rep(rgb(0,1,1,0.5),16),
            rep(rgb(1,0,0,0.6),12),rep(rgb(1,.5,0,0.6),12),
            rep(rgb(1,1,1,0),12),rep(rgb(1,.7,.3,0.6),12))

pca.data <- t(pca.data)

colnames(pca.data) <- paste(grps, colnames(pca.data), sep="-")

pca <- prcomp(pca.data, scale=TRUE)
summary(pca)
plot(pca$x[,1], pca$x[,2], xlab="PCA1", ylab="PCA2", main="PCA", type="p", pch=19, cex=2, col=grpcol)
#text(pca$x[,1], pca$x[,2], rownames(pca.data), cex=0.75)

coord <- data.frame(pca$x[,1:2])
coord$group <- grps
coord$group <- factor(coord$group, levels=c("yBase","y1st","yB4Final","yFinal",
                                            "oBase","o1st","oB4Final","oFinal"))
ggplot(coord, aes(x=PC1, y=PC2, color=group)) + theme_bw() +
  geom_point(size=5, alpha=.5) + coord_fixed()


#########################################################################
design <- cbind(yBase=c(rep(1,15),rep(0,95)),
                y1st=c(rep(0,15), rep(1,16), rep(0,79)),
                yPL=c(rep(0,31), rep(1,15), rep(0,64)),
                yFinal = c(rep(0,46), rep(1,16), rep(0,48)),
                oBase=c(rep(0,62), rep(1,12), rep(0,36)),
                o1st=c(rep(0,74), rep(1,12), rep(0,24)),
                oPL=c(rep(0,86), rep(1,12), rep(0,12)),
                oFinal = c(rep(0,98), rep(1,12)))  

# Differential expression analysis
######################################
######################################
######################################
pval=0.05            
log2FoldChange=log2(1.5)   # "log2(1)" or "0" for not to consider fold-changes
                        
control <- "yBase"
sample <- "yFinal"
######################################
######################################
######################################

cont <- paste(sample,"-",control, sep="")
c <- makeContrasts(cont, levels=design) 

# DE analysis
fit <- lmFit(eset,design) 
fit2<- contrasts.fit(fit,c) 
fit2<- eBayes(fit2) 

n <- length(featureNames(eset))
top <- topTable(fit2, coef=1, n=n, adjust = "fdr")
res <- top[,c(3,10,22:26)]

# aggregate to gene symbols by means
aggdata <-aggregate(res, by=list(res$Gene.symbol),
                    FUN=mean, na.rm=TRUE)
aggdata2 <- aggdata[,-c(1,2,3)]
row.names(aggdata2) <- aggdata[,1]
head(aggdata2)
dim(aggdata2)
#write.csv(aggdata2,paste("DE_whole_genes_",sample,"_vs_",control,".csv",sep=""))

#########################################################################
nrow1 <- nrow(subset(aggdata2, P.Value <= pval & logFC >= log2FoldChange))
nrow2 <- nrow(subset(aggdata2, P.Value <= pval & logFC <= -log2FoldChange))

#aggdata2$logFC[abs(aggdata2) >= 4] <- 4                            # replace (foldChange > 4) to 4
aggdata2$P.Value[aggdata2$P.Value <= 1e-15] <- 1e-15                # replace (pvalue < 1e-15) to 1e-15

with(aggdata2, plot(logFC,-log10(P.Value), pch=19,cex=1,
               col=rgb(.5,.5,.5,.5), 
               main=paste("Whole genes\npval < ",pval,",log2FC > ",round(log2FoldChange,3),"\n",control," (",nrow2,")   < >   ",sample," (",nrow1,")",sep="")))
with(subset(aggdata2, P.Value<=pval & logFC>=log2FoldChange),
     points(logFC,-log10(P.Value),pch=19,cex=.9,col=rgb(1,.2,.2,.5)))
with(subset(aggdata2, P.Value<=pval & logFC<=-log2FoldChange),
     points(logFC,-log10(P.Value),pch=19,cex=.9,col=rgb(0,0,1,.5)))
abline(h=c(-log10(pval)),col="grey30",lty=2)
abline(v=c(-log2FoldChange,log2FoldChange),col="grey30",lty=2)
abline(v=0,col="black")

#########################################################################

#test <- merge(annot,e, by=0)
#dim(test)
#head(test[,1:5])
#aggdata3 <-aggregate(test, by=list(test$Gene.symbol),
#                     FUN=mean, na.rm=TRUE)
#head(aggdata3[,1:5])
#write.csv(annot,"annotation.csv")

#res_sig <- subset(aggdata2, P.Value <= pval & abs(logFC) >= log2FoldChange)
#res_sig <- res_sig[order(-res_sig$logFC),]     
#head(res_sig)
#dim(res_sig)
#write.csv(res_sig,paste("sigDEGs_",sample,"_vs_",control,".csv",sep=""))
#########################################################################
# check if "CLCF1" is included
clcf <- aggdata[grepl("CLCF1",aggdata[,1]),]
#ccl2 <- aggdata[grepl("^CCL2$",aggdata[,1]),]

#res_sig[grepl("CLCF1",res_sig[,1]),]

#########################################################################

#write.csv(res_sig,file=paste("limma_sigDEGs_",sample,"_vs_",control,".csv",sep=""), quote=F, row.names=F)

#########################################################################

secretory.proteins <- read.csv("secretory.proteins.csv", header=F)
res_secretory.proteins <- merge(secretory.proteins, aggdata, by.x="V1", by.y="Group.1")
#write.csv(res_secretory.proteins,paste("DE_secretory.proteins_",sample,"_vs_",control,".csv",sep=""))
#
res_sig <- subset(res_secretory.proteins, P.Value <= pval & abs(logFC) >= log2FoldChange)
res_sig <- res_sig[order(-res_sig$logFC),]    

#write.csv(res_sig,paste("sigDEGs_secretoryProteins_",sample,"_vs_",control,".csv",sep=""))

#res_non_sig <- subset(res_secretory.proteins, P.Value > pval & abs(logFC) < log2FoldChange)
#write.csv(res_non_sig,"non_sigDEGs_o1st_oBase.csv")

##########################################################################
######################################
pval=0.01            
log2FoldChange=log2(1.5)   # "log2(1)" or "0" for not to consider fold-changes

nrow1 <- nrow(subset(res_secretory.proteins, P.Value <= pval & logFC >= log2FoldChange))
nrow2 <- nrow(subset(res_secretory.proteins, P.Value <= pval & logFC <= -log2FoldChange))

res_secretory.proteins$logFC[abs(res_secretory.proteins$logFC) >= 4] <- 4                         # replace (foldChange > 4) to 4
res_secretory.proteins$P.Value[res_secretory.proteins$P.Value <= 1e-8] <- 1e-8                # replace (pvalue < 1e-15) to 1e-15

# Volcano plot
with(res_secretory.proteins, plot(logFC,-log10(P.Value), pch=19,cex=1.3,
               col=rgb(.5,.5,.5,.5), 
               main=paste("Secretory proteins\npval < ",pval,",   log2FC > ",round(log2FoldChange,3),"\n",control," (",nrow2,")   < >   ",sample," (",nrow1,")",sep="")))
with(subset(res_secretory.proteins, P.Value<=pval & logFC>=log2FoldChange),
     points(logFC,-log10(P.Value),pch=19,cex=1.2,col=rgb(1,.2,.2,.6)))
with(subset(res_secretory.proteins, P.Value<=pval & logFC<=-log2FoldChange),
     points(logFC,-log10(P.Value),pch=19,cex=1.2,col=rgb(0,0,1,.6)))
abline(h=c(-log10(pval)),col="grey30",lty=2)
abline(v=c(-log2FoldChange,log2FoldChange),col="grey30",lty=2)
abline(v=0,col="black")

clcf <- aggdata[grepl("CLCF1",aggdata[,1]),]
points(clcf$logFC,-log10(clcf$P.Value),col="green3", pch=19)
text(clcf$logFC,-log10(clcf$P.Value), "CLCF1", col="green3", pos=4, cex = 1.3)

ANGPTL2 <- aggdata[grepl("ANGPTL2",aggdata[,1]),]
points(ANGPTL2$logFC,-log10(ANGPTL2$P.Value),col="green3", pch=19)
text(ANGPTL2$logFC,-log10(ANGPTL2$P.Value), "ANGPTL2", col="green3", pos=3, cex = 1.3)

CCL2 <- aggdata[grepl("CCL2$",aggdata[,1]),]
points(CCL2$logFC,-log10(CCL2$P.Value),col="green3", pch=19)
text(CCL2$logFC,-log10(CCL2$P.Value), "CCL2", col="green3", pos=4, cex = 1.3)

FJX1 <- aggdata[grepl("FJX1$",aggdata[,1]),]
points(FJX1$logFC,-log10(FJX1$P.Value),col="green3", pch=19)
text(FJX1$logFC,-log10(FJX1$P.Value), "FJX1", col="green3", pos=3, cex = 1.3)

LGI1 <- aggdata[grepl("LGI1$",aggdata[,1]),]
points(LGI1$logFC,-log10(LGI1$P.Value),col="green3", pch=19)
text(LGI1$logFC,-log10(LGI1$P.Value), "LGI1", col="green3", pos=2, cex = 1.3)
##########################################################################
##########################################################################
##########################################################################
exp_secretory.proteins <- read.csv("Expression.matrix_secretory.proteins.csv", header = T, row.names = 1)
dim(exp_secretory.proteins)
head(exp_secretory.proteins)

mean_yBase <- rowMeans(exp_secretory.proteins[,1:15])
mean_y1st <- rowMeans(exp_secretory.proteins[,16:31])
mean_yFinal <- rowMeans(exp_secretory.proteins[,32:47])
mean_oBase <- rowMeans(exp_secretory.proteins[,48:59])
mean_o1st <- rowMeans(exp_secretory.proteins[,60:71])
mean_oFinal <- rowMeans(exp_secretory.proteins[,72:83])

y_mean <- data.frame(mean_yBase,mean_y1st,mean_yFinal)
o_mean <- data.frame(mean_oBase,mean_o1st,mean_oFinal)

#################################
fc_y <- y_mean - mean_yBase
fc_o <- o_mean - mean_oBase
fc <- data.frame(fc_y,fc_o)

#################################
all_mean <- data.frame(mean_yBase,mean_y1st,mean_yFinal,mean_oBase,mean_o1st,mean_oFinal)
fc <- all_mean - mean_yBase

head(fc)

##########################################################################
# k-means clustering
logfc = log2(1.5)
nclust = 20
set.seed(5)

cl <- kmeans(fc, nclust)
cluster <- cl$cluster

rNames <- rownames(fc)
df<-data.frame(fc, rNames, cluster) 
df[grepl("CLCF1",rownames(df)),7:8] # cluster number of CLCF1
nCl <- df[grepl("CLCF1",rownames(df)),8]
mdata <- melt(df, id=c("rNames","cluster")) 

table(cluster)
df[grepl("CLCF1",rownames(df)),7:8] # cluster number of CLCF1

jit <- position_jitter (width = .08, height = .08)
ggplot(mdata, aes(x=variable, y=value)) + theme_bw() + ylim(-1.5,1.5) + 
  geom_hline(yintercept = c(-logfc,0,logfc), linetype="dotted", size =.5, color ="blue") +
  geom_line(aes(group=rNames), alpha =.1, position = jit) +
  stat_summary(aes(group = cluster), fun.y = mean, geom = "line", color="red", size = .5) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_wrap(~ cluster) # plots divided by the clusters

keepSamples = (cluster==nCl) # CLCF1 containing cluster
fc_sel = fc[keepSamples, ]
head(fc_sel)
nGenes <- nrow(fc_sel)

exp_ind <- merge(fc_sel, exp_secretory.proteins, by=0)  
exp_ind2 <- exp_ind[,-c(1:7)]
row.names(exp_ind2) <- exp_ind[,1]
exp_ind2

fc_y <- exp_ind2[,1:47] - rowMeans(exp_ind2[,1:15])
fc_o <- exp_ind2[,48:83] - rowMeans(exp_ind2[,48:59])
fc_cluster <- data.frame(fc_y,fc_o)
head(fc_cluster)

nR <- which(rownames(fc_cluster) == "CLCF1") # cluster number of CLCF1
nBelow <- nGenes-nR

my_palette <- colorRampPalette(c("darkblue", "black", "yellow"))
breaks = unique(c(seq(-2,0,length=75), seq(0,2,length=75)))
heatmap.2(as.matrix(fc_cluster), density="none",trace="none",
          breaks=breaks, col=my_palette,
          Rowv=T, Colv=F, dendrogram="row", margins = c(8,8),
          main=paste0("Secretory Proteins (cluster = ",nCl,", ",nGenes," genes)\nlog2 (Rel. to Base)"),
          ColSideColors=c(rep("red",16), rep("green",16), rep("blue",15),
                          rep("orange",12), rep("darkgreen",12), rep("darkblue",12)),
          RowSideColors = c(rep("white",nR-1),rep("red",1),rep("white",nBelow)))


