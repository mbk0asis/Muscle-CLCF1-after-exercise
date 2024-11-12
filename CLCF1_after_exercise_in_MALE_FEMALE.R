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
library(tidyverse)
library(magrittr)

#########################################################################
setwd("/home/bio0/BIO6/DrYang_Muscle_Aging")

# Download GDS file
gse <- getGEO('GSE28422', GSEMatrix=TRUE)
gse <- gse[[1]]

#########################################################################
# To get phenoGender data (Samples)
p_data <- pData(gse)
p_data <- rownames_to_column(p_data, "Sample") %>% 
  dplyr::select(Sample, Age=characteristics_ch1, Gender=characteristics_ch1.1, Trained=characteristics_ch1.2, TimePoint=characteristics_ch1.3)
head(p_data)

p_data$Age <- gsub("age: ","",p_data$Age)
p_data$Gender <- gsub("gender: ","",p_data$Gender)
p_data$Trained <- gsub("training state: ","",p_data$Trained)
p_data$TimePoint <- gsub("time point: ","",p_data$TimePoint)
p_data$TimePoint <- gsub("4hrs post exercise","post_ex",p_data$TimePoint)
p_data %<>% mutate(group = paste(Age,Gender,Trained,TimePoint,sep="_"))

tail(p_data)

p_data %>% count(Age, Gender, Trained, TimePoint) %>% mutate(grp = paste(Trained,TimePoint,sep="_")) %>%
  ggplot(aes(x=Age, y=n, fill=grp)) + theme_bw() +
  geom_bar(stat = "identity") +
  geom_text(aes(label = n), position = position_stack(vjust=0.5), size=7) +
  facet_wrap(~Gender) + labs(x="", y="", title = "Number of Samples") +
  theme(axis.text.x=element_text(size=15, angle=90, vjust=0.5, hjust=1),
        strip.text = element_text(size=12),
        panel.grid = element_blank())

# To get feature data (Genes)
# fData(gse)

#########################################################################
# To generate a expression dataset
exprs_data <- exprs(gse)
feature_data <- fData(gse)
gene_symbols <- feature_data$`Gene Symbol`
names(gene_symbols) <- rownames(feature_data)
rownames(exprs_data) <- gene_symbols[rownames(exprs_data)]
exprs_data_df <- as.data.frame(exprs_data)
exprs_data_df$Symbol <- rownames(exprs_data_df)
aggregated_data <- exprs_data_df %>% group_by(Symbol) %>% summarize(across(everything(), mean, na.rm = TRUE))
agg_data <- column_to_rownames(aggregated_data, "Symbol")
agg_data %<>% dplyr::select(p_data$Sample)
head(agg_data)

#
# PCA
pca.data <- t(agg_data)
pca <- prcomp(pca.data, scale=TRUE)
pca_data <- p_data %>% mutate(group = paste(Age, Gender, Trained, sep="_")) %>% 
  inner_join(rownames_to_column(data.frame(pca$x)), by=c("Sample"="rowname")) 
#dplyr::select(Sample, group, PC1, PC2)
head(pca_data)
ggplot(pca_data, aes(x=PC1, y=PC2, color=TimePoint)) + theme_bw() +
  geom_point() + stat_ellipse() +
  coord_fixed() +
  facet_wrap(~group, ncol=4) +
  theme(panel.grid = element_blank())


#install.packages("rstatix")
library(rstatix)
library(ggpubr)
clcf1 <- agg_data[grep("CLCF1", rownames(agg_data)), ]
clcf1 <- merge(t(clcf1), p_data, by.x=0, by.y="Sample") %>% mutate(group2 = paste(Age, TimePoint, sep="_"))
clcf1$group <- factor(clcf1$group, levels = unique(clcf1$group)[c(1:4,9:12,5:8,13:16)])
clcf1$group2 <- factor(clcf1$group2, levels = unique(clcf1$group2)[c(1,3,2,4)])
clcf1$Trained <- factor(clcf1$Trained, levels = unique(clcf1$Trained))
clcf1$Age <- factor(clcf1$Age, levels = c("Young","Old"))
head(clcf1)
ggplot(clcf1, aes(x=group, y=CLCF1, fill=group2)) + theme_bw() +
  geom_point(position = position_jitterdodge( jitter.width=0.1, dodge.width = 0.75)) +
  geom_boxplot(outlier.shape = NA, alpha=0.6) + xlab("") +
  stat_summary(position = position_dodge(width = 0.75), geom="point", color="red", pch=18, size=3) +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1, size=12), panel.grid.minor.x = element_blank(),
        strip.background = element_blank(), strip.text = element_text(size=12), legend.position = "none") +
  stat_compare_means(aes(group = TimePoint), 
                     label = "p.signif", method = "wilcox.test", label.y.npc=0.85, color="blue", size=15, hide.ns = TRUE) +
  scale_fill_manual(values = c("green","green4","yellow","brown")) +
  facet_wrap(~Gender, scales = "free")

stat.test <- compare_means(CLCF1 ~ group, data=clcf1, method = "wilcox")
stat.test %>% dplyr::filter(p.signif %in% "*") %>% View()

# Differential expression analysis
######################################
######################################
######################################
pval=0.05            
log2FoldChange=log2(1.5)   # "log2(1)" or "0" for not to consider fold-changes

######################################
######################################
######################################
design <- model.matrix(~0 + group, data=p_data)
head(design)

control <- "groupYoung_fast_Untrained"
sample <- "groupYoung_fast_Trained"
cont <- paste0(sample,"-",control)
c <- makeContrasts(cont, levels=design) 

# DE analysis
fit <- lmFit(agg_data, design) 
fit2<- contrasts.fit(fit,c) 
fit2<- eBayes(fit2) 
n <- nrow(agg_data)
top <- topTable(fit2, coef=1, n=n, adjust = "fdr")
top2 <- top %>% mutate(grp = case_when(logFC > log2FoldChange & P.Value < pval ~ "sigUp", logFC < log2FoldChange & P.Value < pval ~ "sigDn", TRUE ~ "nonSig")) %>%
  mutate(logFC2 = case_when(logFC >= 10 ~ 10, logFC <= -10 ~ -10, TRUE~logFC))

ggplot(top2, aes(x=logFC2, y=-log10(P.Value), color=grp)) + theme_bw() +
  geom_point()

#########################################################################
nrow1 <- nrow(subset(top, P.Value <= pval & logFC >= log2FoldChange))
nrow2 <- nrow(subset(top, P.Value <= pval & logFC <= -log2FoldChange))

#agg_data$logFC[abs(agg_data) >= 4] <- 4                            # replace (foldChange > 4) to 4
top$P.Value[top$P.Value <= 1e-15] <- 1e-15                # replace (pvalue < 1e-15) to 1e-15

with(agg_data, plot(logFC,-log10(P.Value), pch=19,cex=1,
                    col=rgb(.5,.5,.5,.5), 
                    main=paste("Whole genes\npval < ",pval,",log2FC > ",round(log2FoldChange,3),"\n",control," (",nrow2,")   < >   ",sample," (",nrow1,")",sep="")))
with(subset(agg_data, P.Value<=pval & logFC>=log2FoldChange),
     points(logFC,-log10(P.Value),pch=19,cex=.9,col=rgb(1,.2,.2,.5)))
with(subset(agg_data, P.Value<=pval & logFC<=-log2FoldChange),
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

#res_sig <- subset(agg_data, P.Value <= pval & abs(logFC) >= log2FoldChange)
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
  geom_hline(yintercept = c(-logfc,0,logfc), lineGender="dotted", size =.5, color ="blue") +
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


