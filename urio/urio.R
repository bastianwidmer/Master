####### Final exam ######
library(affy)
openVignette()
library(stringr)
library(affyPLM)
library(limma)
library(cluster)
library(xtable)
library(gplots)

GSE1561 <- ReadAffy(celfile.path = "~/Desktop/EPFL/2016-2017 2ème semestre master/Projet-Meta-analysis/project_data/GSE1561_RAW/GSE1561",compress = TRUE)
names <- str_sub(sampleNames(GSE1561),1,8)
sampleNames(GSE1561) <- names
rm(names)
# 49 chips

## Quality assessment ##

# Histogram, boxplot and images #
GSE1561.pm <- pm(GSE1561)
rm(GSE1561.pm)
hist(GSE1561,main="Histogram of data GSE1561")
par(mai=c(1.2,1.2,0.1,0.1))
boxplot(GSE1561,xlab="Array", ylab = expression(paste(" ", log[2], " PM signal intensity", sep = "")), col="gray")

# PLM fitting #
Pset1561 <- fitPLM(GSE1561, model.param = list("Huber"))
par(mfrow=c(3,3),pty="s",mai=c(0.1,0.1,0.3,0.1))
for (i in 10:36) {
  image(Pset1561,which=i)
}
rm(i)
par(mfrow=c(1,1), pty="s", mai=c(0.2,0.2,0.2,0.2))
image(Pset1561, which = "GSM26870")
image(Pset1561, which = "GSM26880")
image(Pset1561, which = "GSM26885")
image(Pset1561, which = "GSM26890")
image(Pset1561, which = "GSM26901")
image(Pset1561, which = "GSM26903")
image(Pset1561, which = "GSM26904")
image(Pset1561, which = "GSM26906")
image(Pset1561, which = "GSM26908")

# RLE #
par(mfrow=c(1,1), mai=c(1.2,1.2,0.1,0.1))
RLE(Pset1561, xlab = "Array", ylab = "RLE",col="gray")

# NUSE #
par(mfrow=c(1,1), mai=c(1.2,1.2,0.1,0.1))
couleur <- rep("gray", 49)
couleur[which((NUSE(Pset1561, type="stat")[1,]+(0.5*NUSE(Pset1561, type="stat")[2,])) > 1.05)] <- "orange"
NUSE(Pset1561, xlab="Array",ylab="NUSE", col=couleur, ylim=c(0.93,1.15)) # tous en dessous de 1.05
abline(h=1.05,col="red", lty="dotdash", lwd=2)
legend("topright","Arrays with 75% NUSE values higher than 1.05", fill="orange")
rm(couleur)
colnames(GSE1561[,which((NUSE(Pset1561, type="stat")[1,]+(0.5*NUSE(Pset1561, type="stat")[2,])) > 1.05)])

# Normalization #
GSE1561.rma <- rma(GSE1561)
GSE1561.expr <- exprs(GSE1561.rma)

# initialization of the ER data #

GSE1561.clindata <- read.table("~/Desktop/EPFL/2016-2017 2ème semestre master/Projet-Meta-analysis/project_data/GSE1561.txt", skip = 28, fill = TRUE, sep = "")
GSE1561.clin <- GSE1561.clindata[1,]
row.names(GSE1561.clin) <- GSE1561.clin[1,1]
GSE1561.clin <- t(GSE1561.clin[,-1])
row.names(GSE1561.clin) <- c(1:49)
rm(GSE1561.clindata)
GSE1561.ER <- str_sub(GSE1561.clin,6,7)
for (i in 1:49){
  if (GSE1561.ER[i] == "En"){
    GSE1561.ER[i] <- 0
  }
  else {
    GSE1561.ER[i] <- 1
  }
}
GSE1561.ER <- as.numeric(GSE1561.ER)
rm(i)
rm(GSE1561.clin)

# DE genes #
design1561 <- cbind(rep(1,49),(rep(1,49)-GSE1561.ER))
colnames(design1561) <- c("ER+","ER- vs ER+")
rownames(design1561) <- sampleNames(GSE1561)
xtable(design1561, digits = 0)

fit1561 <- lmFit(GSE1561.rma, design1561)
fit1561 <- eBayes(fit1561)
table1561 <- topTable(fit1561, coef = "ER- vs ER+", adjust="BH", number=50, p.value=0.05)
table1561 <- table1561[,c(1,3,4,5)]
xtable(table1561, caption = "Top of fifty of genes differentially expressed for GSE1561.", label = "topTable1561", align = c("l","c","c","c","c"), display = c("s","f","f","e","e"))

results.BH1561 <- decideTests(fit1561, adjust.method = "BH", p.value = 0.05)
vennDiagram(results.BH1561, main="Venn diagram BH") #rejette le plus
pv1561 <- fit1561$p.value[,2] 
pv.BH1561 <- p.adjust(pv1561, method="BH") #rejects 3647

# Volcano plot #
M1561 <- fit1561$coefficients[,2]
t1561 <- abs(fit1561$t[,2])
par(mfrow=c(1,1), mai=c(1.2,1.2,0.1,0.1))
plot(M1561,t1561,pch=".",ylab="| mod t |", xlab = expression(paste(hat(beta)[2][k])))
rej.BH1561 <- which(pv.BH1561<= 0.05)
points(M1561[rej.BH1561],t1561[rej.BH1561],pch=".",col="green")
abline(h=min(t1561[rej.BH1561]),col="red",lty="dotted",lwd=2)
abline(v=c(-1,1), col="red", lty="dotted", lwd=2)
rm(rej.BH1561)
rm(pv1561)
rm(pv.BH1561)
rm(results.BH1561)
rm(M1561)

###################
### Clustering ####
###################
index <- which(GSE1561.ER==0) 
erneg <- GSE1561[,erneg]
erneg.rma <- GSE1561.rma[,index]

variance <- matrix(nrow=22283,ncol=1)
for(i in 1:22283){
  variance[i,1]=var(as.matrix(erneg.rma)[i,1:3])
}
rm(i)
var.order <- order(variance, decreasing = TRUE)
top100 <- var.order[1:100]
GSE1561.top100 <- erneg.rma[top100,]

## hierarchical clustering ##

# correlation distance #
clust.cor.complete <- hclust(as.dist(1-cor(as.matrix(GSE1561.top100))), method = "complete") 
plot(clust.cor.complete, xlab = "Dissimilarity (1-correlation)", sub = "Complete", main = "")
plot(silhouette(cutree(clust.cor.complete, k=2), dist = as.dist(1-cor(as.matrix(GSE1561.top100)))), main = "") #average 0.52

par(mfrow=c(1,1), mai=c(1.5,1.2,0.1,0.1))
clust.cor.ward <- hclust(as.dist(1-cor(as.matrix(GSE1561.top100))), method = "ward.D") 
plot(clust.cor.ward, xlab = "Dissimilarity (1-correlation)", sub = "Ward's method", main = "")
par(mfrow=c(1,1), mai=c(1.5,0.8,0.4,0.4))
plot(silhouette(cutree(clust.cor.ward, k=2), dist = as.dist(1-cor(as.matrix(GSE1561.top100)))), main = "") #average 0.52

clust.cor.average <- hclust(as.dist(1-cor(as.matrix(GSE1561.top100))), method = "average") 
plot(clust.cor.average, xlab = "Dissimilarity (1-correlation)", sub = "Average", main = "")
plot(silhouette(cutree(clust.cor.average, k=2), dist = as.dist(1-cor(as.matrix(GSE1561.top100)))), main = "") #average 0.52

# k-means #
dissE <- daisy(t(as.matrix(GSE1561.top100)), metric = "euclidean")
dist.e <- dist(t(as.matrix(GSE1561.top100)), method = "euclidean")

clust.km.2 <- kmeans(t(as.matrix(GSE1561.top100)),centers=2) 
clust.km.3 <- kmeans(t(as.matrix(GSE1561.top100)),centers=3)
plot(clust.km.2$cluster) 
plot(t(as.matrix(GSE1561.top100)), col=clust.km.2$cluster)
points(clust.km.2$centers, col = 1:2, pch=8, cex=2)
plot(silhouette(clust.km.2$cluster, dist = dist.e)) #average 0.34
plot(silhouette(clust.km.3$cluster, dist = dist.e)) #average 0.32

## PAM ##
pam.2 <- pam(t(as.matrix(GSE1561.top100)), diss = FALSE,2)
plot(pam.2$clustering)
summary(pam.2)
pam.2$silinfo # average 0.34

## Heatmap ##

samples.corr.ward <- as.dendrogram(hclust(as.dist(1-cor(as.matrix(GSE1561.top100))), method = "ward.D"))
genes.corr.ward <-  hclust(as.dist(1-cor(t(as.matrix(GSE1561.top100)))), method = "ward.D")
par(mfrow=c(1,1), mai=c(1.5,1.2,0.2,0.1))
heatmap.2(as.matrix(GSE1561.top100),genes.corr.ward,samples.corr.ward,col=topo.colors(32), key.title = "Color key", labRow = "", margins = c(5,0.1), dendrogram = "column")

par(mfrow=c(1,1), mai=c(1.5,1.2,0.2,0.1))
plot(genes.corr.ward, main="",xlab = "Dissimilarity (1-correlation)", sub = "Ward's method")
plot(samples.corr.ward)

genes.corr.ward$labels
x <- matrix(genes.corr.ward$labels, ncol=2,nrow=50)
print(xtable(x, caption = "List of 100 genes with largest variance.", label = "topTable"), include.rownames = FALSE)
