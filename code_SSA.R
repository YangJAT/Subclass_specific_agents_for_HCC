library(ggplot2)
library(ggpubr)
library(curl)
library("NMF")
library(sva)
library(CancerSubtypes)
library(clusterProfiler)
library(enrichplot)
library(gridExtra)
library(ISOpureR)
library(pRRophetic)
library(MCPcounter)
set.seed(123)

# removing batch effect
# combat function
data <- read.csv("input.csv",row.names = 1, stringsAsFactors = F)
data <- as.matrix(data)
batchType <- c(rep(1,90),rep(2,90))
data <- ComBat(data, batchType, par.prior=TRUE)
data <- rbind(geneNames=colnames(data),data)
data <- t(data[-1,])
write.table(data, "normlize.txt", sep="\t", quote=F, row.names=T,col.names=NA)
# PCA
data <- read.table("normlize.txt",sep="\t",header=T,check.names=F,row.names=1) 
data.pca <- prcomp(data, scale. = TRUE) 
group <- c(rep("cohort1",90),rep("cohort2",90))
pcaPredict <- predict(data.pca)
PCA <- data.frame(PCA1 = pcaPredict[,1], PCA2 = pcaPredict[,2],group=group)
PCA.mean=aggregate(PCA[,1:2],list(group=PCA$group),mean)

veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) {
    theta <- (0:npoints) * 2 * pi/npoints
    Circle <- cbind(cos(theta), sin(theta))
    t(center + scale * t(Circle %*% chol(cov)))
}
df_ell <- data.frame()
for(g in levels(PCA$group)){
df_ell <- rbind(df_ell, cbind(as.data.frame(with(PCA[PCA$group==g,],
                  veganCovEllipse(cov.wt(cbind(PCA1,PCA2),
                  wt=rep(1/length(PCA1),length(PCA1)))$cov,
                  center=c(mean(PCA1),mean(PCA2))))),group=g))
}

pdf(file="PCA.pdf",height=6,width=7)
ggplot(data = PCA, aes(PCA1, PCA2)) + geom_point(aes(color = group)) +
    geom_path(data=df_ell, aes(x=PCA1, y=PCA2,colour=group), size=1, linetype=2)+
    annotate("text",x=PCA.mean$PCA1,y=PCA.mean$PCA2,label=PCA.mean$group)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()


# NMF clustering
data <- read.csv("EXP_matrix.csv", row.names = 1, stringsAsFactors = F)
data <- as.matrix(data)
data <- FSbyMAD(data, cut.type="topk",value=1500)

if (all(data >= 0 ))
   print("checked") else
   (data <- (data*2)) & (data[data < 0] <- 0)

estim.r <- nmf(data, 2:6, nrun=100, seed=123456)

pdf(file="estim.pdf",width=4.5,height=5)
plot(estim.r)
dev.off()

res <- nmf(data, 2, nrun=100, seed=123456)
consensusmap(res, labCol=NA, labRow=NA)

pdf(file="consensusmap.pdf",width=4.5,height=5)
consensusmap(estim.r, labCol=NA, labRow=NA)
dev.off()

classifier <- as.matrix(apply(coef(res),2,which.max))
row.names(classifier) <- gsub(".", "-", row.names(classifier), fixed=TRUE) 


# prepare input files for subclass mapping
generateInputFileForSubMap <- function(in_gct, gct_file, cls_file, sam_info, type_name = "type"){
  in_gct <- data.frame(GeneID=rownames(in_gct),
                       description="na",
                       in_gct, 
                       stringsAsFactors = F,
                       check.names = F)
  cat("#1.2\n", file = gct_file)
  cat(nrow(in_gct),"\t",ncol(in_gct)-2,"\n", file = gct_file, append = T)
  cat(paste(colnames(in_gct), collapse = "\t"),"\n", file = gct_file, append = T)
  for(i in 1:nrow(in_gct)) cat(paste(in_gct[i,], collapse = "\t"),"\n", file = gct_file, append = T)
  
  cat(nrow(sam_info),length(levels(factor(sam_info$rank))),1, "\n", file = cls_file )
  cat("#", paste0(levels(factor(sam_info[, type_name])), collapse = " " ), "\n", file = cls_file, sep = "", append = T)
  cat(as.numeric(factor(sam_info[, type_name])), file = cls_file, append = T)
}

inputmatrix1 <- read.table("inputmatrix1.txt",sep = "\t",row.names = 1,header = T,stringsAsFactors = F) 
info1 <- read.table("cluster_info1.txt",sep = "\t",row.names = 1,header = T,stringsAsFactors = F)
info1 <- info1[order(info1$type),]
info1$rank <- rep(c(1,2),times=as.character(table(info1$type)))

inputmatrix2 <- read.table("inputmatrix2.txt",sep = "\t",stringsAsFactors = F,header = T,row.names = 1)
info2 <- read.table("cluster_info2.txt",sep = "\t",row.names = 1,header = T,stringsAsFactors = F)
info2 <- info2[order(info2$type),]
info2$rank <- rep(c(1,2),times=as.character(table(info2$type)))

GENELIST <- intersect(rownames(inputmatrix2),rownames(inputmatrix1))

sam_info <- info1
in_gct <- log2(inputmatrix1[GENELIST,rownames(sam_info)]+1)
gct_file <- "inputmatrix1.gct"
cls_file <- "inputmatrix1.cls"
generateInputFileForSubMap(in_gct = in_gct, gct_file = gct_file, cls_file = cls_file, sam_info = sam_info, type_name = "rank")

sam_info <- info2
in_gct <-inputmatrix2[GENELIST,rownames(sam_info)]
gct_file <- "inputmatrix2.gct"
cls_file <- "inputmatrix2.cls"
generateInputFileForSubMap(in_gct = in_gct, gct_file = gct_file, cls_file = cls_file, sam_info = sam_info, type_name = "rank")


# Gene set enrichment analysis
gsym.fc <- read.table("GSEA_input.txt", header = T)
geneset1 <- read.gmt("h.all.v7.0.entrez.gmt")
geneset2 <- read.gmt("c2.cp.kegg.v7.1.entrez.gmt")
geneset3 <- read.gmt("c2.cp.reactome.v7.1.entrez.gmt")
geneset4 <- read.gmt("c5.mf.v7.1.entrez.gmt")

gsym.id <- bitr(gsym.fc$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
gsym.fc.id <- merge(gsym.fc, gsym.id, by="SYMBOL", all=F)
gsym.fc.id.sorted <- gsym.fc.id[order(gsym.fc.id$logFC, decreasing = T),]

id.fc <- gsym.fc.id.sorted$logFC
names(id.fc) <- gsym.fc.id.sorted$ENTREZID

kk1 <- GSEA(id.fc,TERM2GENE =geneset1,pvalueCutoff = 0.05)
kk2 <- GSEA(id.fc,TERM2GENE =geneset2,pvalueCutoff = 0.05)
kk3 <- GSEA(id.fc,TERM2GENE =geneset3,pvalueCutoff = 0.05)
kk4 <- GSEA(id.fc,TERM2GENE =geneset4,pvalueCutoff = 0.05)

kk1.gsym <- setReadable(kk1, 'org.Hs.eg.db', 'ENTREZID')
sortkk1 <- kk.gsym1[order(kk.gsym1$enrichmentScore, decreasing = T),]

kk2.gsym <- setReadable(kk2, 'org.Hs.eg.db', 'ENTREZID')
sortkk2 <- kk.gsym2[order(kk.gsym2$enrichmentScore, decreasing = T),]

kk3.gsym <- setReadable(kk3, 'org.Hs.eg.db', 'ENTREZID')
sortkk3 <- kk.gsym3[order(kk.gsym3$enrichmentScore, decreasing = T),]

kk4.gsym <- setReadable(kk4, 'org.Hs.eg.db', 'ENTREZID')
sortkk4 <- kk.gsym4[order(kk.gsym4$enrichmentScore, decreasing = T),]


# get purified data
normal <- read.csv("Normal_matrix.csv",row.names = 1, stringsAsFactors = F) 
tumor <- read.csv("Tumor_matrix.csv",row.names = 1, stringsAsFactors = F) 
normal <- as.matrix(normal+1)
tumor <- as.matrix(tumor+1)

ISOpureS1model <- ISOpure.step1.CPE(tumor, normal)
ISOpureS2model <- ISOpure.step2.PPE(tumor, normal, ISOpureS1model)

# estimate response
trainExpr  <- read.csv("EXP_matrix.csv",row.names = 1, stringsAsFactors = F)
trainPtype  <-read.csv("response.csv",row.names = 1, stringsAsFactors = F)
testExpr  <- read.csv("EXP_matrix_test.csv",row.names = 1, stringsAsFactors = F)
trainExpr <- t(scale(t(as.matrix(trainExpr))))
testExpr <- t(scale(t(as.matrix(testExpr))))
outTab=data.frame()
for (drug in rownames(trainPtype) ){
trainPtype_1 <- as.numeric(trainPtype[drug,])
ptypeOut <- calcPhenotype(trainExpr, trainPtype_1, testExpr, selection=1)
outTab<-rbind(outTab,ptypeOut)
}
outTab <- as.matrix(outTab)
row.names(outTab) <- row.names(trainPtype)
colnames(outTab) <- colnames(testExpr)

# MCP-counter
data <- read.csv("EXP_matrix.csv",row.names = 1, stringsAsFactors = F)
estimateCell <- MCPcounter.estimate(data, featuresType="HUGO_symbols")
out=rbind(ID=colnames(data),estimateCell)



