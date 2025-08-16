library(limma)
library(Glimma)
library(org.Mm.eg.db)
library(RColorBrewer)
library(gplots)
library(NMF)
library(DESeq2)
#Reading the count matrix that stringtie generated
#countData <- read.delim("count.txt", header=T, row.names ="Geneid", check.names = FALSE)
countData <- as.matrix(read.csv("gene_count_matrix.csv", row.names="gene_id"))

gene_IDS = as.list(rownames(countData))

pheno_data = read.csv("Sample_name.csv")
pheno_data
pheno_data <- pheno_data[-(1:12),]
pheno_data <- pheno_data[,-(4)]
head(pheno_data)
head(countData)
countdata <- countData[,-(1:12)]
head(countdata)
dim(countdata)
head(pheno_data)
pheno_data
ann <- read.delim("Ent11IV18_1_GeneCatalog_20221108_edited.txt",header=T, row.names ="geneID", check.names = FALSE)
head(ann)
# The count matrix had all the conditions so I had delete the ones to only have conditions which is 0micro molar and 10micro molar
colnames(countdata)
x <- DGEList(countdata)
names(x)
x$samples
group <-paste(pheno_data$Al.concentration, pheno_data$Replicate.number, sep=".")
group
x$samples
myCPM <- cpm(countdata)
head(myCPM)
# Made threshold such as cpm> 0.5 (10 counts)
thresh <- myCPM >0.5
head(thresh)
table(rowSums(thresh))
keep <- rowSums(thresh) >=2
summary(keep)
x <- x[keep, keep.lib.sizes=FALSE]
#x$samples$lib.size
dim(x)
plot(myCPM[,1],countdata[,1],ylim=c(0,50),xlim=c(0,3))
x$samples$lib.size
barplot(x$samples$lib.size, names=colnames(x), las=2)
logcounts <- cpm(x,log = TRUE)
boxplot(logcounts,xlab="",ylab = "Log2 counts per million", las=2)
# Normalizing reads based on library size
x <- calcNormFactors(x)
x$samples
logcounts <- cpm(x,log = TRUE)
head(x)
plotMDS(x, gene.selection="common")
cor(logcounts, method="pearson")
heatmap(as.matrix(cor(logcounts, method="pearson")), 
        main="Clustering of Pearson correlations", scale="none")
cor(logcounts, method="spearman")
heatmap(as.matrix(cor(logcounts, method="spearman")), 
        main="Clustering of Spearman correlations", scale="none", col = heat.colors(30), margins=c(5,5))
legend("bottomright", legend=c("Low", "Medium", "High"), fill=heat.colors(3), bty="n")
plotMD(logcounts,column = 6)
abline(h=0,col="grey")
group <-paste(pheno_data$Al.concentration)
# Edge R and Limma
design <- model.matrix(~ 0 + group)
design
group

v <- voom(x, design, plot=TRUE)
v
names(v)
boxplot(v$E, xlab="", ylab="Log2 counts per million", las=2, main="Voom transformed logCPM")
abline(h=median(v$E),col="blue")
fit <- lmFit(v)
names(fit)
design

cont.matrix <- makeContrasts(extremes=group0 - group10, levels = design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
head(fit.cont)
view(fit.cont)
dim(fit.cont)
summa.fit <- decideTests(fit.cont, p.value = 0.05, lfc=0.5)
summa.fit[,1]
summa.fit <- decideTests(fit.cont)
summary(summa.fit)

fit.cont$gene_IDS = rownames(fit.cont)
fit.cont

library(tibble)


res <- rownames_to_column(as.data.frame(fit.cont),"geneIDS")
res <- res[,(-1)]
res
colnames(res)[13] <- "geneID"
res
ann <- rownames_to_column(as.data.frame(ann),"geneID")
res_final <- merge(res,ann, "geneID")
view(res_final)

#write.csv(res_final, 'controlextreme_featurecounts_proteinID.csv')
head(fit.cont)
volcanoplot(fit.cont, coef=1, highlight = 100, names=fit.cont$geneID$SYMBOL, main="10 micro molar vs 0 micro molar")
group2 <- group
levels(group2) <- c("0","10")



glXYPlot(x=fit.cont$coefficients[,1], y=fit.cont$lods[,1],
         xlab="logFC", ylab="B", main="B.PregVsLac",
         counts=v$E, groups=group2, status=summa.fit[,1],
         ,folder="volcano")
glXYPlot(x=fit.cont$lods[,1], y=-log(fit.cont$F.p.value[,1]),
         xlab="logFC", ylab="B", main="extreme",
         counts=v$E, groups=group2, status=summa.fit[,1],
         anno=fit.cont$genes, side.main="ENTREZID", folder="volcano")


library("DESeq2")
countDataD <- as.matrix(read.csv("gene_count_matrix.csv", row.names="gene_id"))
countData <- as.matrix("ramyafeatcounts_genelevel.txt", row.names ="Geneid", check.names = FALSE)
colData <- read.csv("Sample_name.csv", row.names=1)
colData <- colData[-(1:12),]

### Ignore DESeq part wasnt working
countDataD <- countDataD[,-(1:12)]
all(rownames(colData) %in% colnames(countDataD))
countDataD <- countDataD[, rownames(colData)]
all(rownames(colData) == colnames(countDataD))
head(countDataD)
dds <- DESeqDataSetFromMatrix(countData = countDataD,
                              colData = colData, design  = ~ Al.concentration)
dds <- DESeq(dds)
fit.cont
res <- results(dds, alpha = 0.05, cooksCutoff = FALSE)
plotMA(res, ylim=c(-2,2))           
plotMA(resLFC, ylim=c(-2,2))

