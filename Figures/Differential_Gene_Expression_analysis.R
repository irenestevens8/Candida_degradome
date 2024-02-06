ctx<-as.matrix(read.delim("/Users/vip/Documents/Candida/RNA-seq/candida_RNAseq_CDS_counts.txt", sep="\t", row.names = 1, header=T, check.names=F,stringsAsFactors = TRUE))
coldata<-read.delim("/Users/vip/Documents/Candida/RNA-seq/metadata.txt",row.names=1,stringsAsFactors = TRUE)
coldata$Treatment<-factor(coldata$Treatment)
coldata$Strain<-factor(coldata$Strain)
coldata$Treatment=relevel(coldata$Treatment, 'Control')
coldata$Strain=relevel(coldata$Strain, 'sc5314')

dds <- DESeqDataSetFromMatrix(countData=ctx, colData = coldata, design <- ~Strain + Treatment + Strain:Treatment)
dds<-DESeq(dds)
#sd <- rlog(dds,blind=TRUE)          
sd <- vst(dds,blind=TRUE)
plotPCA(sd,intgroup=c("Treatment","Strain"))
plot_data <- plotPCA(sd,intgroup=c("Treatment","Strain"),returnData=TRUE)

#read in counts and metadata
ctx<-as.matrix(read.delim("/Users/vip/Documents/Candida/RNA-seq/counts_sc5314.txt", sep="\t", row.names = 1, header=T, check.names=F,stringsAsFactors = TRUE))
coldata<-read.delim("/Users/vip/Documents/Candida/RNA-seq/metadata_sc5314.txt",row.names=1,stringsAsFactors = TRUE)

coldata<-coldata[,c("Treatment","Replicate")]
coldata$Treatment<-factor(coldata$Treatment)
coldata$Treatment=relevel(coldata$Treatment, 'Control')


dds <- DESeqDataSetFromMatrix(countData=ctx,  colData = coldata, design <- ~Treatment)
dds<-DESeq(dds)
sd <- rlog(dds,blind=TRUE)


plotPCA(sd,intgroup=c("Treatment"))

res = results(dds, contrast=c("Treatment","Fluconazole","Control"))
ix = which.min(res$padj)
res <- res[order(res$padj),] 	
summary(res)

write.csv(as.data.frame(res),file='/Users/vip/Documents/Candida/RNA-seq/DEG_sc5314.txt')

ctx<-read.delim("/Users/vip/Documents/Candida/RNA-seq/DEG_sc5314_geneid.txt", sep="\t", header=T)
vol_plot <- ctx %>% ggplot(aes(x = log2FoldChange, y = -log10(padj))) + geom_point()
vol_plot

#gene names
ctx<-read.delim("/Users/vip/Documents/Candida/RNA-seq/DEG_sc5314.txt", sep="\t", header=T)
anno<-read.delim("/Users/vip/Documents/CUG_reprocessed/Figure 2/annotation_pretty.txt", sep="\t", header=T, check.names=F)
common_col<-'cds_id'
merged.sc<-merge(ctx,anno,by=common_col)
write.table(merged.sc, file="/Users/vip/Documents/Candida/RNA-seq/DEG_sc5314_geneid.txt", sep="\t", quote=F, row.names=T, col.names=NA)
