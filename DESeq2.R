#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("DESeq2")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("BiocParallel")

#install.packages("pheatmap")
#install.packages("ggplot2")


#引用包
library(limma)
library(DESeq2)
library(pheatmap)
library(ggplot2)
library("BiocParallel")
register(MulticoreParam(5))

expFile="expMatrix.txt"      #表达输入文件
logFCfilter=1             #logFC过滤条件
pvalueFilter=0.05           #矫正后的p值过滤条件
conFile="s1.txt"          #对照组样品
treatFile="s2.txt"        #实验组样品
setwd("C:\\Users\\lexb\\Desktop\\yzng2dd\\02.DESeq2")      #设置工作目录

#读取输入文件，并对输入文件进行整理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]
data=round(data,0)

#读取样品信息
sample1=read.table(conFile,sep="\t",header=F,check.names=F)
sample2=read.table(treatFile,sep="\t",header=F,check.names=F)
conData=data[,as.vector(sample1[,1])]
treatData=data[,as.vector(sample2[,1])]
data=cbind(conData,treatData)
conNum=ncol(conData)
treatNum=ncol(treatData)

#差异分析
group=c(rep("normal",conNum), rep("tumor",treatNum))
coldata=data.frame(condition=group, row.names=colnames(data))
dds=DESeqDataSetFromMatrix(countData=data, colData=coldata, design=~condition)
dds=DESeq(dds, parallel=TRUE)
diff=results(dds, parallel=TRUE)
diff=as.data.frame(diff)

#输出差异结果
diff=diff[is.na(diff$pvalue)==FALSE,]
diff=diff[order(diff$pvalue),]

Sig=ifelse(diff$log2FoldChange>logFCfilter & diff$pvalue<pvalueFilter,"Up",ifelse(diff$log2FoldChange< -logFCfilter & diff$pvalue<pvalueFilter,"Down","No-sig"))
diffOut=cbind(diff, Sig)
diffOut=rbind(id=colnames(diffOut), diffOut)
write.table(diffOut, file="all.txt", sep="\t", quote=F, col.names=F)
diffSig=diff[(diff$pvalue<pvalueFilter & (diff$log2FoldChange>logFCfilter | diff$log2FoldChange<(-logFCfilter))),]

Sig=ifelse(diffSig$log2FoldChange>logFCfilter & diffSig$pvalue<pvalueFilter,"Up",ifelse(diffSig$log2FoldChange< -logFCfilter & diffSig$pvalue<pvalueFilter,"Down","No-sig"))
diffSigOut=cbind(diffSig, Sig)
diffSigOut=rbind(id=colnames(diffSigOut),diffSigOut)
write.table(diffSigOut, file="diff.xls",sep="\t",quote=F,col.names=F)
write.table(diffSigOut, file="diff.txt",sep="\t",quote=F,col.names=F)

#输出矫正后的表达量
newData=counts(dds, normalized=TRUE)
normalizeExp=rbind(id=colnames(newData), newData)
write.table(normalizeExp,file="normalize.txt",sep="\t",quote=F,col.names=F)
diffExpOut=rbind(id=colnames(newData), newData[rownames(diffSig),])
write.table(diffExpOut,file="diffExp.txt",sep="\t",quote=F,col.names=F)

#绘制差异基因热图
geneNum=100
diffSig=diffSig[order(as.numeric(as.vector(diffSig$log2FoldChange))),]
diffGeneName=as.vector(rownames(diffSig))
diffLength=length(diffGeneName)
hmGene=c()
if(diffLength>(2*geneNum)){
    hmGene=diffGeneName[c(1:geneNum,(diffLength-geneNum+1):diffLength)]
}else{
    hmGene=diffGeneName
}
hmExp=newData[hmGene,]
hmExp=log2(hmExp+1)

Type=c(rep("Control",conNum),rep("Treat",treatNum))
names(Type)=colnames(newData)
Type=as.data.frame(Type)
pdf(file="heatmap.pdf",height=7,width=8)
pheatmap(hmExp, 
         annotation=Type, 
         color = colorRampPalette(c("blue", "white", "red"))(50),
         cluster_cols =F,
         show_colnames = F,
         scale="row",
         fontsize = 8,
         fontsize_row=3,
         fontsize_col=8)
dev.off()

#定义显著性
Significant=ifelse((diff$pvalue<pvalueFilter & abs(diff$log2FoldChange)>logFCfilter), ifelse(diff$log2FoldChange>logFCfilter,"Up","Down"), "Not")
#绘制火山图
p = ggplot(diff, aes(log2FoldChange, -log10(pvalue)))+
    geom_point(aes(col=Significant))+xlim(-6,6)+
    scale_color_manual(values=c("green", "grey", "red"))+
    labs(title = " ")+
    theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))
p=p+theme_bw()
#保存为图片
pdf("vol.pdf",width=5.5,height=5)
print(p)
dev.off()




