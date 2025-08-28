#gsea 做kegg   GO_BP（生物过程）和GO_MF（分子功能）
BiocManager::install("org.Sc.sgd.db")

#library
suppressMessages({ 
    library(dplyr)
    library(org.Sc.sgd.db)
    library(clusterProfiler)
    library(DOSE)
    library(ggplot2)
    })

#function
save_plot <- function(od, filename, plot, width = 6, height = 6, scale = 1, dpi = 300){
  png <- paste0(filename, ".png")
  pdf <- paste0(filename, ".pdf")
  tryCatch({#ggplot object
    ggsave(filename = file.path(od, png), plot = plot, width = width, height = height, scale = scale, dpi = dpi)
    ggsave(filename = file.path(od, pdf), plot = plot, width = width, height = height, scale = scale, dpi = dpi)},
    error = function(e) {tryCatch({
      png(file.path(od, png), width = width, height = height, units = "in", res = dpi, type = "cairo")
      print(plot)
      dev.off()

      pdf(file.path(od, pdf), width = width, height = height)
      print(plot)
      dev.off()
    },
    error = function(e){
      print(paste0("sample: ", filename, " save plot error"))
    }
    )
    }
  )
}
GSEA_result_out <- function(od,result,out_prefix,topn = 10){
    if(nrow(result) > 0 ){
    write.table(as.data.frame(result),file = file.path(od,paste0(out_prefix,".GSEA.txt")),row.names = F, quote = F, sep = "\t")
    #dotplot
    (dotplot(result, showCategory=topn, split=".sign") + facet_grid(.~.sign)) %>%
        save_plot(od = od,filename = paste0(out_prefix,".dotplot.top",topn),width = 12,plot = .)
    #GSEA plot
    if(nrow(result) < topn){
        #pathway_draw <- result$pathway
        for(x in 1:nrow(result)){
        (gseaplot(result, by = "all", title = result$Description[x], geneSetID = x)) %>%
            save_plot(od = od,filename = paste0(out_prefix,".GSEA.top",x),plot = .)
        }
    }else{
        #pathway_draw <- result$pathway[1:topn]
        for(x in 1:topn){
        (gseaplot(result, by = "all", title = result$Description[x], geneSetID = x)) %>%
            save_plot(od = od,filename = paste0(out_prefix,".GSEA.top",x),plot = .)
        }
    }
    }
    else{
    print("No GSEA result")
    }
}

#标准富集分析

expr = read.table("Yeast_all.txt",header=T,sep="\t")
bitr_gene <- bitr(geneID=geneID,fromType="ENSEMBL",toType="ENTREZID",OrgDb=org.Sc.sgd.db)
expr %<>% left_join(.,bitr_gene,by=c("id"="ENSEMBL"))
gene_list <-  dplyr::arrange(expr,desc(log2FoldChange)) %>%
  .[,"log2FoldChange"]
names(gene_list) <- expr %>%
  dplyr::arrange(.,desc(log2FoldChange)) %>%
  .[,"ENTREZID"]

#
go_BP <- gseGO(
      geneList  = gene_list,
      OrgDb  = "org.Sc.sgd.db",
      ont  = "BP",
      keyType = "ENTREZID",
      minGSSize  = 10,
      maxGSSize  = 500,
      pvalueCutoff = 0.05,
      verbose  = F)

GSEA_result_out(od="./",result=go_BP,out_prefix="Yeast_BP")

go_MF <- gseGO(
      geneList  = gene_list,
      OrgDb  = "org.Sc.sgd.db",
      ont  = "MF",
      keyType = "ENTREZID",
      minGSSize  = 10,
      maxGSSize  = 500,
      pvalueCutoff = 0.05,
      verbose  = F)

GSEA_result_out(od="./",result=go_MF,out_prefix="Yeast_MF")

#--> Expected input gene ID: YOR136W,YJL121C,YMR083W,YMR110C,YGL062W,YJL097W
names(gene_list) <- expr %>%
  dplyr::arrange(.,desc(log2FoldChange)) %>%
  .[,"id"]

kegg <- gseKEGG(
  geneList  = gene_list,
  keyType  = 'kegg',
  organism = 'sce',
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  verbose = F
)

GSEA_result_out(od="./",result=kegg,out_prefix="Yeast_KEGG")