#
library(DESeq2)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(M3C)

library(msigdf)
library(edgeR)
library(corrplot)

# OE
#library(matrixStats)
#library(survival)
#library(ROCR)
#library(Hmisc)
#library(lme4)
#library(lmerTest);attach(mtcars)

#library(scde)
#library(ppcor)
#library(rms)
#library(mixtools)
#library(plotrix)

#library(Seurat)



## some from https://github.com/SimonLab-RU/DESeq2/blob/master/DESeq2.Rmd

plotPCA_DR <- function (object, intgroup = "condition", ntop = 500, PCs = c(1,2),
                        size = 7, alpha = 1, colors = NULL, returnData = FALSE) 
{
  PC_x <- paste("PC", PCs[1], sep = "")
  PC_y <- paste("PC", PCs[2], sep = "")
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
  percentVar <- pca$sdev^2 / sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
  
  #return(pca$x)
  #return(pca$x[, c(PCs[1], PCs[2])])
  
  # To use different shapes and fill using different colors
  d <- data.frame(PCx = pca$x[, PCs[1]], PCy = pca$x[, PCs[2]],
                  to_fill = factor(intgroup.df[,1]), to_shape = factor(intgroup.df[,2]),
                  name = colnames(object))
  
  if (returnData) {
    attr(d, "percentVar") <- percentVar[PCs[1]:PCs[2]]
    return(d)
  }
  
  ggplot(data = d, aes_string(x = "PCx", y = "PCy", fill = "to_fill", shape = "to_shape")) +
    labs(fill = intgroup[1], shape = intgroup[2]) + geom_point(size = size, alpha = alpha) +
    scale_fill_manual(values = nice_colors, guide = guide_legend(override.aes = aes(shape = 21))) +
    scale_shape_manual(values = c(24, 21, 22, 25)) + coord_fixed() + theme_bw() +
    xlab(paste0(PC_x, ": ", round(percentVar[PCs[1]] * 100), "% variance")) +
    ylab(paste0(PC_y, ": ", round(percentVar[PCs[2]] * 100), "% variance"))
}


## OE
# https://github.com/livnatje/ImmuneResistance/blob/master/Code/

discretize<-function(v,n.cat){
  q1<-quantile(v,seq(from = (1/n.cat),to = 1,by = (1/n.cat)))
  u<-matrix(nrow = length(v))
  for(i in 2:n.cat){
    u[(v>=q1[i-1])&(v<q1[i])]<-i
  }
  return(u)
}

get.semi.random.OE <- function(r,genes.dist.q,b.sign,num.rounds = 1000,full.flag = F){
  # Previous name: get.random.sig.scores
  
  # sign.q : count signature genes located in some bins
  sign.q<-as.matrix(table(genes.dist.q[b.sign]))
  # q : located bins
  q<-rownames(sign.q)
  idx.all<-c()
  B<-matrix(data = F,nrow = length(genes.dist.q),ncol = num.rounds)
  Q<-matrix(data = 0,nrow = length(genes.dist.q),ncol = num.rounds)   # Q has nothing to do here
  
  # B each col is an index for same number of genes randomly selected in same bins 
  for (i in 1:nrow(sign.q)){
    num.genes<-sign.q[i]
    if(num.genes>0){
      # index of all genes in that bin (q[i])
      idx<-which(is.element(genes.dist.q,q[i]))
      for (j in 1:num.rounds){
        idxj<-sample(idx,num.genes) 
        Q[i,j]<-sum(B[idxj,j]==T)    # stupid Q, always zero matrix, waste of time to doubt it
        B[idxj,j]<-T
      }  
    }
  }
  rand.scores<-apply(B,2,function(x) colMeans(r$zscores[x,]))   # get mean of 'zscore's of one round
  if(full.flag){return(rand.scores)}
  rand.scores<-rowMeans(rand.scores)  # get mean of num.rounds rounds
  return(rand.scores)
}

## inlocal
zscore_mat <- function(mat){
  Mean <- rowMeans(mat)
  Sd <- apply(mat,1,sd)
  mat <- sweep(mat,1,Mean,FUN='-')
  mat <- sweep(mat,1,Sd,FUN='/')
  return(mat)
}

##
##

# write a function to quickly do this   

# input   
# mat_e : expression matrix(CPM/TPM)  
# cells_s : cells_selected(character vector)  
# path_n : pathway_names(character list,path_o="pathwah way")    
# gene_sign : list(path_o=path_g)  
#    (path_o : pathway_names in short)   
#    (path_g : genes in this pathway, character vector)   
# seed_r : random seed  

# output list:  
#   list$stat : table of pathway/expressed genes  
# list$OE : OE of sorted cells  
# list$mat_z : Zscore of sorted cells/genes
# list$bar : bar plot  
# list$heat : heat map  

# mod and debug: 
#   unit var names and align dimensions  
# pheat(modified by UncleY with additional class 'pheatmap' and par 'silent=T')  
# but still pheat object can't plot in rmd, just use 'silent=F' ~  


# when function is done, copy into analysis.r  

easy_OE <- function(mat_e,cells_s,path_n,gene_sign,seed_r=7788){
  
  ret <- list()
  ret$tpm <- log2(mat_e[,cells_s]+1)
  ret$tpm <- ret$tpm[rowSums(ret$tpm)>0,]
  ret$genes <- rownames(ret$tpm)
  
  #
  set.seed(seed_r)
  
  ret$genes.mean <- rowMeans(ret$tpm)
  ret$genes.sd <- apply(ret$tpm,1,sd)
  ret$zscores <- sweep(ret$tpm,1,ret$genes.mean,FUN='-')
  ret$zscores <- sweep(ret$zscores,1,ret$genes.sd,FUN='/')
  
  ret$genes.dist <- ret$genes.mean
  ret$genes.dist.q <- discretize(ret$genes.dist, n.cat=50)
  ret$sig.scores <- matrix(data=0,nrow=ncol(ret$tpm),ncol=length(gene_sign))
  
  ret$sig.names <- names(gene_sign)   # path_o
  colnames(ret$sig.scores) <- ret$sig.names
  rownames(ret$sig.scores) <- colnames(ret$tpm)
  
  ret$sig.scores.raw <- ret$sig.scores
  ret$sig.rand.scores <- ret$sig.scores
  
  ret$mat_z <- list()
  ret$heat <- list()
  ret$bar <- list()
  
  ret$stat <- list()
  
  for(i in ret$sig.names){
    b.sign <- is.element(ret$genes, gene_sign[[i]])
    
    # scores
    ret$sig.rand.scores[,i] <- get.semi.random.OE(ret,ret$genes.dist.q,b.sign,num.rounds=100)
    ret$sig.scores.raw[,i] <- colMeans(ret$zscores[b.sign,])
    ret$sig.scores[,i] <- ret$sig.scores.raw[,i]-ret$sig.rand.scores[,i]
    ret$sig.scores[,i] <- round(ret$sig.scores[,i],3)
    # ret$sig.scores[,i] <- sort(ret$sig.scores[,i],decreasing=TRUE)
    # here can't sort, could only sort numbers but no names sorted, sort in OE barplot
    new_order <- order(ret$sig.scores[,i],decreasing = T)
    
    # OE barplot    
    ret$bar[[i]] <- ggplot(data=cbind.data.frame(Score=(ret$sig.scores[,i])[new_order],
                                                 Name=factor(names(ret$sig.scores[,i])[new_order],levels=(names(ret$sig.scores[,i]))[new_order])),
                           mapping=aes(x=Score,y=Name)) +
      geom_bar(stat='identity') +
      #coord_flip() +
      labs(y="",x=paste0("Overall Expression of geneset:\n",path_n[[i]]))
    
    # mat_z
    ret$mat_z[[i]] <- zscore_mat(ret$zscores[b.sign,])
    
    # sort genes by mean value distance: mean(OE>0) - mean(OE<0) 
    idx_cells.up <- names(ret$sig.scores[,i][ret$sig.scores[,i]>0])
    idx_cells.down <- names(ret$sig.scores[,i][ret$sig.scores[,i]<0])
    
    idx_genes <- rowSums(ret$mat_z[[i]][,idx_cells.up])-rowSums(ret$mat_z[[i]][,idx_cells.down])
    idx_genes <- sort(idx_genes,decreasing=TRUE)
    
    ret$mat_z[[i]] <- ret$mat_z[[i]][names(idx_genes),rev((names(ret$sig.scores[,i]))[new_order])]
    
    
    # mat_z heatmap
    ret$heat[[i]] <- pheatmap::pheatmap(t(t(ret$mat_z[[i]])), cluster_cols=FALSE,cluster_rows=FALSE,fontsiize_row=7.5,
                                        main=paste0("Zscore of genes in geneset: ",path_n[[i]]),
                                        color=colorRampPalette(c("blue","white","red"))(100),
                                        breaks=seq(-2,2,0.04))
    
    # stat 
    ret$stat[[i]] <- rbind("*** Stat Table ***",
                           paste0("Pathway: ",path_n[[i]]),
                           paste0("total genes: ",length(gene_sign[[i]])),
                           paste0("expressed genes: ",sum(b.sign)),
                           #paste(ret$genes[b.sign],collapse=" ")
                           paste(rownames(ret$mat_z[[i]]),collapse=" ")
    )
  }
  
  
  # output
  rett <- list()
  
  rett$stat <- ret$stat
  rett$OE <- ret$sig.scores
  rett$mat_z <- ret$mat_z
  rett$bar <- ret$bar
  rett$heat <- ret$heat
  
  return(rett)
}

#### additional functions to add  

# filter matrix for minimum expression and minimum expressed cells  
filt_mat <- function(mat,min_expr=2,min_cells=3){
  mat <- mat[rowSums(mat > min_expr) >= min_cells,]
  return(mat)
}



#define a function to plot correlation  
plotCor <- function(Mat){
  library(edgeR)
  library(corrplot)
  Cor <- cor(log2(edgeR::cpm(Mat)+1))
  par(cex=0.8, pin=c(8,8))
  corrplot(Cor,method="number",title = "pearson correlation of log2(CPM+1)",mar = c(0, 0, 1, 0))
}

# function for easy volcano plot and heatmap on DEGs
finalplot <- function(Dat, Ret, Cnt, Pval, FC, Font=8, Sign=FALSE, Sign_dn = 10, Sign_up = 10, Label=NULL, padjust = TRUE, heatmin=0){
  
  library(ggpubr)
  library(ggthemes)
  library(ComplexHeatmap)
  
  if(padjust==FALSE){
    Ret$padj <- Ret$pvalue
    volcanoy <- "- log10 pvalue"
  }else{
    volcanoy <- "- log10 p.adjust"
  }
  
  Ret$pp <- -log10(Ret$padj)
  Ret$Group <- "not_significant"
  Ret$Group[which((Ret$padj < Pval) & (Ret$log2FoldChange > log2(FC)))] <- "up-regulated"
  Ret$Group[which((Ret$padj < Pval) & (Ret$log2FoldChange < -log2(FC)))] <- "down-regulated"
  
  Ret$Genes <- ""
  Genes_up <- rownames(Ret)[which((Ret$padj < Pval) & (Ret$log2FoldChange > log2(FC)))]
  Genes_down <- rownames(Ret)[which((Ret$padj < Pval) & (Ret$log2FoldChange < -log2(FC)))]
  
  if(Sign==TRUE){
    Ret$Genes[match(c(Genes_up,Genes_down),rownames(Ret))] <- c(Genes_up,Genes_down)
    Marker <- c(Genes_up,Genes_down)
  }else{
    #Genes_up_top <- rownames(Ret[Ret$log2FoldChange>0,])[head(order(Ret[Ret$log2FoldChange>0,]$padj,decreasing = F),Sign_up)]
    #Genes_down_top <- rownames(Ret[Ret$log2FoldChange<0,])[head(order(Ret[Ret$log2FoldChange<0,]$padj,decreasing = F),Sign_dn)]
    #Genes_up_top <- rownames(Ret %>% filter(log2FoldChange > log2(FC)) %>% arrange(padj))[1:Sign_up]
    #Genes_down_top <- rownames(Ret %>% filter(log2FoldChange < -log2(FC)) %>% arrange(padj))[1:Sign_dn]
    Genes_up_top <- Genes_up[1:Sign_up]
    Genes_down_top <- Genes_down[1:Sign_dn]
    
    Marker <- c(Genes_up_top,Genes_down_top)
    for(M in Marker){
      Ret$Genes[rownames(Ret)==M] <- M
    }
  }
  
  # for custom labeled genes, 20210323 added
  if(!is.null(Label)){
    Ret$Genes <- ""
    for(M in Label){
      Ret$Genes[rownames(Ret)==M] <- M
    }
  }
  
  
  
  volplot <- ggscatter(Ret, x = "log2FoldChange", y = "pp", color = "Group",
                       palette = c("#2f5688","#BBBBBB","#CC0000"),size = 1.5,
                       #label = "Genes",
                       font.label = Font, repel = F,
                       xlab = "log2 fold change", ylab=volcanoy) + theme_classic() +
    geom_text_repel(mapping = aes(x = Ret$log2FoldChange, y = Ret$pp, label = Ret$Genes), max.overlaps = 500) +
    guides(color=guide_legend(reverse=TRUE,title = Cnt)) +
    geom_vline(xintercept = c(-log2(FC),log2(FC)), linetype="dashed") +
    geom_hline(yintercept = -log10(Pval), linetype="dashed") 
  #volplot
  
  Ret_mat <- Dat[c(Genes_up,Genes_down),]
  Ret_cor <- cor(Ret_mat)
  
  lower = heatmin
  upper = 1
  pal <- "Reds"
  
  ht1 <- Heatmap(Ret_cor, col = circlize::colorRamp2(seq(lower, upper, ((upper-lower)/7)),RColorBrewer::brewer.pal(8, pal)),
                 heatmap_legend_param = list(
                   color_bar = "continuous",
                   legend_direction = "horizontal",
                   legend_width = unit(5, "cm"),
                   title_position = "topcenter"),
                 name = "Pearson correlation",
                 column_names_gp = gpar(fontsize = 10),
                 row_names_gp = gpar(fontsize = 10),
                 top_annotation = NULL)
  heatplot <- draw(ht1, heatmap_legend_side = "top")
  
  
  Genes_len <- cat("P: ",Pval,", ","FC: ",FC,"\n\n",
                   "up: ",length(Genes_up),"\n","down: ",length(Genes_down),"\n","total: ",length(c(Genes_up,Genes_down)),
                   "\n",sep="")
  
  
  
  final_sig <- list(up= Genes_up, down=Genes_down, len=Genes_len, vol=volplot, heat=heatplot, marker=Marker)
  return(final_sig)  
}

#
# rename edger result to match deseq2 function
edg2des <- function(retr){
  colnames(retr) <- c("log2FoldChange","logCPM","LR","pvalue","padj")
  return(retr)
}


# write a function to recover gene column of DEGseq2 result
# and add new FC column to recover log2FC with "+/-"
rec_gcol <- function(res){
  res <- cbind(rownames(res),res)
  colnames(res)[1] <- "gene"
  
  res$FC <- res$log2FoldChange
  res$FC[res$FC<0] <- -(2^(-res$FC[res$FC<0]))
  res$FC[res$FC>0] <- (2^(res$FC[res$FC>0]))
  return(res)
}


## gsea part
library(msigdf)
library(dplyr)
library(clusterProfiler)

gsea_mm <- function(res,padj,FC){
  #res: DESeq2 output DEG res, with FC column.
  de <- res$log2FoldChange[res$padj < padj & abs(res$FC)>FC]
  names(de) <- res$gene[res$padj < padj & abs(res$FC)>FC]
  
  de <- sort(de,decreasing = TRUE)
  
  # de: gsea input, sorted log2fc, names are related gene names
  hallmark <- msigdf.mouse %>%
    filter(category_code == "hallmark") %>% select(geneset, mouse.symbol) %>% as.data.frame
  if(length(intersect(names(de),hallmark[,"mouse.symbol"])) == 0){
    gse_ret <- NULL
    gse_num <- NULL
  }
  else{
    gse <- GSEA(de, TERM2GENE = hallmark, pvalueCutoff = 0.05)
    gse_ret <- gse@result[order(gse@result$NES,decreasing = TRUE),]
    gse_num <- dim(gse@result)[1]
  }
  return(list(SUM = paste0("GSEA pathways: ", gse_num), gsea = gse_ret))
}

##

# old edger script
#   20210304 add TMM counts
run_MEedgeR <- function(MAT,ED,WT,n1,n2,lcm=30,padj=0.05,lfc=0.4,gsea=TRUE){
  
  # load requied packages
  library(edgeR)
  library(msigdf)
  library(dplyr)
  library(clusterProfiler)
  
  # load matrix
  mat <- MAT
  
  # filter genes with low expression
  #mat <- mat[rowMeans(mat)>2,]
  
  # build group_list
  group_list <- c(rep(ED,n1),rep(WT,n2))
  group_list <- factor(group_list,levels = c(ED,WT))
  
  # build design
  design <- model.matrix(~0+group_list)
  rownames(design) <- colnames(mat)
  colnames(design) <- levels(group_list)
  
  # DGE (Digital Gene Expression data)
  DGElist <- DGEList(counts = mat, group = group_list)
  # keep genes with CPM cutoff
  # here ignore this par
  #keep_gene <- rowSums(cpm(DGElist)>lcm) >= 3
  #DGElist <- DGElist[keep_gene,, keep.lib.sizes = FALSE]
  DGElist <- calcNormFactors(DGElist, method = "TMM")
  
  #png(paste0(ED,"vs","WT",".MDS.png"))
  #plotMDS(DGElist,method="bcv", col=as.numeric(DGElist$samples$group))
  #dev.off()
  
  d1 <- estimateGLMCommonDisp(DGElist, design)
  d1 <- estimateGLMTrendedDisp(d1, design)
  d1 <- estimateGLMTagwiseDisp(d1, design)
  
  fit <- glmFit(d1, design)
  ret_edgeR <- glmLRT(fit, contrast = c(1,-1))
  DEG_edgeR <- topTags(ret_edgeR, n=nrow(d1))
  DEG_edgeR <- as.data.frame(DEG_edgeR)   
  # output1-edgeR
  
  #
  deGenes <- decideTestsDGE(ret_edgeR, adjust.method="BH",p.value = padj,lfc = lfc)
  deGenes <- rownames(ret_edgeR)[as.logical(deGenes)]
  
  png(paste0(ED,"vs",WT,".Smear.png")) 
  plotSmear(ret_edgeR, de.tags = deGenes)
  abline(h=c(-lfc,lfc),col=2)
  dev.off()
  
  DEG_edgeR_sig <- DEG_edgeR[DEG_edgeR$FDR < padj & DEG_edgeR$logCPM > 1 &abs(DEG_edgeR$logFC) > lfc, ]
  DEG_edgeR_sig <- DEG_edgeR_sig[order(DEG_edgeR_sig$logFC, decreasing = TRUE), ] 
  # output2-dedgeR
  
  ## for GSEA
  # DGEs sorted by logFC
  de <- DEG_edgeR_sig$logFC
  names(de) <- rownames(DEG_edgeR_sig)
  
  
  
  if(gsea==TRUE){
    # get catalogue
    hallmark <- msigdf.mouse %>%
      filter(category_code == "hallmark") %>% select(geneset, mouse.symbol) %>% as.data.frame
    #filter(category_code == "c5") %>% select(geneset, mouse.symbol) %>% as.data.frame
    if(length(intersect(names(de),hallmark[,"mouse.symbol"])) == 0){
      gse_ret <- NULL
      gse_num <- NULL
    }
    else{
      gse <- GSEA(de, TERM2GENE = hallmark, pvalueCutoff = 0.05)
      gse_ret <- gse@result[order(gse@result$NES,decreasing = TRUE),]
      gse_num <- dim(gse@result)[1]
    }
  } else{
    gse_ret <- NULL
    gse_num <- NULL
  }
  # summary
  SUM <- rbind(paste0(ED," vs ",WT),
               paste0("DEGs: ",length(de)),
               paste0("  up: ",sum(de>0)),
               paste0("down: ",sum(de<0)),
               paste0("GSEA pathways: ",gse_num))
  
  # rename edger result to match deseq2 function
  edg2des <- function(retr){
    colnames(retr) <- c("log2FoldChange","logCPM","LR","pvalue","padj")
    #retr <- as.numeric(retr)
    return(retr)
  }
  
  DEG_edgeR <- edg2des(DEG_edgeR)
  DEG_edgeR_sig <- edg2des(DEG_edgeR_sig)
  
  # whole output results
  return(list(DEG_edgeR = DEG_edgeR,DEG_edgeR_sig = DEG_edgeR_sig, SUM = SUM,
              de = de, gsea = gse_ret, TMM = DGElist))
}

##  


## to calculate rowratio for a matrix
#
rowRatio <- function(df){
  for(i in 1:dim(df)[1]){
    df[i,] <- df[i,]/rowSums(df)[i]
  }
  return(df)
}


## modify heatmap with standard deepblue/white/deepred
#

color.test <- colorRampPalette(
  c(
    "#03047F", # deep blue
    "white",
    "#CC2627"  # deep red
  )
)(100)

#

## proc_DEG()
#       to process edgeR result for DEGs' comparison
#              mat_cut is a matrix after advanced filtering, 'like TPM > n in at least one condition'
# 

proc_DEG <- function(deg, p.cut=0.05, FC.cut = 2, padj=TRUE, abs=TRUE, mat_cut=NULL, gene_cut=NULL){
  rownames(deg) <- deg$gene
  
  if(padj==TRUE){
    deg <- deg %>% filter(padj < p.cut)
  }else{
    deg <- deg %>% filter(pvalue < p.cut)
  }
  
  if(abs==TRUE){
    deg <- deg %>% filter(abs(FC) > FC.cut)
  }else if(FC.cut >0){
    deg <- deg %>% filter(FC > FC.cut)
  }else{
    deg <- deg %>% filter(FC < FC.cut)
  }
  
  if(!is.null(mat_cut)){
    deg <- deg[rownames(deg) %in% rownames(mat_cut),]
  }
  if(!is.null(gene_cut)){
    deg <- deg[rownames(deg) %in% gene_cut,]
  }
  return(deg)
}






