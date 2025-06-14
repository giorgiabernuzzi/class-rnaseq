library(DESeq2)
library(tximport)
library(tidyverse)
library(pheatmap)
library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db)

###################################
## PREPARE DATASET CONDITIONS #####
###################################

setwd("/workspaces/class-rnaseq/analysis")

dataset <- tibble(
  sample = c("sample_01",
             "sample_02",
             "sample_03",
             "sample_04",
             "sample_05",
             "sample_06"),
  condition = c(rep("control", 3),
                rep("case", 3))
)
tx2gene <- read_tsv("/workspaces/class-rnaseq/datasets_reference_only/trascriptome/gencode.v29.transcripts_no-vers_chr21_tx2gene.txt")


###################################
#### READ LOCAL FILES IN ##########
###################################

files <- file.path("/workspaces/class-rnaseq/analysis/reads/", paste0(dataset$sample,".quant"), "quant.sf")
names(files) <- dataset$sample

txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

colnames(txi$counts)
rownames(dataset) <- colnames(txi$counts)

dds <- DESeqDataSetFromTximport(txi, dataset, ~condition)

###################################
## PREFILTER MIN COUNTS >10 #####
###################################

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

### make sure base level is control
dds$condition <- relevel(dds$condition, ref = "control")


###################################
##### DIFFERENTIAL EXPRESSION #####
###################################

dds <- DESeq(dds)


###################################
## EXTRACT ANALYSIS RESULTS #####
###################################

res <- results(dds)
resOrdered <- res[order(res$pvalue),]

## writeLines(summary(res), "differential_expression_summary.txt")

plotMA(res, ylim=c(-3,3))

plotDispEsts(dds)

plotCounts(dds, gene=which.min(res$padj), intgroup="condition")

###################################
## WRITE RESULTS OF ANALYSIS #####
###################################

resdata <- as_tibble(resOrdered)
resdata$gene <- rownames(resOrdered)
write_tsv(resdata, "analysis_results.tsv")

save.image("deseq2_analysis.RData")


#### on bash
# cd /workspaces/class-rnaseq/analysis
# git add deseq2_analysis.RData 
# git commit -m "saving analysis"
# git push

library(DESeq2)
library(tximport)
library(tidyverse)
library(pheatmap) ### this might have to be installed
library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db)

setwd("/workspaces/class-rnaseq/analysis")
load("deseq2_analysis.RData")


############################################
## CLUSTERING ##############################
############################################

ntd <- normTransform(dds)
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("condition")])

pheatmap(assay(ntd)[select,],
         cluster_cols=FALSE, annotation_col=df$condition)

plotPCA(ntd, intgroup=c("condition"))

save.image("deseq2_analysis.RData")

library(DESeq2)
library(tximport)
library(tidyverse)
library(pheatmap)
library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db)

setwd("/workspaces/class-rnaseq/analysis")
load("deseq2_analysis.RData")


###################################
## EXTRACT SIGNIFICANT GENES #####
###################################

universe <- AnnotationDbi::select(org.Hs.eg.db,
                                  keys = keys(org.Hs.eg.db),
                                  columns = c('ENTREZID','SYMBOL','ENSEMBL','ENSEMBLTRANS'),
                                  keytype = 'ENTREZID')

sig_genes <- resdata$gene[which(resdata$padj<0.05)]
entrez_genes_sig <- unique(universe[which(universe$ENSEMBL %in% sig_genes),]$ENTREZID)

pvalue_ens_genes <- resdata$padj[which(resdata$padj<0.05)]
names(pvalue_ens_genes)<-sig_genes

pvalue_entrez_genes <- resdata$padj[which(resdata$padj<0.05)]
names(pvalue_entrez_genes) <- entrez_genes_sig


###################################
## ENRICH GO ANALYSIS #####
###################################

ego <- enrichGO( gene = sig_genes,
                 universe = unique(tx2gene$GENEID),
                 OrgDb = org.Hs.eg.db,
                 keyType = 'ENSEMBL',
                 ont = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05)

## writeLines(summary(ego), "enrich_GO_results.txt")

pdf("plots_ego.pdf")
dotplot(ego, showCategory=30)
dev.off()

pdf("plots_network-ego.pdf")
cnetplot(ego, foldChange=resdata$log2FoldChange[which(resdata$padj<0.5)])
dev.off()

###################################
## DISGNET ANALYSIS #####
###################################

## we need to unpack the file we have to read first
## use the terminal
## cd /home/rstudio/data/datasets_class/reference/trascriptome/
## gunzip all_gene_disease_associations.tsv.gz

gda <- read_tsv(gzfile("/workspaces/class-rnaseq/datasets_reference_only/trascriptome/all_gene_disease_associations.tsv.gz"))

disease2gene=gda[, c("diseaseId", "geneId")]
disease2name=gda[, c("diseaseId", "diseaseName")]

disgnet = enricher(entrez_genes_sig, TERM2GENE=disease2gene, TERM2NAME=disease2name)

## writeLines(summary(disgnet), "summary_disgnet.txt")

cnetplot(disgnet, foldChange=resdata$log2FoldChange[which(resdata$padj<0.5)])


save.image("deseq2_analysis.RData")
