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

getwd()
setwd("/home/gitpod/datiesame/analysis")

#####creazione di una tibble che mi rappresenti il mio design sperimentale

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

####con salmon quantifichiamo i trascritti, adesso vogliamo calcolare i livelli di espressione a livello genico
##abbiamo bisogno di una tabella di corrispondenza gene-trascritto

tx2gene <- read_tsv("/home/gitpod/datiesame/datasets_reference_only/trascriptome/gencode.v29.transcripts_no-vers_chr21_tx2gene.txt")


###################################
#### READ LOCAL FILES IN ##########
###################################
######adesso dobbiamo importare i dati quantificati da salmon e sommare le conte per ottenere la quantificazione a livello di gene
###si combinano le conte (cioè i livelli di espressione) di tutti  i trascritti che appartengono a quel gene
### otteniamo le conte totali per ogni gene
####creiamo un vettore coi percorsi dei file da leggere (percorso ai file di output di Salmon)
##a ciascun elemento del vettore assegno un nome corrispondente al campione

files <- file.path("/home/gitpod/datiesame/analysis/reads/", paste0(dataset$sample,".quant"), "quant.sf")
names(files) <- dataset$sample


######lettura dei file di quantificazione di salmon e aggregazione delle conte per gene
###fornisco l'oggetto tx2gene al software che permette di aggregare le conte a livello di gene partendo dai dati per trascritto


txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

#ogni riga ho un gene, ogni colonna un campione (conte per ciascun campione)

colnames(txi$counts)
rownames(dataset) <- colnames(txi$counts)

########funzione per creare un datasetDeSeq che combina le conte per gene (cioè i dati ricalcolati ed aggregati per gene) e il design dell’esperimento (a quale gruppo appartiene ciascun campione)
####IMPORTANTE: condition definisce il modello sperimentale che stiamo usando—> a sinistra variabile dipendente
####a destra la variabile indipendente che DISTINGUE I DUE GRUPPI SPERIMENTALI
###QUESTA FORMULA DICE CHE SUPPONIAMO CHE L’ESPRESSIONE DIFFERENZIALE (DE) DIPENDA DAL FATTORE CONDITION

dds <- DESeqDataSetFromTximport(txi, dataset, ~condition)


###################################
## PREFILTER MIN COUNTS >10 #####
###################################
###per ogni gene (riga della matrice delle conte) sommiamo i valori di tutte le colonne (somma di tutte delle conte di tutti i campioni)
###keep contiene un vettore logico true o false che indica quali geni soddisfano il criterio
###filtro l’oggetto dds, mantenendo solo i geni che soddisfano il criterio

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

### make sure base level is control
###specificare chiaramente quale gruppo è il riferimento 
### funzione relevel: imposta il livello di riferimento per il fattore, questo caso controllo

dds$condition <- relevel(dds$condition, ref = "control")


###################################
##### DIFFERENTIAL EXPRESSION #####
###################################
####per fare l'analisi di espressione differenziale si utilizza la funzione DESeq

dds <- DESeq(dds)

###questi dati li vogliamo ora in forma di tabella
###################################
## EXTRACT ANALYSIS RESULTS #####
###################################
###con la funzione results stiamo estraendo i risultati della DE e li stiamo mettendo nell'oggetto res
###li ordiniamo poi in base alla p-value in resOrdered

res <- results(dds)
resOrdered <- res[order(res$pvalue),]

head(resOrdered)

## writeLines(summary(res), "differential_expression_summary.txt")

plotMA(res, ylim=c(-3,3))

###questo è il plot che ci fa vedere la fit sulle nostre dispersioni e lo shrinkage

plotDispEsts(dds)

###utilizzando il plot delle counts posso fare un plot per ciascun gene
##qui sto scegliendo quello con la p-value piu bassa

plotCounts(dds, gene=which.min(res$padj), intgroup="condition")

###################################
## WRITE RESULTS OF ANALYSIS #####
###################################
#####converte l'oggetto di resOrdered in una tibble ed aggiunge una colonna coi nomi dei geni

resdata <- as_tibble(resOrdered)
resdata$gene <- rownames(resOrdered)
write_tsv(resdata, "analysis_results.tsv")

##salvo i risultati in un file in formato tsv, in cui campi sono separati da tabulazioni
save.image("deseq2_analysis.RData")




############################################
## CLUSTERING ##############################
############################################
###questa funzione applica una trasformazione di normalizzazione ai dati del dds (deseq-dataset)
####I DATI DI ESPRESSIONE VENGONO NORMALIZZATI PER CORREGGERE LE DIFFERENZE TRA CAMPIONI
###SELEZIONE DEI 20 GENI CON MEDIA DELLE CONTE PIU ALTA --> GENI PIU ESPRESSI
###################ORDINA I GENI IN BASE ALLA LORO MEDIA DELLE CONTE, DAL PIU ALTO AL PIU BASSO
######SELEZIONA I GENI CON LA MEDIA DELLE CONTE PIU ALTA 

ntd <- normTransform(dds)
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("condition")])

pheatmap(assay(ntd)[select,],
         cluster_cols=FALSE, annotation_col=df$condition)


###esegue una component analysis sui dati trasformati (LOG2FC delle conte normalizzate)
###la PCA aiuta a vedere la separazione tra i gruppi, riducendo la dimensionalita
###le componenti principali spiegano la maggior parte della varianza tra campioni

plotPCA(ntd, intgroup=c("condition"))

save.image("deseq2_analysis.RData")





###################################
## EXTRACT SIGNIFICANT GENES #####
###################################
##ESTRAZIONE DELL'UNIVERSO DA DATABASE
###sono presenti tutti i geni del genoma umano ed estraiamo ID ENTREZID,SYMBOL, ID ENSEMBL E ID TRASCRITTO


universe <- AnnotationDbi::select(org.Hs.eg.db,
                                  keys = keys(org.Hs.eg.db),
                                  columns = c('ENTREZID','SYMBOL','ENSEMBL','ENSEMBLTRANS'),
                                  keytype = 'ENTREZID')

head(universe)


#####a questo punto facciamo un test di overrappresentazione, serve:
###l'universo
###lista di geni differential expressed
###funzione (che è il pacchetto cluster profile) che prende tutte le liste categoria-specifiche di Gene Onthology
###le testa una ad una contro la lista di geni


###lista di geni differential expressed

sig_genes <- resdata$gene[which(resdata$padj<0.05)]

##il comando che segue serve solo per convertire i nomi dei geni e prendere quelli significativi

entrez_genes_sig <- unique(universe[which(universe$ENSEMBL %in% sig_genes),]$ENTREZID)
pvalue_ens_genes <- resdata$padj[which(resdata$padj<0.05)]
names(pvalue_ens_genes)<-sig_genes

pvalue_ens_genes

pvalue_entrez_genes <- resdata$padj[which(resdata$padj<0.05)]
names(pvalue_entrez_genes) <- entrez_genes_sig

pvalue_entrez_genes


###################################
## ENRICH GO ANALYSIS #####
###################################
######la funzione enrichGO va a testare la lista di geni facendo un test ipergeometrico contro tutte le funzioni della categoria che indico
###con ONT scelgo la categoria di Gene Onthology
###si ha la conta di geni che sono membri di quella funzione, poi ho il rapporto 
####(ovvero quanti geni della mia lista sono membri di quella funzione)
###se quest'analisi non è significativa avremo 0 geni

ego <- enrichGO( gene = sig_genes,
                 universe = unique(tx2gene$GENEID),
                 OrgDb = org.Hs.eg.db,
                 keyType = 'ENSEMBL',
                 ont = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05)

ego

## writeLines(summary(ego), "enrich_GO_results.txt")

#pdf("plots_ego.pdf")
dotplot(ego, showCategory=10)
#dev.off()

######alcuni geni della mia analisi di enrichment saranno annotati per piu funzioni
###vado a vedere quali geni partecipano a quali funzioni
#####LE CONTE SI RIFERISCONO AI NOMI DEI PROCESSI
####OSSERVARE QUALI PROCESSI HANNO PIU GENI CHE PARTECIPANO TRA QUELLI SIGNIFICATIVI

#pdf("plots_network-ego.pdf")
cnetplot(ego, foldChange=resdata$log2FoldChange[which(resdata$padj<0.5)])
#dev.off()

###################################
## DISGNET ANALYSIS #####
###################################

## we need to unpack the file we have to read first
## use the terminal
## cd /home/rstudio/data/datasets_class/reference/trascriptome/
## gunzip all_gene_disease_associations.tsv.gz



#########GDA - DATABASE GENE-DISEASE ASSOCIATION: RACCOLTA DEI RISULTATI IN LETTERATURA
###DAL PUNTO DI VISTA GENE-MALATTIA
##avere un idea su che e quanti geni tra quelli significativi sono gia stati riportati come 
###associati a malattie
gda <- read_tsv(gzfile("/home/gitpod/datiesame/datasets_reference_only/trascriptome/all_gene_disease_associations.tsv.gz"))

disease2gene=gda[, c("diseaseId", "geneId")]
disease2name=gda[, c("diseaseId", "diseaseName")]

disgnet = enricher(entrez_genes_sig, TERM2GENE=disease2gene, TERM2NAME=disease2name)

## writeLines(summary(disgnet), "summary_disgnet.txt")

cnetplot(disgnet, foldChange=resdata$log2FoldChange[which(resdata$padj<0.5)])


save.image("deseq2_analysis.RData")
