#### Decoding the molecular landscape of the developing spatial processing system and production of entorhinal stellate cell-like cells by directed programming approach
#### Bergmann, TB; Liu, Y; Mogus, L; Lee, J; Pfister, U; Hendfield, LF; Ansenjo-Maretinez, A; Lisa-Vargas, I; Seemann, SE; Lee, JTH; Patikas, N; Kornum, RB; Denham, M; Hyttel, P; Witter, MP; Gorodkin, J; Pers, TH; Hemberg, M; Khodosevich, K; Hall, VJ

#### R code used for analysis of single-cell RNA sequencing data for Figure 4 and 5C
#### Dependent on object created in analysis of parent dataset in Figure 1


#Required functions and libraries
library(Seurat) #v. 2.3.4
library(SingleCellExperiment) #v. 1.2.0
library(Matrix) #v. 1.2-14
library(scater) #v. 1.8.0
library(DropletUtils) #v. 1.2.2
library(scran) #v. 1.8.2
library(scfind) #v. 3.7.0


makeSeuratFromSingleCellExperiment <- function(sce, celllist =c(), assay.raw.data.name="counts", makeSparse=T) {
  library("SingleCellExperiment");
  library("Seurat");
  library("Matrix");
  if (length(celllist) == 0){
    rawcounts <- Matrix(assays(sce)[[assay.raw.data.name]], sparse=makeSparse)
    if (is.null(rownames(rawcounts))){
      rownames(rawcounts) <- rowData(sce)$feature_symbol
    }else{
      rownames(rawcounts) <- rownames(sce)
    }
    if (is.null(colnames(rawcounts))){
      colnames(rawcounts) <- colData(sce)$feature_symbol
    }else{
      colnames(rawcounts) <- colnames(sce)
    }
  }
  sro <- CreateSeuratObject( rawcounts )
  toadd <-data.frame(row.names = colnames(sce))
  for(i in names(colData(sce))){
    toadd[[i]] <- colData(sce)[[i]]
  }
  if (ncol(toadd) > 0) sro <- AddMetaData(sro, toadd)
  
  for( drn in reducedDimNames(sce)){
    tmp <- reducedDim(sce,drn)
    colnames(tmp) <- paste(drn,1:ncol(tmp),sep="_")
    rownames(tmp) <- rownames(sro@meta.data)
    sro <- SetDimReduction(sro, reduction.type = drn, slot = "cell.embeddings", new.data = tmp)
    sro <- SetDimReduction(sro, reduction.type = drn, slot = "key", new.data = drn)	
  }
  return(sro)
}


#Clustering of oligodendrocyte subset

sro_to_subset <- sro_from_parent_dataset #Object from analysis of parent dataset for Figure 1

pathlist <- c(E60_3, E70_3, E50_1, E60_2, E60_1, E70_1, E70_2, Adult_1, E70_nc, Adult_nc)

bgtable # Defined in analysis of parent dataset for Figure 1

clusters <- c(0, 5, 10, 15, 17, 19) # Clusters including excitatory and IPs
cellnamesfilter <- c()

for (i in 1:length(clusters)){
  cellnamesfilter <- c(cellnamesfilter, sro_to_subset@cell.names[sro_to_subset@meta.data$res.1.6 == clusters[i]])
}

outersect <- function(x, y) {
  sort(c(setdiff(x, y),
         setdiff(y, x)))
}

dofilter <- T

bgtable 

genebgfreq <- rowSums(bgtable[,c("Adult_soupfactor","E70_soupfactor")])/2
names(genebgfreq) <- rownames(bgtable)

genes <- read.table(file = "Adult_nc_r1_p/genes.tsv" , sep = '\t', header = F)  

dalist <- list()
davar <- list()
batch <- c()
aliases <- c()
if (is.null(aliases)) aliases <- sub(".*/","", pathlist)

for( i in 1:length(pathlist)){
  newite <- read10xCounts(pathlist[i], col.names=T)
  colnames(newite) <- paste(aliases[i],colnames(newite),sep="_")
  rownames(newite) <- genes$V2
  
  filter <- !is.na(match(colnames(counts(newite)), cellnamesfilter))
  newite <- SingleCellExperiment(list(counts =counts(newite)[,filter]))
  
  
  logcounts(newite) <- counts(newite)
  logcounts(newite)@x <- log2(logcounts(newite)@x + 1)
  dec <- tryCatch(rownames(decomposeVar(newite, trendVar(newite, parametric=T, use.spikes = FALSE))),error=function(cond){return(c())})
  
  dalist <- c(dalist, newite)
  
  davar <- c(davar, dec)
  batch <- c(batch, rep(aliases[i], ncol(counts(newite))))
}

dalist <- do.call(multiBatchNorm, dalist)
davar <- sort(table(as.character(davar)))
if(dofilter){
  nlist <- rownames(bgtable)[(genebgfreq < 0.3)]
  nlist <- nlist [!is.na(match(nlist, names(davar)))]
  davar <- davar[nlist]
}

if (length(davar) > 5000) davar <- names(davar)[1:5000]

mnninput <- list();
data <- list()
for( i in 1:length(pathlist)){
  mnninput <- c(mnninput, (logcounts(dalist[[i]]))[davar,] )
  data <- c(data, (logcounts(dalist[[i]])))
}
mnnout <- do.call(fastMNN, c(mnninput, list(k=20,d=50,approximate=T)))

sce <- SingleCellExperiment(list(logcounts=do.call(cbind,data)))
reducedDim(sce, "MNN") <- mnnout$corrected

sce$Batch <- batch


sub_counts <- c()
for( i in 1:length(pathlist)){
  sub_counts <- cbind(sub_counts, (counts(dalist[[i]])))
}
counts(sce) <- sub_counts


set.seed(100)
csce_ip.neuron <- runTSNE(sce, use_dimred="MNN")

sro_ip_neuron <- makeSeuratFromSingleCellExperiment(csce_ip.neuron)
sro_ip_neuron@data <- logcounts(csce_ip.neuron)

sro_ip_neuron <- FindClusters(object = sro_ip_neuron, reduction.type = "MNN", dims.use = 1:50, resolution = 1)

# Identification of markers using scfind
sce_ip.neuron <- csce_ip.neuron[!duplicated(rownames(csce_ip.neuron)), ]
rowData(csce_ip.neuron)$feature_symbol <- rownames(csce_ip.neuron)
csce_ip.neuron <- csce_ip.neuron[!duplicated(rownames(csce_ip.neuron)), ]

csce_ip.neuron$cell_type1 <- sro_ip_neuron@meta.data$res.1

geneIndex <- buildCellTypeIndex(csce_ip.neuron, "cell_type1")

sc_found_genes_top5 <-c()
for (i in 1:length(cellTypeNames(geneIndex))){
  sc_found_genes_top5 <- rbind(sc_found_genes_top5, cellTypeMarkers(geneIndex, cellTypeNames(geneIndex)[i], top.k=5, sort.field = "f1"))
}

# Identification of enriched transcription factors in RELN6 and RELN7 populations
tfs <- read.csv("TFcheckpoint3_csv.csv", header = F) # Curated list of transcription factors published by Chawla et al., 2013; PMID 23933972
pig_tfs <- Reduce(intersect, list(tfs[,1], rownames(sro_ip_neuron@data)))

c6_tfs <- FindMarkers(sro_ip_neuron, genes.use = pig_tfs, ident.1 = 6, min.pct = 0.25)
c6_tf_genes <- rownames(c6_tfs[rev(order(c6_tfs$avg_logFC)),])[1:20]

c7_tfs <- FindMarkers(sro_ip_neuron, genes.use = pig_tfs, ident.1 = 7, min.pct = 0.25)
c7_tf_genes <- rownames(c7_tfs[rev(order(c7_tfs$avg_logFC)),])[1:20]






