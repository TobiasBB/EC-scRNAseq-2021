
#### Decoding the molecular landscape of the developing spatial processing system and production of entorhinal stellate cell-like cells by directed programming approach
#### Bergmann, TB; Liu, Y; Mogus, L; Lee, J; Pfister, U; Hendfield, LF; Ansenjo-Maretinez, A; Lisa-Vargas, I; Seemann, SE; Lee, JTH; Patikas, N; Kornum, RB; Denham, M; Hyttel, P; Witter, MP; Gorodkin, J; Pers, TH; Hemberg, M; Khodosevich, K; Hall, VJ

#### R code used for analysis of single-cell RNA sequencing data for Figure 3
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



sro_to_subset <- sro_from_parent_dataset #Object from analysis of parent dataset for Figure 1

pathlist <- c(E60_3, E70_3, E50_1, E60_2, E60_1, E70_1, E70_2, Adult_1, E70_nc, Adult_nc)

bgtable # Defined in analysis of parent dataset for Figure 1

clusters <- c(12, 27, 3, 24)
cellnamesfilter <- c()

for (i in 1:length(clusters)){
  cellnamesfilter <- c(cellnamesfilter, sro_to_subset@cell.names[sro_to_subset@meta.data$res.1.6 == clusters[i]])
}

outersect <- function(x, y) {
  sort(c(setdiff(x, y),
         setdiff(y, x)))
}

dofilter <- T

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
csce_oligo <- runTSNE(sce, use_dimred="MNN")

sro_oligo <- makeSeuratFromSingleCellExperiment(csce_oligo)
sro_oligo@data <- logcounts(csce_oligo)

sro_oligo <- FindClusters(object = sro_oligo, reduction.type = "MNN", dims.use = 1:50, resolution = 1.2)

# Identification of markers using scfind
csce_oligo <- csce_oligo[!duplicated(rownames(csce_oligo)), ]
rowData(csce_oligo)$feature_symbol <- rownames(csce_oligo)
csce_oligo <- csce_oligo[!duplicated(rownames(csce_oligo)), ]

csce_oligo$cell_type1 <- sro_oligo@meta.data$res.1.2

geneIndex <- buildCellTypeIndex(csce_oligo, "cell_type1")

sc_found_genes_top5 <-c()
for (i in 1:length(cellTypeNames(geneIndex))){
  sc_found_genes_top5 <- rbind(sc_found_genes_top5, cellTypeMarkers(geneIndex, cellTypeNames(geneIndex)[i], top.k=5, sort.field = "f1"))
}



