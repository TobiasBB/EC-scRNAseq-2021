
#### Decoding the molecular landscape of the developing spatial processing system and production of entorhinal stellate cell-like cells by directed programming approach
#### Bergmann, TB; Liu, Y; Mogus, L; Lee, J; Pfister, U; Hendfield, LF; Ansenjo-Maretinez, A; Lisa-Vargas, I; Seemann, SE; Lee, JTH; Patikas, N; Kornum, RB; Denham, M; Hyttel, P; Witter, MP; Gorodkin, J; Pers, TH; Hemberg, M; Khodosevich, K; Hall, VJ

#### R code used for analysis of single-cell RNA sequencing data for Figure 1


#Required functions and libraries
library(Seurat) # Seurat_2.3.4
library(SingleCellExperiment) # SingleCellExperiment_1.2.0
library(Matrix) # Matrix_1.2-14 
library(scater) # scater_1.8.0 
library(scran) # scran_1.8.2 
library(dplyr) # dplyr_0.7.6


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




#Initial clustering of all cells

pathlist <- c(E60_3, E70_3, E50_1, E60_2, E60_1, E70_1, E70_2, Adult_1, E70_nc, Adult_nc)

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
  logcounts(newite) <- counts(newite)
  logcounts(newite)@x <- log2(logcounts(newite)@x + 1)
  
  dec <- tryCatch(rownames(decomposeVar(newite, trendVar(newite, parametric=T, use.spikes = FALSE))),error=function(cond){return(c())})
  
  dalist <- c(dalist, newite)
  
  davar <- c(davar, dec)
  batch <- c(batch, rep(aliases[i], ncol(counts(newite))))
  
}

dalist <- do.call(multiBatchNorm, dalist)
davar <- sort(table(as.character(davar)))
nlist <- rownames(bgtable)[(genebgfreq < 0.3)]
nlist <- nlist [!is.na(match(nlist, names(davar)))]

davar <- davar[nlist]

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
csce <- runTSNE(sce, use_dimred="MNN")

sro_all <- makeSeuratFromSingleCellExperiment(csce)
sro_all@data <- logcounts(csce)
sro_all <- FindClusters(object = sro_all, reduction.type = "MNN", dims.use = 1:50, resolution = 0.8)


#Exclude mesodermal cells and cluster final parental dataset
clusters <- c(0, 1, 2, 3, 5, 6, 7,  9, 10, 11, 12, 13, 14, 15,  17, 18, 19, 21, 22, 23, 24)
cellnamesfilter <- c()

for (i in 1:length(clusters)){
  cellnamesfilter <- c(cellnamesfilter, sro_all@cell.names[sro_all@meta.data$res.0.8 == clusters[i]])
}

dofilter <- T

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
csce <- runTSNE(sce, use_dimred="MNN")

sro <- makeSeuratFromSingleCellExperiment(csce)
sro@data <- logcounts(csce)
sro <- FindClusters(object = sro, reduction.type = "MNN", dims.use = 1:50, resolution = 1.6)






