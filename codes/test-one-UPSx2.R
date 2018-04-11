##############################
# Datasets 
##############################

print("data loading")
setwd("../preproc_data")

### UPS1x25 or UPS1x2
#dat <- readRDS("UPSpep25-preproc.Rdata")
dat <- readRDS("UPSPepx2Check-preproc.Rdata")

#list of the 48 DA protein
da <- c(2,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56)
#data <- makeYYpeptideFusion(dat, Nshare*5)
data <- FuseListOfPepHY(dat,da,Nshare)
setwd("../codes")
print("shared peptides added")

dato <- data
pData(dato)$Label <- as.factor(c("C1","C1","C1","C2","C2","C2"))
fData(dato)$Protein.group.IDs <- as.factor(fData(dato)$Protein.group.IDs)

##############################
# Preprocessing
##############################

print("preprocessing starts")
print("adjacency matrix")
## Build design matrix X
X <- as.matrix(BuildAdjacencyMatrix(data, 'Protein.group.IDs', unique = FALSE))
X.spec <- as.matrix(BuildAdjacencyMatrix(data, 'Protein.group.IDs', unique = TRUE))

## Expression data
de <- 1*(as.numeric(colnames(X)) %in% da)
y <- exprs(data)
y <- y[rownames(X), ]
n1 <- 3
n2 <- 3
n <- n1+n2

## Keep track of proteins that are lost in aggregation step
unsup.prot <- !(colnames(X) %in% colnames(X.spec))

q <- nrow(X) # Number of peptides
p <- ncol(X) # Number of proteins

print("preprocessing done")

##############################
# Aggregation
##############################

print("two aggregations start")
# sum all specific peptides
tmp <- (pepAgregate(data, "Protein.group.IDs", "sum overall",  matAdj = X.spec))
y.prot.sum.specif <- exprs(tmp); rm(tmp); gc()
print("first aggregation done")
# sum top 3 specific peptides
tmp <- pepAgregate(data, "Protein.group.IDs", "sum on top n",  matAdj = X.spec, n=3)
y.prot.sum.stop3 <- exprs(tmp); rm(tmp,data); gc()
print("second aggregation done")

##############################
# Tests
##############################

print("tests starts")

dd <- c(rep(1, n1), rep(0, n2))
des <- cbind(dd, rev(dd)); colnames(des) <- c('c1', 'c2')

## Student t-test
tt.sum.stop3.tmp <- apply(y.prot.sum.stop3, 1, FUN=function(v) t.test(x=v[1:n1], y=v[-(1:n1)], var.equal=TRUE)$p.value)
tt.sum.stop3.pv <- rep(1, ncol(X))
tt.sum.stop3.pv[!unsup.prot] <- tt.sum.stop3.tmp
tt.sum.stop3.res <- -log10(tt.sum.stop3.pv)
print("Top3-TT done")

# SAM (autofudge)
sam.sum.specif.tmp <- sam(data=y.prot.sum.specif, cl=dd, var.equal=TRUE)@p.value
sam.sum.specif.pv <- rep(1, ncol(X)) 
sam.sum.specif.pv[!unsup.prot] <- sam.sum.specif.tmp
sam.sum.specif.res <- -log10(sam.sum.specif.pv)
print("Allspec-SAM done")

# MSqRob
print("MSqRob: data prep")
msqrob.proteins <- MSnSet2protdata(dato, accession="Protein.group.IDs")
msqrob.fixed <- c("Label")
msqrob.random <- c("id","Experiment")
msqrob.L <- makeContrast(contrasts=c("LabelC2-LabelC1"), levels=c("LabelC1", "LabelC2"))
print("MSqRob: model fitting")
msqrob.models <- fit.model(protdata=msqrob.proteins, response="quant_value", fixed=msqrob.fixed, random=msqrob.random)
print("MSqRob: test")
msqrob.results <- test.contrast_adjust(msqrob.models, msqrob.L, level=0.05)
print("MSqRob: output")
# remove from 'msqrob.results' proteins ID that indicates merging
colXsp <- colnames(X.spec)
msqrob.res <- rep(0, ncol(X))
for(i in 1:ncol(X.spec)){
  msqrob.res[which(colnames(X) %in% colXsp[i])] <- -log10(msqrob.results[as.numeric(which(rownames(msqrob.results) %in% colXsp[i])),5])
}
print("MSqRob done")

# Peptide Spec-based model
peptide.spec.based.tmp <- groupttest(X.spec,y)
peptide.spec.based.pv <- rep(1, ncol(X))
peptide.spec.based.pv[!unsup.prot] <- peptide.spec.based.tmp
peptide.spec.based.res <- -log10(peptide.spec.based.pv)
print("PeptideModel done")

# LRT-ML and LRT-MAP
pepa.res <- pepa.test(X, y, n1, n2)
print("Pepa tests done")

print("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx")
print("All tests done")
print("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx")

gc()
