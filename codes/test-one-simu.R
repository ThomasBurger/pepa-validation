##############################
# Datasets 
##############################

print("simulation starts")
res.SimulData <- dataGenerator(n1, n2, q, p, NumberDiffAb, DiffAbRatio, pepcv, propShared, save=0)
de <- res.SimulData[[3]]
##transfo SimData into an MsnSet
data <- Sim2MSnSet(res.SimulData[[1]],res.SimulData[[2]])

dato <- data
pData(dato)$Label <- as.factor(c("C1","C1","C1","C2","C2","C2"))
fData(dato)$Protein.group.IDs <- as.factor(fData(dato)$Protein.group.IDs)

print("simulation done")







##############################
# Preprocessing
##############################

print("preprocessing starts")
print("adjacency matrix")
## Build design matrix X
X <- as.matrix(BuildAdjacencyMatrix(data, 'Protein.group.IDs', unique = FALSE))
X.spec <- as.matrix(BuildAdjacencyMatrix(data, 'Protein.group.IDs', unique = TRUE))

rm(res.SimulData); gc()
y <- exprs(data)
y <- y[rownames(X), ]

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
msqrob.res <- rep(1, ncol(X))
msqrob.remainingProt <- as.numeric(rownames(msqrob.results))
msqrob.remainingProt.index <- which(!is.na(msqrob.remainingProt))  
msqrob.remainingProt.filtered <- msqrob.remainingProt[msqrob.remainingProt.index]
msqrob.res[msqrob.remainingProt.filtered] <- -log10(msqrob.results[msqrob.remainingProt.index,5])
msqrob.res[which(is.na(msqrob.res))] <- 0
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
