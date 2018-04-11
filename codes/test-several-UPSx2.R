# # to get the devel version of DAPAR
# source("http://bioconductor.org/biocLite.R")
# useDevel()
# biocLite("Biobase")
# biocLite("MSnbase")
# biocLite("DAPAR")
# ###################


library(Biobase)
library(MSnbase) # readRDS
library(DAPAR) # requires version 1.9.3 or newer
library(siggenes)
library(limma) # lmFit
library(graph) # prot-pep connected components
library(parallel) # multi-core apply
library(doParallel)
library(lme4)
library(MSqRob)

rm(list = ls())
CEX <- 1.5
LWD <- 2

# local
setwd("D:/Recherche/Projets/regpen")
# server
#setwd("/home/tb232477/") ###!!!###

## Load functions
setwd("./20170619/codes")
source("addSharedPep.R")
source("groupttest.R")


n.simu <- 10

## dataset parameters
n1 <- n2 <- 3
n <- n1+n2 #number of replicates
q<- 5e3
p <- 1e3 #numbers of peptides and proteins
DiffAbRatio <- 15
NumberDiffAb<- 50
pepcv <- 0.1

NshareList <- c(0, 40, 80, 120, 160, 200, 240, 280)

for(Nshare in NshareList){          
  print(c("Number of shared peptides: ", Nshare))
  for(nn in 1:n.simu){
    print(c("iteration: ", nn))
    source('test-one-UPSx2.R')
    
    
    if(nn == 1){
      p <- length(pepa.res$llr.map)
      mat.msqrob <- mat.peptide.spec.based.res <- mat.sam.sum.specif.res <- mat.tt.sum.stop3.res <- mat.de <- matrix(NA, nrow=n.simu, ncol=p)
      mat.llr.ml <- mat.llr.map <- mat.llr.ml.pv <- mat.llr.map.pv <- matrix(NA, nrow=n.simu, ncol=p)
      
    }          
    
    mat.llr.map[nn, ] <- pepa.res$llr.map
    mat.llr.ml[nn, ] <- pepa.res$llr
    mat.msqrob[nn, ] <- msqrob.res
    mat.peptide.spec.based.res[nn, ] <- peptide.spec.based.res
    mat.sam.sum.specif.res[nn, ] <- sam.sum.specif.res
    mat.tt.sum.stop3.res[nn, ] <- tt.sum.stop3.res
    mat.de[nn, ] <- de 
  }
  
  list.scores <- list();
  list.scores <- list('PEPA-MAP' = mat.llr.map,
                      'PEPA-ML'= mat.llr.ml,
                      'MSqRob'= mat.msqrob,
                      #'AllPeptide-test'= mat.all.peptide.based.res,
                      'PeptideModel'= mat.peptide.spec.based.res,
                      'AllSpec-SAM'= mat.sam.sum.specif.res,
                      'Top3Spec-TT'= mat.tt.sum.stop3.res)
  
  setwd("../comparaisons")
  cols <- c('red', 'pink', 'blue', 'lightblue', 'seagreen', 'green')
  types <- c("solid", "solid", "dashed", "dashed", "dotted", "dotted")
  widths <- c(2,2,1.5,1.5,3,3)
  
  ## PR 
  for(ii in 1:length(list.scores)){
    tpr <- pr <- list()
    for(nn in 1:n.simu){
      sc <- list.scores[[ii]][nn, ]
      thr <- sort(sc, decreasing=TRUE)
      tpr[[nn]] <- sapply(thr, FUN=function(tt) mean(sc[mat.de[nn, ] == 1] > tt))
      pr[[nn]] <- sapply(thr, FUN=function(tt) mean(mat.de[nn, sc > tt] == 1))
    }
    tpr.grid <- sort(unique(unlist(tpr)))
    aligned.pr <- matrix(NA, n.simu, length(tpr.grid))
    for(nn in 1:n.simu){
      aligned.pr[nn, ] <- approx(tpr[[nn]], pr[[nn]], tpr.grid)$y
    }
    mean.pr <- colMeans(aligned.pr)
    if(ii == 1){
      pdf(file=sprintf('pr-%.0fNshare-UPSx2.pdf', Nshare))
      par(mar=c(5,5,4,2))
      plot(tpr.grid, mean.pr, type='s', col=cols[ii], xlab='Recall', ylab='Precision', lty = types[ii],
           lwd=widths[ii], cex.lab=CEX, cex.axis=CEX, cex=CEX)
    }else{            
      lines(tpr.grid, mean.pr, type='s', col=cols[ii],lty = types[ii],
            lwd=widths[ii])
    }
    if(ii == length(list.scores)){
      legend('bottomleft', legend=names(list.scores), col=cols[1:length(list.scores)], lty=types, lwd=widths)
      dev.off()
    }        
  } 
  
  
  ## Calibration plots
  s1 <- pepa.res$s1
  
  sam.pv <- 10^-(sam.sum.specif.res)
  peptide.spec.based.pv <- 10^-(peptide.spec.based.res)
  msqrob.pval <- 10^-(msqrob.res)
  llr.ml.pv <- pepa.res$llr.pv
  llr.map.nocalib.pv <- 1 - pchisq(pepa.res$llr.map, 1)
  llr.map.calib.pv <- pepa.res$llr.map.pv
  
  sam.fpr <- sapply(sam.pv, FUN=function(tt) mean((sam.pv<=tt)[de == 0]))
  peptide.spec.based.fpr <- sapply(peptide.spec.based.pv, FUN=function(tt) mean((peptide.spec.based.pv<=tt)[de == 0]))
  msqrob.fpr <- sapply(msqrob.pval, FUN=function(tt) mean((msqrob.pval<=tt)[de == 0]))
  llr.ml.fpr <- sapply(llr.ml.pv, FUN=function(tt) mean((llr.ml.pv<=tt)[de == 0]))  
  llr.map.nocalib.fpr <- sapply(llr.map.nocalib.pv, FUN=function(tt) mean((llr.map.nocalib.pv<=tt)[de == 0]))
  llr.map.calib.fpr <- sapply(llr.map.calib.pv, FUN=function(tt) mean((llr.map.calib.pv<=tt)[de == 0]))
  
  pdf(paste('calibration - s=',s1,', upsx2 - Nshare=', Nshare,'.pdf', sep=""))
  
  
  plot(log10(sam.pv), log10(sam.fpr), main=paste('Calibration plot, s=',s1,', Exp1_R2_pept - Nshare=', Nshare, sep=""), xlab='log10 p-value threshold', ylab='Actual test level (log10)', 
       xlim=c(-5,0), col=cols[5], pch=6)
  points(log10(peptide.spec.based.pv), log10(peptide.spec.based.fpr), col=cols[4], pch=5)
  points(log10(msqrob.pval), log10(msqrob.fpr), col=cols[3], pch=4)
  points(log10(llr.ml.pv), log10(llr.ml.fpr), col=cols[2], pch=3)
  points(log10(llr.map.nocalib.pv), log10(llr.map.nocalib.fpr), col='grey', pch=1)
  points(log10(llr.map.calib.pv), log10(llr.map.calib.fpr), col=cols[1], pch=2)
  
  abline(a=0, b=1, col='black')
  legend('topleft', legend=c('PEPA-MAP', 'PEPA-MAP-RW', 'PEPA-ML', 'MSqRob','PeptideModel', 'AllSpec-SAM'), 
         col=c('grey', cols[1:5]), pch=1:6)
  dev.off()
  
  setwd("../codes")
}


# useDevel(FALSE)
# biocLite() # to reinstall the release version of packages
# biocValid()



