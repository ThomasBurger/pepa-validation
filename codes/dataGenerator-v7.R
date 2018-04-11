##Simulated data generation
#param =
# n1 and n2 : nb of replicates per condition
# q : peptides nb
# p : proteins nb
# NumberDiffAb : nb of differentially abundant proteins
# DiffAbRatio : ratio of abundance difference 
# pepcv : cv for the peptide intensities
# propShared : There are propShared*(q-p)*p non-zero indicator values for the (q-p) others peptides
# save : flag


dataGenerator<- function(n1=3,n2=3,q=5e3,p=1e3,NumberDiffAb=150,DiffAbRatio=5,pepcv=0.1,propShared=0.05, save=0)
{  
# #   parameter initialization to test the function in script mode
#   n1=3
#   n2=3
#   q=5e3
#   p=1e3
#   NumberDiffAb=150
#   DiffAbRatio=5
#   pepcv=0.05
#   propShared=0.05
  
  ## Hidden parameters:
  ######
  #the beta distribution is tuned so that peptides reponse spans over 3 orders of magnitudes
  offset <- 0.001
  alpha <- 1.1
  beta <- 3
  # precision on the simulation
  precPropShared <- 0.05
  
  
  # check shared peptide proportion
  flagBadProp <- T
  while(flagBadProp){
    ## Peptide/protein membership matrix initialization
    X <- matrix(0, nrow=q, ncol=p)
    diag(X) <- (offset+rbeta(p, alpha, beta)) 
    #The p first peptides are protein specific to avoid empty proteins.
    
    ## The other peptides are also protein-specific, yet randomly distributed, and with arbitrary response coefficients
    meanPepProtRatio <- (q-p)/p
    nPepPerProt <- rpois(p,meanPepProtRatio) 
    ProtMembership <- rep(1:p,nPepPerProt)
    delta <- length(ProtMembership) - (q-p) 
    if(delta>0){ProtMembership <- ProtMembership[1:(q-p)]}
    if(delta<0){
      addProtRef <- sample(1:p,-delta)
      ProtMembership <- c(ProtMembership,addProtRef)
    }
    RespCoefs <- (offset+rbeta((q-p), alpha, beta))
    for(i in 1:(q-p)){
      X[p+i,ProtMembership[i]] <- RespCoefs[i]
    }
    
    ## now shared peptides
    if(propShared){
      nSharedPep <- round(propShared*q)
      sharedPep <- sample(1:(q-p),nSharedPep)+p
      for(i in sharedPep){
        emptyCells <- which(X[i,]==0)
        NpepConnect <- 1+rpois(1,0.5)
        if(NpepConnect){
          protToConnect <- sample(emptyCells, NpepConnect)
          X[i,protToConnect] <- (offset+rbeta(NpepConnect, alpha, beta))
        }  
      }
      effectiveSharePeptideRatio <- mean(rowSums(ceiling(X))>1) ## effectiveSharePeptideRatio <- sum(ceiling(X))/q-1
      deltaRatio <- abs(propShared - effectiveSharePeptideRatio)/propShared
      print(c("effective shared-peptide ratio:",effectiveSharePeptideRatio))
      print(c("delta-Ratio:",deltaRatio))
      if(deltaRatio < precPropShared){
        flagBadProp <- F
      }
    } else{
      flagBadProp <- F
    }   
  }
    
  ## Condition-wise protein abundances
  protAb1 <- rlnorm(p,meanlog=17.8,sdlog=0.45)
  #lognorm dist. is tuned so that final peptides log2-intensities are normaly distributed around 23.5
  protAb2<-protAb1
  vectOfRatio <- rnorm(NumberDiffAb,mean=DiffAbRatio,sd=DiffAbRatio*0.05) #5% cv on the ratio
  protAb2[1:NumberDiffAb] <- vectOfRatio*protAb1[1:NumberDiffAb]

  ## Sample-wise protein abundances
  genCV <- rep(0.001, p)
  m1 <- rep(protAb1,n1); m2 <- rep(protAb2,n2)
  sd1 <- rep(genCV*protAb1,n1); sd2 <- rep(genCV*protAb2,n2)
  theta1 <- matrix(rnorm(p*n1, mean=m1, sd=sd1), ncol=n1, byrow=F)
  theta2 <- matrix(rnorm(p*n2, mean=m2, sd=sd2), ncol=n2, byrow=F)
  theta1[which(theta1 <= 0)] <- 0; theta2[which(theta2 <= 0)] <- 0
  
  ## Observed peptide intensities
  Xtheta1 <- log2(X %*% theta1); Xtheta2 <- log2(X %*% theta2)    
  Xtheta1[which(is.infinite(Xtheta1))] <- 0; Xtheta2[which(is.infinite(Xtheta2))] <- 0
  meanPepInt1 <- apply(Xtheta1,1,mean); meanPepInt2 <- apply(Xtheta2,1,mean)
  pepsd1 <- meanPepInt1*rep(pepcv,q); pepsd2 <- meanPepInt2*rep(pepcv,q)
  y1 <- (Xtheta1) + matrix(rnorm(n1*q, mean=rep(0, n1*q), sd=rep(pepsd1,n1)), nrow=q, ncol=n1)
  y2 <- (Xtheta2) + matrix(rnorm(n2*q, mean=rep(0, n2*q), sd=rep(pepsd2,n2)), nrow=q, ncol=n2)
  y <- cbind(y1, y2)
  de <- 1*((1:ncol(X)) <= NumberDiffAb)
  
  if(save==1){
      save(list=c("y","X","de"),
           file=sprintf("Simul_%d_%d_%d_%d_%d--%d_%.1f_%.5f.RData",n1,n2,q,p,NumberDiffAb,DiffAbRatio,pepcv,propShared))
  }
  return(list(y,X,de))
}

