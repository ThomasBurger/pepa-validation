
fusePeptidesHY<-function(obj,pepNum1,pepNum2,index, ListYeastPepToRemov){
  pep1 <- which(rownames(exprs(obj))== as.numeric(pepNum1))
  pep2 <- which(rownames(exprs(obj))== as.numeric(pepNum2)) 
  exprs(obj)[pep1,] <- log2(2^(exprs(obj)[pep1,])+2^(exprs(obj)[pep2,]))
  fData(obj)[pep1,]$Protein.group.IDs <- paste0(fData(obj)[pep1,]$Protein.group.IDs, ";", fData(obj)[pep2,]$Protein.group.IDs)
  fData(obj)[pep1,]$Protein.group.ID
  fData(obj)[pep1,]$Proteins <- "mixed-Human-Yeast"
  ListYeastPepToRemov[index] <- pep2
  res <- list()
  res$list <- ListYeastPepToRemov
  res$data <- obj
  #obj <- obj[-pep2]
  return(res)
}

FuseListOfPepHY<-function(data,da,Nshare){
  if(Nshare){
    ListYeastPepToRemov <- rep(0,Nshare)
    HumanPeptidesLines <- which(data@featureData$Protein.group.IDs %in% da)
    HumanPeptidesNum <-  rownames(exprs(data)[HumanPeptidesLines,])
    YeastPeptidesNum <- rownames(exprs(data)[-HumanPeptidesLines,])
    
    print(c("Number of human peptides:",length(HumanPeptidesNum)))
    print(c("Number of yeast peptides:",length(YeastPeptidesNum)))
    
    HSelPep <- sample(HumanPeptidesNum, size = Nshare, replace = TRUE)
    YPelPep <- sample(YeastPeptidesNum, size = Nshare, replace = FALSE)
    
    pairs <- cbind(HSelPep,YPelPep)
    L <- dim(pairs)[1]
    for(i in 1:L){
      print(c("HY:", i))
      res <- fusePeptidesHY(data,pairs[i,1],pairs[i,2], i, ListYeastPepToRemov)
      ListYeastPepToRemov <- res$list
      data <- res$data
      #print(c(pairs[i,1],pairs[i,2]))
    }
    data <- data[-ListYeastPepToRemov]
  }
  return(data)
}


makeYYpeptideFusion<-function(data,NSYY){
  if(NSYY){
    for(i in 1:NSYY){
      HumanPeptidesLines <- which(data@featureData$Protein.group.IDs %in% 1:60)
      YeastPeptidesNum <- rownames(exprs(data)[-HumanPeptidesLines,])
      if(i==1){
        print(c("Initial number of peptides:", length(rownames(exprs(data)))))
        print(c("Number of fusable yeast peptides:",length(YeastPeptidesNum)))
      }
      #print(c("YY:", i))
      
      pair2fuse <- sample(YeastPeptidesNum, size = 2, replace = FALSE)
      pep1 <- which(rownames(exprs(data))== as.numeric(pair2fuse[1]))
      pep2 <- which(rownames(exprs(data))== as.numeric(pair2fuse[2])) 
      
      exprs(data)[pep1,] <- log2(2^(exprs(data)[pep1,])+2^(exprs(data)[pep2,]))
      fData(data)[pep1,]$Protein.group.IDs <- paste0(fData(data)[pep1,]$Protein.group.IDs, ";", fData(data)[pep2,]$Protein.group.IDs)
      fData(data)[pep1,]$Proteins <- "Yeast (fused)"
      data <- data[-pep2,]
    } 
  }
  print(c("Final number of peptides:", length(rownames(exprs(data)))))
  return(data)
}