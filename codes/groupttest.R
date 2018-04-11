groupttest <- function(MatAdj, expr){
  nProt <- dim(MatAdj)[2]
  res <- rep(0,nProt)
  for(i in 1:nProt){
    #print(i)
    index <- names(which(MatAdj[,i]==1))
    if(length(index)== 0){
      res[i] <- 1
    } else{
      peptidestotest <- expr[which(rownames(expr)%in% index),]
      if(length(index)== 1){
        res[i] <- t.test(x=peptidestotest[1:3], y=peptidestotest[-(1:3)], var.equal=TRUE)$p.value
      } else{
        res[i] <- t.test(x=peptidestotest[,1:3], y=peptidestotest[,-(1:3)], var.equal=TRUE)$p.value
      }
    }
  }
  return(res)
}  