# Corrected for Msnbase package version 2.


Sim2MSnSet<-function(y,X)
{
  #id = pept id vector
  id<-c()
  Protein.group.IDs<-c()

  for(i in 1:dim(X)[1]) {

    tmp<-X[i,]
    id[i]<-i
    Protein.group.IDs[i]<-paste(which(tmp!=0), collapse=";")
 
} #end for i

  #fD<-data.frame(cbind(pept, prot), stringsAsFactors = FALSE)

    ## SimData.MSnSet<-new("MSnSet", exprs=y, fData=data.frame(cbind(id, Protein.group.IDs), stringsAsFactors = FALSE))
    SimData.MSnSet<-new("MSnSet", exprs=y)  
  SimData.MSnSet@experimentData@other$typeOfData <- "peptide"
  fData(SimData.MSnSet)<-data.frame(cbind(id, Protein.group.IDs), stringsAsFactors = FALSE)
  ## exprs(SimData.MSnSet)<-y
  
  lab<-c()
  for (j in 1:dim(y)[2]) {
    lab[j]<-paste("sim",j,sep="")
  }
  
  tmp<-data.frame(lab,lab,c(1:dim(y)[2]), rep("", dim(y)[2]), rep("", dim(y)[2]))
  colnames(tmp)<-c("Experiment", "Label", "Bio.Rep", "Tech.Rep", "Analyt.Rep")
  rownames(tmp)<-lab
  #pData(SimData.MSnSet)$Experiment<-lab
  #pData(SimData.MSnSet)$Label<-lab
  #rownames(pData(SimData.MSnSet))<-lab
  
  pData(SimData.MSnSet)<-tmp
  colnames(exprs(SimData.MSnSet))<-lab
  rownames(exprs(SimData.MSnSet))<-fData(SimData.MSnSet)$id

  
  return(SimData.MSnSet)
  
}

# MSnSet with sim data -> MSnSet construction from "scratch"
