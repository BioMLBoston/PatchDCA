
#get distance between two atoms
# pdb1 and pdb2 are two lines of pdb file which is an object of bio3d package
getDistance=function(pdb1,pdb2){
  return(  sqrt((pdb2$x-pdb1$x)^2+(pdb2$y-pdb1$y)^2+(pdb2$z-pdb1$z)^2))
  
}

#convert an MSA to binary matrix representation. every amino acid represented by a vector of size 21. Then the covariance of the binary
#matrix is calculated
convertToNumericalMSAMOdifies=function(X)
{
  LetterToNumberTable=c('-','A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y')
  # get length of sequence after triming white Spaces
  seq=strsplit(trimws(as.character( X[1,])),split="")[[1]]
  numberofSeq=dim(X)[1]
  SequenceLength=length(seq)
  NumericalX=matrix(0, ncol =21*SequenceLength ,nrow=numberofSeq)
  for(i in 1:numberofSeq){
    
    tmpseq=strsplit(toupper(trimws(as.character( X[i,]))),split="")[[1]]
    for(j in 1:SequenceLength){
      index=which(LetterToNumberTable==tmpseq[j])
      binvectory=rep(0,21)
      binvectory[index]=1
      replaceMentIndex=seq(((j-1)*21+1),j*21,1)
      NumericalX[i,replaceMentIndex]=binvectory
    }
    
  }
  
  
  return (cov(NumericalX))
}

#calculte euclidean distance between every two residues within a protein as a matrix
GetDistanceMatrix=function(pdbFile1,pdbFile2,P1Distance,P2Distance)
{
  for(r in 1:length(p1Index)){
    pdb1=pdbFile1[ pdbFile1$resno==as.numeric(rownames(P1Distance)[r]) &pdbFile1$elety=="CB" & pdbFile1$type=="ATOM",]
    if(nrow(pdb1)==0)
    {
      pdb1=pdbFile1[  pdbFile1$resno==as.numeric(rownames(P1Distance)[r] ) &pdbFile1$elety=="CA" & pdbFile1$type=="ATOM",]
    }
    for(cz in 1:length(p1Index)){
      pdb2=pdbFile1[pdbFile1$resno==as.numeric(rownames(P1Distance)[cz]) &pdbFile1$elety=="CB" & pdbFile1$type=="ATOM",]
      if(nrow(pdb2)==0)
      {
        pdb2=pdbFile1[pdbFile1$resno==as.numeric(rownames(P1Distance)[cz]) &pdbFile1$elety=="CA" & pdbFile1$type=="ATOM",]
        
      }
      P1Distance[r,cz]=getDistance(pdb1,pdb2)
    }
  }
  write.table(P1Distance, file = paste0("data/DistanceMatrix/",p1Name,"-",files[z],".distance"))
  #protein2
  for(r in 1:length(p2Index)){
    pdb1=pdbFile2[pdbFile2$resno==as.numeric(rownames(P2Distance)[r]) &pdbFile2$elety=="CB" & pdbFile2$type=="ATOM",]
    if(nrow(pdb1)==0)
    {
      pdb1=pdbFile2[pdbFile2$resno==as.numeric(rownames(P2Distance)[r]) &pdbFile2$elety=="CA" & pdbFile2$type=="ATOM",]
    }
    for(cz in 1:length(p2Index)){
      pdb2=pdbFile2[ pdbFile2$resno==as.numeric(rownames(P2Distance)[cz]) &pdbFile2$elety=="CB" & pdbFile2$type=="ATOM",]
      if(nrow(pdb2)==0)
      {
        pdb2=pdbFile2[ pdbFile2$resno==as.numeric(rownames(P2Distance)[cz]) &pdbFile2$elety=="CA" & pdbFile2$type=="ATOM",]
        
      }
      P2Distance[r,cz]=getDistance(pdb1,pdb2)
    }
  }
  write.table(P2Distance, file = paste0("data/DistanceMatrix/",p2Name,"-",files[z],".distance"))
}

#get the Gaussian Filter
GaussianFilter=function(sigma = 0.5,W = 3){
  a1=c(0.0113437, 0.0838195 ,0.0113437)
  a2=c(0.0838195 ,0.619347, 0.0838195)
  a3=c(0.0113437, 0.0838195, 0.0113437)
  kernel =rbind(a1,a2,a3)
  return(kernel)
}


#calculate the binding site of given complex
getContacts=function(fileName,cutoff=12){
  pdbObj=read.pdb(paste0("data/complexes/",substr(fileName,1,4),".pdb"),multi=T)
  pdbObj=pdbObj$atom
  pdb1=pdbObj[pdbObj$type=="ATOM" & pdbObj$elety=="CA" & pdbObj$chain==substr(fileName,6,6),]
  pdb2=pdbObj[pdbObj$type=="ATOM" & pdbObj$elety=="CA" & pdbObj$chain==substr(fileName,13,13),]
  ContactMatrix=matrix(100,nrow=nrow(pdb1),ncol=nrow(pdb2))
  contactDF=data.frame(pdbID=substr(fileName,1,4),Chains=paste0(substr(fileName,6,6),"_",substr(fileName,13,13)),residue1=NA,residue2=NA,distance=NA,r1Type=NA,r2Type=NA,stringsAsFactors = F)
  rownames(ContactMatrix)=pdb1$resid
  colnames(ContactMatrix)=pdb2$resid
  
  for(j in 1:nrow(pdb1)){
    k=j+1
    while(k<=nrow(pdb2)){
      distance=getDistance(pdb1[j,],pdb2[k,])
      if(distance<=cutoff){
        tmpDF=data.frame(pdbID=substr(fileName,1,4),Chains=paste0(substr(fileName,6,6),"_",substr(fileName,13,13)),residue1=pdb1$resno[j],residue2=pdb2$resno[k],distance=distance,r1Type= pdb1$resid[j],r2Type= pdb2$resid[k],stringsAsFactors = F)
        contactDF=rbind(contactDF,tmpDF)
      }
      k=k+1
    }
  }
  
  contactDF=contactDF[-1,]
  write.table(contactDF,file=paste0("data/complexes/contacts/",fileName,".contact"))
  
}

getModules=function(fileName,dockingfiles){
  moduleDF=data.frame(model=NA,pdbID=NA,Chains=NA,residue1=1,residue2=1,distance=0,stringsAsFactors = F)
  for(i in 1:length(dockingfiles)){
    pdbObj=read.pdb(paste0("data/docking/",fileName,"/",dockingfiles[i]),multi=T)
    pdbObj=pdbObj$atom
    pdb1=pdbObj[pdbObj$type=="ATOM" & pdbObj$elety=="CA" & pdbObj$chain==substr(fileName,6,6),]
    pdb2=pdbObj[pdbObj$type=="ATOM" & pdbObj$elety=="CA" & pdbObj$chain==substr(fileName,13,13),]
    for(j in 1:nrow(pdb1)){
      k=j+1
      while(k<=nrow(pdb2)){
        distance=getDistance(pdb1[j,],pdb2[k,])
        if(distance<=12){
          tmpDF=data.frame(model=dockingfiles[i],pdbID=substr(fileName,1,4),Chains=paste0(substr(fileName,6,6),"_",substr(fileName,13,13)),
                           residue1=pdb1$resno[j],residue2=pdb2$resno[k],distance=distance,stringsAsFactors = F)
          
          moduleDF=rbind(moduleDF,tmpDF)
        }
        k=k+1
      }
    }
  }
  moduleDF=moduleDF[-1,]
  write.table(moduleDF,file=paste0(fileName,".Voting"))
}


#convert prior matrix to list as
# residue1 residue2 penalty
# i j 0.012
# M is the maximum value in the penalty Matrix
# m is the minimum value in the penalty Matrix
ConvertPenaltyMatrixToList=function( M,m,Chain1,Chain2,mapFile,InitialPenalty){
  maxm=m*30
  size=nrow(mapFile)
  Findex=grep(paste0("_",Chain1),mapFile$proteins)
  
  cutoff=max(Findex)
  
  rhoMatrix=data.frame(iIndex=sort(rep((1:(size)),size)),jIndex=rep((1:(size)),size),score=(M+0.02),stringsAsFactors = F)
  for(r in 1:nrow(InitialPenalty))
  {
    j=r+1
    while(j < ncol(InitialPenalty)){
      rhoMatrix$score[rhoMatrix$iIndex==r & rhoMatrix$jIndex==(j+cutoff)]=rhoMatrix$score[rhoMatrix$jIndex==r & rhoMatrix$iIndex==(j+cutoff)]=m+(1-(InitialPenalty[r,j]-min(InitialPenalty))/(max(InitialPenalty)-min(InitialPenalty)))*maxm
      j=j+1
    }
  }
  rhoMatrix$score[rhoMatrix$iIndex<=cutoff & rhoMatrix$jIndex<=cutoff]=M
  rhoMatrix$score[rhoMatrix$iIndex>cutoff & rhoMatrix$jIndex>cutoff]=M
  
  return(rhoMatrix)
  
}


GuassianFilterOnDocking=function(dockingFiles){

  for(i in 1:length(dockingFiles)){
    
    dockingObj=read.table(file=dockingFiles[i],header = T,stringsAsFactors = F)
    smoothObj=read.table(file=paste0("data/docking/smooth/",dockingFiles[i],".voting"),header = T,stringsAsFactors = F)
    m=mapply(dockingObj,FUN=as.numeric)
    m[is.na(m)]=0
   
    kernel=GaussianFilter()
    
    n=ncol(dockingObj)
    m=nrow(dockingObj)
    for(j in 2:(n-1)){
      smoothObj[1,j]=sum(dockingObj[1:2,(j-1):(j+1)]*kernel[2:3,])
      smoothObj[m,j]=sum(dockingObj[(m-1):m,(j-1):(j+1)]*kernel[1:2,])
      
    }
    smoothObj[1,1]=sum(dockingObj[1:2,1:2]*kernel[2:3,2:3])
    smoothObj[1,n]=sum(dockingObj[1,(n-1):(n)]*kernel[2:3,1:2])
    
    smoothObj[m,n]=sum(dockingObj[(m-1):m,(n-1:n)]*kernel[1:2,1:2])
    smoothObj[m,1]=sum(dockingObj[(m-1):m,1:1]*kernel[1:2,2:3])
    
    
    for(j in 2:(m-1)){
      smoothObj[j,1]=sum(dockingObj[(j-1):(j+1),1:2]*kernel[1:3,2:3])
      smoothObj[j,n]=sum(dockingObj[(j-1):(j+1),(n-1):(n)]*kernel[1:3,1:2])
    }
    write.table(smoothObj,file=paste0("data/docking/smooth/",dockingFiles[i],".smooth2"))
    
  }
  
  
}


