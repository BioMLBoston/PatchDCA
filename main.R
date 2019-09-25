

##############################
#MAIN
#######################

library(Biostrings)
library(glasso)
library(bio3d)
library(seqinr)
library(corrplot)
library(network)

######################
#read the Docking and Prior information from an input and build the penalty matrix. 
#######################

WorkDirectory="PatchDCA/"

setwd(WorkDirectory)

###variables:
#p1Name: first protein's PDB with chain ID: example 1CLL_A
#p2Name: second protein's PDB with chain ID
#map is matrix which maps every position from MSA to the residue ID from pdb file.
#ComplexSet: Input protein files. See example/proteinsSettings.txt
files=list.files(path="data/map/", pattern=".map")
ComplexSet=read.table("proteinsSettings.txt",stringsAsFactors = F,header = T)

# Read Amino Acid Propesnity Matrix( join probability)
ContacPropensity=read.csv("data/Propensity", header = T,sep = " ")
colnames(ContacPropensity)=rownames(ContacPropensity)=a(colnames(ContacPropensity))
ContacPropensity=ContacPropensity/(sum(ContacPropensity)/2)

#For each row in number of input complexes, it builds one penalty matrix. The penlaty matrix of a complex is saved under data/penalty/ComplexName.penalty.
for(z in 1:nrow(ComplexSet)){
  
  ComplexName=paste0(ComplexSet[i,],collapse = "_")
   
  p1Name=paste0(ComplexSet[i,1:2],collapse = "_")
  p2Name=paste0(ComplexSet[i,3:4],collapse = "_")
  mapFile=read.table(paste0("data/map/",ComplexName,".map"),stringsAsFactors = F,header = T)
  
    
  
  pdbFile1=read.pdb(paste0("data/pdb/",p1Name,".pdb"))
  pdbFile1=pdbFile1$atom
  pdbFile2=read.pdb(paste0("data/pdb/",p2Name,".pdb"))
  pdbFile2=pdbFile2$atom
  
  ComplexNameCov=paste0("data/covariance/",ComplexName,".fas.cov")
  covariance=read.table(ComplexNameCov,header = T,stringsAsFactors = F)
  maxLambda=max(covariance)
  
  
  p1Index=grep(p1Name,mapFile$proteins)
  p2Index=grep(p2Name,mapFile$proteins)
  
  #get distanceMatrix for each protein
  # P1Distance=matrix(-1,nrow=length(p1Index),ncol = length(p1Index)):length(mapFiles)){
  # rownames(P1Distance)=colnames(P1Distance)=mapFile$index[p1Index]
  # P2Distance=matrix(-1,nrow=length(p2Index),ncol = length(p2Index))
  # rownames(P2Distance)=colnames(P2Distance)=mapFile$index[p2Index]
  #
  #caclculate euclidan distance for each protein between all the proteins
  #  GetDistanceMatrix(pdbFile1,pdbFile2,P1Distance,P2Distance)
  
   
  #calculate joint probability
  
  InitialPenalty2=matrix(0.0,nrow=length(p1Index),ncol = length(p2Index))
  InitialPenalty3=matrix(0.0,nrow=length(p1Index),ncol = length(p2Index))
  

  #performe voting between docking complexes from ClusPro
  # first for a complex 
  dockingFile=list.files(paste0("data/docking/",ComplexName,"/"),pattern = ".pdb")
  getModules(ComplexName,dockingFile)
  dockingDF=read.table(paste0("data/docking/smooth/",ComplexName,".docking"),header = T,stringsAsFactors = F)
  dockingDF=dockingDF[dockingDF$distance<=12,]
  for(r in 1:nrow(dockingDF)){
    r1=which(rownames(InitialPenalty2)==dockingDF$residue1[r])
    c1=which(colnames(InitialPenalty2)==dockingDF$residue2[r])
    if(length(r1)>=1 & length(c1)>=1)
    {
      InitialPenalty2[r1,c1]=InitialPenalty2[r1,c1]+1
    }
  }
  write.table(InitialPenalty2, file=paste0("data/docking/smooth/",ComplexName,".distance"))
  
  smoothObj=read.table(file=paste0("data/docking/smooth/",ComplexName,"_8.tab",".smooth2"),header = T,stringsAsFactors = F)
  InitialPenalty2=mapply(smoothObj,FUN=as.numeric)
  InitialPenalty2[is.na(InitialPenalty2)]=0
  InitialPenalty2=(InitialPenalty2)/sort(InitialPenalty2,decreasing = T)[3]
  InitialPenalty2[InitialPenalty2>1]=1
  
  rownames(InitialPenalty2)=rownames(InitialPenalty3)=mapFile$index[p1Index]
  colnames(InitialPenalty2)=colnames(InitialPenalty3)=mapFile$index[p2Index]
  
  
  
#build the propensity prior based on propensity score.
  for(r in 1:nrow(InitialPenalty3))
  {
    for(cz in 1:ncol(InitialPenalty3)){
      aa1=which(rownames(ContacPropensity)==toupper(mapFile$seq[p1Index[r]]))
      aa2=which(colnames(ContacPropensity)==toupper(mapFile$seq[p2Index[cz]]))
      InitialPenalty3[r,cz]=ContacPropensity[aa1,aa2]
        }
  }

  #the coefficients for combining the prior information 
  InitialPenalty=(4*InitialPenalty2+1*InitialPenalty3)/5
 
  #get the penalty matrix

  rhoMatrix=ConvertPenaltyMatrixToList(0.16,0.0001,ComplexSet[i,2],ComplexSet[i,4],mapFile,InitialPenalty)
  write.table(rhoMatrix,paste0("data/penalty/",ComplexName,".penalty"),row.names = F,col.names = F)
  
  
}


########################
#patch building 
##building patches using https://github.com/ACRMGroup/IntPred/blob/master/bin/getPatches.pl
#######################
#after running the PSICOV and getting the output


setwd("/data/benchmark/")
input="normalized/"
# allFiles=list.files(paste0(input),pattern=".psicov.distance$")
# 
# TruthStat=read.table("statAllCont.tab",header = T,stringsAsFactors = F)
bindingCutoff=8
for(i in 1:nrows(ComplexSet)){
  
  patchDF=data.frame(protein=NA,chain=NA,patchMember=0,score=0,center=0,centerDistance=0,stringsAsFactors = F)
  pnameIndex=grep(substr(ComplexSet[i,1]),TruthStat$proteins)
  topObj=read.table(paste0(input,allFiles[i]),stringsAsFactors = F,header = T)
  CLUSTERS=ifelse(TruthStat$Radius8[pnameIndex]*2.5<nrow(topObj),floor(TruthStat$Radius8[pnameIndex]*2.5),nrow(topObj))
  
  colnames(topObj)=c("Residue1","Residue2","V5","psicovDistance")
  for(j in 1:CLUSTERS){
    atom1=topObj$Residue1[j]
    
    chain1=substr(allFiles[i],6,6)
    chain2=substr(allFiles[i],13,13)
    pdbName=substr(allFiles[i],1,4)
    atom2=topObj$Residue2[j]
    command1=paste0("~/bioptools/src/pdbmakepatch -r 4 -m 1 ",chain1,".",atom1, " CA ../patchBuilding/pdbConverted/",substr(allFiles[i],1,6),".pdb.asa.pdb > ",input,"/patch/", pdbName,"_",chain1,"_",atom1,"_",atom2,".pdb")
    command2=paste0("~/bioptools/src/pdbmakepatch -r 4 -m 1 ",chain2,".",atom2, " CA ../patchBuilding/pdbConverted/",substr(allFiles[i],1,4),"_",chain2,".pdb.asa.pdb > ",input,"/patch/", pdbName,"_",chain2,"_",atom1,"_",atom2,".pdb2")
    
    system(command2,wait = T)
    system(command1,wait = T)
    
    pdbObj=read.pdb(paste0(input,"/patch/", pdbName,"_",chain1,"_",atom1,"_",atom2,".pdb"))
    pdbObj=pdbObj$atom
    members=unique(pdbObj$resno[pdbObj$b==1])
    tmp=data.frame(protein=pdbName,chain=chain1,patchMember=paste0(members,collapse = ","),score=topObj$V5[j],center=topObj$Residue1[j],centerDistance=topObj$psicovDistance[j],stringsAsFactors = F)
    patchDF=rbind(patchDF,tmp)
    pdbObj=read.pdb(paste0(input,"/patch/", pdbName,"_",chain2,"_",atom1,"_",atom2,".pdb2"))
    pdbObj=pdbObj$atom
    
    members=unique(pdbObj$resno[pdbObj$b==1])
    tmp=data.frame(protein=pdbName,chain=chain2,patchMember=paste0(members,collapse = ","),score=topObj$V5[j],center=topObj$Residue2[j],centerDistance=topObj$psicovDistance[j],stringsAsFactors = F)
    patchDF=rbind(patchDF,tmp)
    
  }
  
  
  pachDF=patchDF[-1,]
  
  #calculte the patch distance based on modification to jaccard distance
  connectionDF=matrix(0,ncol=CLUSTERS,nrow=CLUSTERS)
  
  tmpPatch=patchDF[patchDF$protein==substr(allFiles[i],1,4) & patchDF$chain %in% c(substr(allFiles[i],6,6),substr(allFiles[i],13,13)),]
  colors=rep("purple",CLUSTERS)
  index=which(topObj$psicovDistance[1:CLUSTERS]<=12)
  colors[index]="green"
  for(j in 1:CLUSTERS){
    atom1=topObj$Residue1[j]
    chain1=substr(allFiles[i],6,6)
    chain2=substr(allFiles[i],13,13)
    pdbName=substr(allFiles[i],1,4)
    atom2=topObj$Residue2[j]
    setA=paste0(chain1,"_",strsplit(tmpPatch$patchMember[2*j-1],split = ",")[[1]])
    setA2=paste0(chain2,"_",strsplit(tmpPatch$patchMember[2*j],split = ",")[[1]])
    
    k=j+1
    while(k<=CLUSTERS){
      setB=paste0(chain1,"_",strsplit(tmpPatch$patchMember[2*k-1],split = ",")[[1]])
      setB2=paste0(chain2,"_",strsplit(tmpPatch$patchMember[2*k],split = ",")[[1]])
      distance=length(intersect(setA,setB))*length(intersect(setA2,setB2))/(length(union(setA,setB))+length(union(setA2,setB2)))
      connectionDF[j,k]= connectionDF[k,j]=distance
      
      k=k+1
    }
    
  }
}

pathMember=apply(connectionDF, 2, sum)

tmp=cbind(topObj[1:(CLUSTERS),],pathMember)
write.table(tmp,paste0("data/patch/",allFiles[i],".patch"),row.names = F)
