
##Function to compute the covariance between the correlation between x1 and x2 and that between x3 and x4##
CCcov = function(x1,x2,x3,x4,n){
  r12 = cor(x1,x2)
  r13 = cor(x1,x3)
  r14 = cor(x1,x4)
  r23 = cor(x2,x3)
  r24 = cor(x2,x4)
  r34 = cor(x3,x4)
  
  (r12 *r34 *(r13^2 +r14^2 +r23^2 +r24^2)/2 + r13 *r24 +r14 *r23 -
      (r12 *r13 *r14 +r12 *r23 *r24 +r13 *r23 *r34 +r14 *r24 *r34))/n
}


Studies_Summary = function(#This function compute the gene-gene correlations for each study 
  #as well as their covariance matrix if needed
  
    G, #The number of genes consider in each study, we assume it is the same foe all the studies
    Omics_list, #The list of merged multi-omics data of each study. 
    #The multi-omics for each gene should be stack in the same order for all study
    #The genes for each multi-omics should be in the same order
    #Each study multi-omics data should be scaled
    
    numb_omic, #The number of omics per gene
    covM=FALSE, #Whether you want the correlation covariance matrix to be computed for the CC, MOC and the MSOC
    #Change it to true if you aim to implemente the full correlation covariance matrix methods
){
  S = length(Omics_list) # # of studies
  total_individuals = G * numb_omic # total # of omics
  n = rep(NA,S) #The sample size for each study data
  Studies = list()
  
  for (s in 1:S){
    omic = Omics_list[[s]]
    n[s] = dim(omic)[1]
    corM = Diagonal(G) #CC matrix for each study
    corM2 = Diagonal(G) #MOC matrix for each study
    corM3 = Diagonal(G) #SOC matrix for each study
    
    pstar = G*(G-1)/2
    rho1 = rep(NULL, pstar) #The column vector of the CC
    rho2 = rep(NULL, pstar) #The column vector of the MOC
    rho3 = rep(NULL, pstar) #The column vector of the MSOC
    
    if (covM){
      CCdata = list()
      MOCdata = list()
      MSOCdata = list()
    }

    M1 = abs(cor(omic[,seq(1,total_individuals,numb_omic)])) #Cor matrix from omic 1
    M2 = abs(cor(omic[,seq(2,total_individuals,numb_omic)])) #Cor matrix from omic 2
    M3 = abs(cor(omic[,seq(3,total_individuals,numb_omic)])) #Cor matrix from omic 3
    k = 1
    temp_t = Sys.time()
    for (i in 1:(G-1)){
      ##Computing the Maximum Same Omic Correlation
      #This is only the amplitude#
      corM3[i,(i+1):G] = apply(rbind(M1[i,(i+1):G],M2[i,(i+1):G],M3[i,(i+1):G]),2,max)
      Vindex = apply(rbind(M1[i,(i+1):G],M2[i,(i+1):G],M3[i,(i+1):G]),2,which.max)
      corM3[(i+1):G,i] = corM3[i,(i+1):G]
      for (j in (i+1):G){
        ##Computing the canonical correlations
        CC = cancor(omic[,(1:numb_omic)+numb_omic*(i-1)],omic[,(1:numb_omic)+numb_omic*(j-1)])
        corM[i,j] = CC$cor[1]
        corM[j,i] = corM[i,j]
        rho1[k] = CC$cor[1]
        
        if (covM){
          #Keeping the vectors whose Pearson correlation gives the CC
          CCdata[[k]] = list(x1 = omic[,(1:numb_omic)+numb_omic*(i-1)]%*%CC$xcoef[,1],
                             x2 = omic[,(1:numb_omic)+numb_omic*(j-1)]%*%CC$ycoef[,1])
        }
 
        
        ##Computing the Maximum Omic Correlation
        dat = cbind(omic[,(1:numb_omic)+numb_omic*(i-1)],omic[,(1:numb_omic)+numb_omic*(j-1)])
        temp1 = abs(cor(dat)[1:numb_omic,numb_omic+(1:numb_omic)])
        temp_max = max(c(temp1))
        #Keeping the vectors whose Pearson correlation gives the MOC
        index = which.max(c(temp1))
        if (index %% numb_omic == 0){
          temp_x1 = dat[,numb_omic]
          temp_x2 = dat[,numb_omic+index/numb_omic]
        } else {
          temp_x1 = dat[,index %% numb_omic]
          temp_x2 = dat[,(index %/% numb_omic+numb_omic)+1]
        }
        
        if (cor(temp_x1,temp_x2) < 0){ # temp_x1 and temp_x2 are negatively correlated
          temp_max = -temp_max
        }
        corM2[i,j] = temp_max
        corM2[j,i] = corM2[i,j]
        rho2[k] = corM2[i,j]
        
        if (covM){
          MOCdata[[k]] = list(x1 = temp_x1,
                              x2 = temp_x2)
        }
        
        
        #Keeping the vectors whose Pearson correlation gives the MSOC
        rm(temp_x1,temp_x2)
        temp_x1 = dat[,Vindex[j-i]]
        temp_x2 = dat[,numb_omic+Vindex[j-i]]
        if (cor(temp_x1,temp_x2) < 0){ # temp_x1 and temp_x2 are negatively correlated
          corM3[i,j] = -corM3[i,j]
          corM3[j,i] = corM3[i,j]
        }
        
        rho3[k] = corM3[i,j]
        
        if (covM){
          MSOCdata[[k]] = list(x1 = temp_x1,
                               x2 = temp_x2)
        }
        
        k = k+1
      }
      print(i)
    }
    temp_t = Sys.time() - temp_t
    
    if (covM){
      # Computing the cov matrix of the CC, MOC and MSOC for each study
      temp_t2 = Sys.time()
      M1 = diag((1-rho1 ^2)^2 / n[s]) #For CC
      M2 = diag((1-rho2 ^2)^2 / n[s]) #For MOC
      M3 = diag((1-rho3 ^2)^2 / n[s]) #For MSOC

      for (i in (1:(pstar-1))){
        for (j in (i+1):pstar){
          M1[i,j] = CCcov(CCdata[[i]]$x1,CCdata[[i]]$x2,
                          CCdata[[j]]$x1,CCdata[[j]]$x2,n[s])
          M1[j,i] = M1[i,j]

          M2[i,j] = CCcov(MOCdata[[i]]$x1,MOCdata[[i]]$x2,
                          MOCdata[[j]]$x1,MOCdata[[j]]$x2,n[s])
          M2[j,i] = M2[i,j]

          M3[i,j] = CCcov(MSOCdata[[i]]$x1,MSOCdata[[i]]$x2,
                          MSOCdata[[j]]$x1,MSOCdata[[j]]$x2,n[s])
          M3[j,i] = M3[i,j]
        }
        # print(i)
      }
      temp_t2 = Sys.time() - temp_t2
    }
    
    print(s)
    
    if (covM){
      Studies[[s]] = list(omic = omic,
                          CCdata = CCdata,
                          MOCdata = MOCdata,
                          MSOCdata = MSOCdata,
                          rho1 = rho1,
                          rho2 = rho2,
                          rho3 = rho3,
                          cancorM = corM,
                          Max_Omic_Cor = corM2,
                          Max_Same_Omic_Cor = corM3,
                          rho1_cov = M1,
                          rho2_cov = M2,
                          rho3_cov = M3)
    } else {
      Studies[[s]] = list(omic = omic,
                          rho1 = rho1,#
                          rho2 = rho2,
                          rho3 = rho3,
                          cancorM = corM,
                          Max_Omic_Cor = corM2,
                          Max_Same_Omic_Cor = corM3)
    }
    
  }
  
  return(Studies)
}


Combined_rho = function(
    list, #The output of the Studies_summary() function
    G, #The number of genes consider in each study, we assume it is the same foe all the studies
    diag_cov=TRUE, #Whether we are computing the the methods using the diagonal of the correlation covariance matrix or not
    covM=TRUE #Whether, the correlation covariance matrices are available in list.
    #If they are available, it will also compute the methods using the whole correlation covariance matrices 
    ){ 
  S = length(list) #the # of studies
  pstar = G*(G-1)/2
  Sum_Wr1 = rep(0,pstar)
  Sum_Wr2 = rep(0,pstar)
  Sum_Wr3 = rep(0,pstar)
  if (diag_cov){ #We want to consider the diagonal of the correlations covariance matrix
    #We only keep the diagonal of the matrix as a vector
    SumW1 = rep(0,pstar)
    SumW2 = rep(0,pstar)
    SumW3 = rep(0,pstar)
    for (s in 1:S){
      if (covM){
        W1 = round(1/diag(list[[s]]$rho1_cov),digits=5)
        W2 = round(1/diag(list[[s]]$rho2_cov),digits=5)
        W3 = round(1/diag(list[[s]]$rho3_cov),digits=5)
      } else { #We don't have the covariance matrix
        W1 = round(n[s] / (1-list[[s]]$rho1 ^2)^2,digits=5)
        W2 = round(n[s] / (1-list[[s]]$rho2 ^2)^2,digits=5)
        W3 = round(n[s] / (1-list[[s]]$rho3 ^2)^2,digits=5)
      }
      SumW1 = SumW1 + W1
      SumW2 = SumW2 + W2
      SumW3 = SumW3 + W3
      Sum_Wr1 = Sum_Wr1 + W1*list[[s]]$rho1
      Sum_Wr2 = Sum_Wr2 + W2*list[[s]]$rho2
      Sum_Wr3 = Sum_Wr3 + W3*list[[s]]$rho3
      
      print(s)
    }
    comb_rho1 = round(1/SumW1,digits=5)*Sum_Wr1
    comb_rho2 = round(1/SumW2,digits=5)*Sum_Wr2
    comb_rho3 = round(1/SumW3,digits=5)*Sum_Wr3
    
  } else { #We want to consider the whole correlations covariance matrix
    SumW1 = matrix(rep(0,pstar^2), nrow = pstar)
    SumW2 = matrix(rep(0,pstar^2), nrow = pstar)
    SumW3 = matrix(rep(0,pstar^2), nrow = pstar)
    for (s in 1:S){
      #We are rounding the inverses to insure we numerically get positive definite matrix as output
      #This will insure we theoretically get the combined correlation's value in [0,1]
      W1 = round(solve(list[[s]]$rho1_cov, tol = 10^(-20)),digits=5)
      W2 = round(solve(list[[s]]$rho2_cov, tol = 10^(-20)),digits=5)
      W3 = round(solve(list[[s]]$rho3_cov, tol = 10^(-20)),digits=5)
      SumW1 = SumW1 + W1
      SumW2 = SumW2 + W2
      SumW3 = SumW3 + W3
      Sum_Wr1 = Sum_Wr1 + W1%*%list[[s]]$rho1
      Sum_Wr2 = Sum_Wr2 + W2%*%list[[s]]$rho2
      Sum_Wr3 = Sum_Wr3 + W3%*%list[[s]]$rho3
      
      print(s)
    }
    comb_rho1 = round(solve(SumW1, tol = 10^(-20)),digits=5)%*%Sum_Wr1
    comb_rho2 = round(solve(SumW2, tol = 10^(-20)),digits=5)%*%Sum_Wr2
    comb_rho3 = round(solve(SumW3, tol = 10^(-20)),digits=5)%*%Sum_Wr3
  }
  
  #Making the combiened CC, MOC and MSOC matrices#
  comb_CC = diag(rep(1,G))
  comb_MOC = diag(rep(1,G))
  comb_MSOC = diag(rep(1,G))
  k = 1
  for (i in 1:(G-1)){
    for (j in (i+1):G){
      comb_CC[i,j] = comb_rho1[k]
      comb_CC[j,i] = comb_rho1[k]
      comb_MOC[i,j] = comb_rho2[k]
      comb_MOC[j,i] = comb_rho2[k]
      comb_MSOC[i,j] = comb_rho3[k]
      comb_MSOC[j,i] = comb_rho3[k]
      k = k+1
    }
    #print(i)
  }
  
  list(comb_rho1 = comb_rho1,
       comb_CC = comb_CC,
       comb_rho2 = comb_rho2,
       comb_MOC = comb_MOC,
       comb_rho3 = comb_rho3,
       comb_MSOC = comb_MSOC)
} 


##Function that selects the power based on the second derivative to the mean connectivity#
pickPower = function(x){#x is the change rate of the mean connectivity change rate
  trend = x[1] < x[2]
  new_trend = trend
  i=2
  while (trend == new_trend & i < length(x)){
    new_trend = x[i] < x[i+1]
    if (trend == new_trend){
      i = i+1
    }
  }
  if(trend != new_trend){
    output = i+2
  } else{
    if(trend==T){ #the curve is increasing
      output = length(x)+2
    } else { #the curve is decreasing
      output = 3
    }
  }
  output
}



###############################################################################################
#####################CODE FOR RUNNING THE REAL DATA############################################
###############################################################################################
powers = 1:40
minClusterSize=5

###Clustering using CC###

sft = pickSoftThreshold.fromSimilarity(similarity=as.matrix(CC_M), powerVector = powers, verbose = 0)
connectivity = sft$fitIndices$mean.k.
k = length(connectivity)
percentChange = - (connectivity[-1]-connectivity[-k])/connectivity[-k]
percentChange2 = - (percentChange[-1]-percentChange[-(k-1)])/percentChange[-(k-1)]

softPower = pickPower(percentChange2)

# par(mfrow = c(1,2));
# cex1 = 0.9;
# # Scale-free topology fit index as a function of the soft-thresholding power
# plot(sft$fitIndices[,1], sft$fitIndices[,4], #Why -sign of the slope?
#      xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed truncated R^2",type="n",
#      main = paste("Scale independence"));
# text(sft$fitIndices[,1],sft$fitIndices[,4],labels=powers,cex=cex1,col="red");
# # this line corresponds to using an R^2 cut-off of h
# abline(h=0.6,col="red")
# # Mean connectivity as a function of the soft-thresholding power
# plot(sft$fitIndices[,1], sft$fitIndices[,5],
#      xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
#      main = paste("Mean connectivity"))
# text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
# 
# #plot(powers[-1], -(sft$fitIndices$mean.k.[-1]-sft$fitIndices$mean.k.[-39])/sft$fitIndices$mean.k.[-39], type = "l",ylab = "Mean connectivity relative change", xlab = "Powers")
# plot(powers, connectivity, type = "l",lwd=2,ylab = "Mean connectivity", xlab = "Powers")
# abline(v=softPower,col="red",lwd=2)
# plot(powers[-1], percentChange, type = "l",lwd=2,ylab = "Mean connectivity percentage change", xlab = "Powers")
# abline(v=softPower,col="red",lwd=2)
# plot(powers[-(1:2)], percentChange2, type = "l",lwd=2,ylab = "Percentage change of Mean connectivity percentage change", xlab = "Powers")
# abline(v=softPower,col="red")

adjacency = CC_M^softPower
TOM = TOMsimilarity(as.matrix(adjacency))
dissTOM = 1-TOM

geneTree = hclust(as.dist(dissTOM), method = "average")
minModuleSize = minClusterSize; #minimum number of clusters
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)

###Considering each unassigned gene as a cluster####
names(dynamicMods) = NULL
unclass_indexes = dynamicMods == 0
n_unclass = sum(unclass_indexes)
if (n_unclass != 0){
  dynamicMods[unclass_indexes] = max(dynamicMods) + (1:n_unclass)
}

CC_True_Class = dynamicMods


###Clustering using MOC###

sft = pickSoftThreshold.fromSimilarity(similarity=as.matrix(abs(MOC_M)), powerVector = powers, verbose = 0)
connectivity = sft$fitIndices$mean.k.
k = length(connectivity)
percentChange = - (connectivity[-1]-connectivity[-k])/connectivity[-k]
percentChange2 = - (percentChange[-1]-percentChange[-(k-1)])/percentChange[-(k-1)]

softPower = pickPower(percentChange2)

# par(mfrow = c(1,2));
# cex1 = 0.9;
# # Scale-free topology fit index as a function of the soft-thresholding power
# plot(sft$fitIndices[,1], sft$fitIndices[,4], #Why -sign of the slope?
#      xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed truncated R^2",type="n",
#      main = paste("Scale independence"));
# text(sft$fitIndices[,1],sft$fitIndices[,4],labels=powers,cex=cex1,col="red");
# # this line corresponds to using an R^2 cut-off of h
# abline(h=0.6,col="red")
# # Mean connectivity as a function of the soft-thresholding power
# plot(sft$fitIndices[,1], sft$fitIndices[,5],
#      xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
#      main = paste("Mean connectivity"))
# text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
# 
# #plot(powers[-1], -(sft$fitIndices$mean.k.[-1]-sft$fitIndices$mean.k.[-39])/sft$fitIndices$mean.k.[-39], type = "l",ylab = "Mean connectivity relative change", xlab = "Powers")
# plot(powers, connectivity, type = "l",lwd=2,ylab = "Mean connectivity", xlab = "Powers")
# abline(v=softPower,col="red",lwd=2)
# plot(powers[-1], percentChange, type = "l",lwd=2,ylab = "Mean connectivity percentage change", xlab = "Powers")
# abline(v=softPower,col="red",lwd=2)
# plot(powers[-(1:2)], percentChange2, type = "l",lwd=2,ylab = "Percentage change of Mean connectivity percentage change", xlab = "Powers")
# abline(v=softPower,col="red")

adjacency = abs(MOC_M)^softPower
TOM = TOMsimilarity(as.matrix(adjacency))
dissTOM = 1-TOM

geneTree = hclust(as.dist(dissTOM), method = "average")
minModuleSize = minClusterSize; #minimum number of clusters
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)

###Considering each unassigned gene as a cluster####
names(dynamicMods) = NULL
unclass_indexes = dynamicMods == 0
n_unclass = sum(unclass_indexes)
if (n_unclass != 0){
  dynamicMods[unclass_indexes] = max(dynamicMods) + (1:n_unclass)
}

MOC_True_Class = dynamicMods


###Clustering using MSOC###

sft = pickSoftThreshold.fromSimilarity(similarity=as.matrix(abs(MSOC_M)), powerVector = powers, verbose = 0)
connectivity = sft$fitIndices$mean.k.
k = length(connectivity)
percentChange = - (connectivity[-1]-connectivity[-k])/connectivity[-k]
percentChange2 = - (percentChange[-1]-percentChange[-(k-1)])/percentChange[-(k-1)]

softPower = pickPower(percentChange2)

# par(mfrow = c(1,2));
# cex1 = 0.9;
# # Scale-free topology fit index as a function of the soft-thresholding power
# plot(sft$fitIndices[,1], sft$fitIndices[,4], #Why -sign of the slope?
#      xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed truncated R^2",type="n",
#      main = paste("Scale independence"));
# text(sft$fitIndices[,1],sft$fitIndices[,4],labels=powers,cex=cex1,col="red");
# # this line corresponds to using an R^2 cut-off of h
# abline(h=0.6,col="red")
# # Mean connectivity as a function of the soft-thresholding power
# plot(sft$fitIndices[,1], sft$fitIndices[,5],
#      xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
#      main = paste("Mean connectivity"))
# text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
# 
# #plot(powers[-1], -(sft$fitIndices$mean.k.[-1]-sft$fitIndices$mean.k.[-39])/sft$fitIndices$mean.k.[-39], type = "l",ylab = "Mean connectivity relative change", xlab = "Powers")
# plot(powers, connectivity, type = "l",lwd=2,ylab = "Mean connectivity", xlab = "Powers")
# abline(v=softPower,col="red",lwd=2)
# plot(powers[-1], percentChange, type = "l",lwd=2,ylab = "Mean connectivity percentage change", xlab = "Powers")
# abline(v=softPower,col="red",lwd=2)
# plot(powers[-(1:2)], percentChange2, type = "l",lwd=2,ylab = "Percentage change of Mean connectivity percentage change", xlab = "Powers")
# abline(v=softPower,col="red")

adjacency = abs(MSOC_M)^softPower
TOM = TOMsimilarity(as.matrix(adjacency))
dissTOM = 1-TOM

geneTree = hclust(as.dist(dissTOM), method = "average")
minModuleSize = minClusterSize; #minimum number of clusters
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)

###Considering each unassigned gene as a cluster####
names(dynamicMods) = NULL
unclass_indexes = dynamicMods == 0
n_unclass = sum(unclass_indexes)
if (n_unclass != 0){
  dynamicMods[unclass_indexes] = max(dynamicMods) + (1:n_unclass)
}

MSOC_True_Class = dynamicMods

adjustedRandIndex(CC_True_Class,MOC_True_Class)
adjustedRandIndex(CC_True_Class,MSOC_True_Class)
adjustedRandIndex(MSOC_True_Class,MOC_True_Class)

# save.image(paste0("/scratch/u/ul55964/Real_data_analysis_Aim2/Aim2Meta_true_data_",G,"_genes.RData"))





##Chopping the data into 5 to get the 5 studies##
# Original_omic = omic
# set.seed = 12345
# ss = sample(1:5, size = nrow(Original_omic),replace = T, prob = rep(1/5,5))
# save.image(paste0("/scratch/u/ul55964/Real_data_analysis_Aim2/Aim2Meta_true_data.RData"))

##Chopping the data into 3 to get the 3 studies##
# load("/scratch/u/ul55964/Real_data_analysis_Aim2/Aim2Meta_true_data.RData")
Original_omic = omic
set.seed = 12345
ss = sample(1:3, size = nrow(Original_omic),replace = T, prob = rep(1/3,3))
S = 3
save.image(paste0("/scratch/u/ul55964/Real_data_analysis_Aim2/Aim2Meta_true_data_",S,"Studies_",G,"_genes.RData"))

##Chopping the data into 2 to get the 2 studies##
# load("/scratch/u/ul55964/Real_data_analysis_Aim2/Aim2Meta_true_data.RData")
# # Original_omic = omic
# set.seed = 12345
# ss = sample(1:2, size = nrow(Original_omic),replace = T, prob = rep(1/2,2))
# S = 2
# save.image(paste0("/scratch/u/ul55964/Real_data_analysis_Aim2/Aim2Meta_true_data_",S,"Studies.RData"))


G=123
S=3
# load(paste0("/scratch/u/ul55964/Real_data_analysis_Aim2/Aim2Meta_true_data_",S,"Studies.RData"))
load(paste0("/scratch/u/ul55964/Real_data_analysis_Aim2/Aim2Meta_true_data_",S,"Studies_",G,"_genes.RData"))
n1 = rep(NA,S)
Studies = list()

for (s in 1:S){
  # load(paste0("/scratch/u/ul55964/Real_data_analysis_Aim2/Aim2Meta_Study_",s,".RData"))
  # load(paste0("/scratch/u/ul55964/Real_data_analysis_Aim2/Aim2Meta_",S,"_Studies_Study_",s,".RData"))
  # load(paste0("/scratch/u/ul55964/Real_data_analysis_Aim2/Aim2Meta_",S,"_Studies_Study_",s,"_",G,"_genes.RData"))
  load(paste0("/scratch/u/ul55964/Real_data_analysis_Aim2/Aim2Meta_",S,"_Studies_Study_",s,"_",G,"_genes_All.RData"))
  Studies[[s]]=Study
  n1[s] = n
}

n= n1
rm(n1)

pstar = G*(G-1)/2
##For diagonal Methods##
Sum_Wr1 = rep(0,pstar)
Sum_Wr2 = rep(0,pstar)
Sum_Wr3 = rep(0,pstar)

SumW1 = rep(0,pstar)
SumW2 = rep(0,pstar)
SumW3 = rep(0,pstar)

for (s in 1:S){
  W1 = round(n[s] / (1-Studies[[s]]$rho1 ^2)^2,digits=5)
  W2 = round(n[s] / (1-Studies[[s]]$rho2 ^2)^2,digits=5)
  W3 = round(n[s] / (1-Studies[[s]]$rho3 ^2)^2,digits=5)
  SumW1 = SumW1 + W1
  SumW2 = SumW2 + W2
  SumW3 = SumW3 + W3
  Sum_Wr1 = Sum_Wr1 + W1*Studies[[s]]$rho1
  Sum_Wr2 = Sum_Wr2 + W2*Studies[[s]]$rho2
  Sum_Wr3 = Sum_Wr3 + W3*Studies[[s]]$rho3
  
  print(s)
}

comb_rho1 = round(1/SumW1,digits=5)*Sum_Wr1
comb_rho2 = round(1/SumW2,digits=5)*Sum_Wr2
comb_rho3 = round(1/SumW3,digits=5)*Sum_Wr3

#Some for epsilon machine are greater than 1
# Old_comb_rho1 = comb_rho1
# Old_comb_rho2 = comb_rho2
# Old_comb_rho3 = comb_rho3
# 
# #Making values greater than 1, one
# comb_rho1[comb_rho1 > 1] = 1
# comb_rho2[comb_rho2 > 1] = 1
# comb_rho3[comb_rho3 > 1] = 1


#Making the combiened CC, MOC and MSOC matrices#
comb_CC = diag(rep(1,G))
comb_MOC = diag(rep(1,G))
comb_MSOC = diag(rep(1,G))
k = 1
for (i in 1:(G-1)){
  for (j in (i+1):G){
    comb_CC[i,j] = comb_rho1[k]
    comb_CC[j,i] = comb_rho1[k]
    comb_MOC[i,j] = comb_rho2[k]
    comb_MOC[j,i] = comb_rho2[k]
    comb_MSOC[i,j] = comb_rho3[k]
    comb_MSOC[j,i] = comb_rho3[k]
    k = k+1
  }
  print(i)
}

# Comb_Cor = list(comb_rho1 = comb_rho1,
#                 comb_CC = comb_CC,
#                 comb_rho2 = comb_rho2,
#                 comb_MOC = comb_MOC,
#                 comb_rho3 = comb_rho3,
#                 comb_MSOC = comb_MSOC)



ARI = matrix(nrow=3, ncol=10+3*S)
names(ARI) = c(paste0(c("Meta_CC","Meta_MOC","Meta_MSOC"),"_Diag"),
               paste0("CC-",1:S),paste0("MOC-",1:S),paste0("MSOC-",1:S),"WGCNA(1)-1","WGCNA(2)-1",
               "K-M(1)-1", "K-M(2)-1","Meta_CC","Meta_MOC","Meta_MSOC")

rownames(ARI) = c("CC_True_Class","MOC_True_Class","MSOC_True_Class")


# For diagonal methods only#
# ARI = matrix(nrow=3, ncol=7+3*S)
# colnames(ARI) = c(paste0(c("Meta_CC","Meta_MOC","Meta_MSOC"),"_Diag"),
#                   paste0("CC-",1:S),paste0("MOC-",1:S),paste0("MSOC-",1:S),"WGCNA(1)-1","WGCNA(2)-1",
#                   "K-M(1)-1", "K-M(2)-1")
# rownames(ARI) = c("CC_True_Class","MOC_True_Class","MSOC_True_Class")
Clust_results = list()


##The Meta-CC method##
sft = pickSoftThreshold.fromSimilarity(similarity=abs(comb_CC), powerVector = powers, verbose = 0)

connectivity = sft$fitIndices$mean.k.
k = length(connectivity)
percentChange = - (connectivity[-1]-connectivity[-k])/connectivity[-k]
percentChange2 = - (percentChange[-1]-percentChange[-(k-1)])/percentChange[-(k-1)]

softPower = pickPower(percentChange2)
adjacency = abs(comb_CC)^softPower
TOM = TOMsimilarity(as.matrix(adjacency))
dissTOM = 1-TOM

geneTree = hclust(as.dist(dissTOM), method = "average")
minModuleSize = minClusterSize; #minimum number of clusters
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)

###Considering each unassigned gene as a cluster####
names(dynamicMods) = NULL
unclass_indexes = dynamicMods == 0
n_unclass = sum(unclass_indexes)
if (n_unclass != 0){
  dynamicMods[unclass_indexes] = max(dynamicMods) + (1:n_unclass)
}

Clust_results[[1]] = dynamicMods
ARI["CC_True_Class",1] = adjustedRandIndex(CC_True_Class,dynamicMods)
ARI["MOC_True_Class",1] = adjustedRandIndex(MOC_True_Class,dynamicMods)
ARI["MSOC_True_Class",1] = adjustedRandIndex(MSOC_True_Class,dynamicMods)

##The Meta-MOC method##
rm(sft,connectivity,k,percentChange,percentChange2,softPower,adjacency,TOM,dissTOM,geneTree,dynamicMods)
sft = pickSoftThreshold.fromSimilarity(similarity=abs(comb_MOC), powerVector = powers, verbose = 0)

connectivity = sft$fitIndices$mean.k.
k = length(connectivity)
percentChange = - (connectivity[-1]-connectivity[-k])/connectivity[-k]
percentChange2 = - (percentChange[-1]-percentChange[-(k-1)])/percentChange[-(k-1)]

softPower = pickPower(percentChange2)
adjacency = abs(comb_MOC)^softPower
TOM = TOMsimilarity(as.matrix(adjacency))
dissTOM = 1-TOM

geneTree = hclust(as.dist(dissTOM), method = "average")
minModuleSize = minClusterSize; #minimum number of clusters
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)

###Considering each unassigned gene as a cluster####
names(dynamicMods) = NULL
unclass_indexes = dynamicMods == 0
n_unclass = sum(unclass_indexes)
if (n_unclass != 0){
  dynamicMods[unclass_indexes] = max(dynamicMods) + (1:n_unclass)
}

Clust_results[[2]] = dynamicMods
ARI["CC_True_Class",2] = adjustedRandIndex(CC_True_Class,dynamicMods)
ARI["MOC_True_Class",2] = adjustedRandIndex(MOC_True_Class,dynamicMods)
ARI["MSOC_True_Class",2] = adjustedRandIndex(MSOC_True_Class,dynamicMods)

##The Meta-MSOC method##
rm(sft,connectivity,k,percentChange,percentChange2,softPower,adjacency,TOM,dissTOM,geneTree,dynamicMods)
sft = pickSoftThreshold.fromSimilarity(similarity=abs(comb_MSOC), powerVector = powers, verbose = 0)

connectivity = sft$fitIndices$mean.k.
k = length(connectivity)
percentChange = - (connectivity[-1]-connectivity[-k])/connectivity[-k]
percentChange2 = - (percentChange[-1]-percentChange[-(k-1)])/percentChange[-(k-1)]

softPower = pickPower(percentChange2)
adjacency = abs(comb_MSOC)^softPower
TOM = TOMsimilarity(as.matrix(adjacency))
dissTOM = 1-TOM

geneTree = hclust(as.dist(dissTOM), method = "average")
minModuleSize = minClusterSize; #minimum number of clusters
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)

###Considering each unassigned gene as a cluster####
names(dynamicMods) = NULL
unclass_indexes = dynamicMods == 0
n_unclass = sum(unclass_indexes)
if (n_unclass != 0){
  dynamicMods[unclass_indexes] = max(dynamicMods) + (1:n_unclass)
}

Clust_results[[3]] = dynamicMods
ARI["CC_True_Class",3] = adjustedRandIndex(CC_True_Class,dynamicMods)
ARI["MOC_True_Class",3] = adjustedRandIndex(MOC_True_Class,dynamicMods)
ARI["MSOC_True_Class",3] = adjustedRandIndex(MSOC_True_Class,dynamicMods)


##CC for each study##
rm(sft,connectivity,k,percentChange,percentChange2,softPower,adjacency,TOM,dissTOM,geneTree,dynamicMods)
for (s in 1:S){
  corM = Studies[[s]]$cancorM #Considering the sth study
  sft = pickSoftThreshold.fromSimilarity(similarity=as.matrix(corM), powerVector = powers, verbose = 0)
  
  connectivity = sft$fitIndices$mean.k.
  k = length(connectivity)
  percentChange = - (connectivity[-1]-connectivity[-k])/connectivity[-k]
  percentChange2 = - (percentChange[-1]-percentChange[-(k-1)])/percentChange[-(k-1)]
  
  softPower = pickPower(percentChange2)
  adjacency = corM^softPower
  TOM = TOMsimilarity(as.matrix(adjacency))
  dissTOM = 1-TOM
  
  geneTree = hclust(as.dist(dissTOM), method = "average")
  minModuleSize = minClusterSize; #minimum number of clusters
  # Module identification using dynamic tree cut:
  dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                              deepSplit = 2, pamRespectsDendro = FALSE,
                              minClusterSize = minModuleSize)
  
  ###Considering each unassigned gene as a cluster####
  names(dynamicMods) = NULL
  unclass_indexes = dynamicMods == 0
  n_unclass = sum(unclass_indexes)
  if (n_unclass != 0){
    dynamicMods[unclass_indexes] = max(dynamicMods) + (1:n_unclass)
  }
  
  Clust_results[[3+s]] = dynamicMods
  ARI["CC_True_Class",3+s] = adjustedRandIndex(CC_True_Class,dynamicMods)
  ARI["MOC_True_Class",3+s] = adjustedRandIndex(MOC_True_Class,dynamicMods)
  ARI["MSOC_True_Class",3+s] = adjustedRandIndex(MSOC_True_Class,dynamicMods)
  
  print(paste0(s,"th study CC"))
  rm(corM,sft,connectivity,k,percentChange,percentChange2,softPower,adjacency,TOM,dissTOM,geneTree,dynamicMods)
}

##Maximum omic correlation for each study##

for (s in 1:S){
  corM2 = Studies[[s]]$Max_Omic_Cor
  sft = pickSoftThreshold.fromSimilarity(similarity=as.matrix(abs(corM2)), powerVector = powers, verbose = 0)
  
  connectivity = sft$fitIndices$mean.k.
  k = length(connectivity)
  percentChange = - (connectivity[-1]-connectivity[-k])/connectivity[-k]
  percentChange2 = - (percentChange[-1]-percentChange[-(k-1)])/percentChange[-(k-1)]
  
  softPower = pickPower(percentChange2)
  adjacency = abs(corM2)^softPower
  TOM = TOMsimilarity(as.matrix(adjacency))
  dissTOM = 1-TOM
  
  geneTree = hclust(as.dist(dissTOM), method = "average")
  minModuleSize = minClusterSize; #minimum number of clusters
  # Module identification using dynamic tree cut:
  dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                              deepSplit = 2, pamRespectsDendro = FALSE,
                              minClusterSize = minModuleSize)
  
  ###Considering each unassigned gene as a cluster####
  names(dynamicMods) = NULL
  unclass_indexes = dynamicMods == 0
  n_unclass = sum(unclass_indexes)
  if (n_unclass != 0){
    dynamicMods[unclass_indexes] = max(dynamicMods) + (1:n_unclass)
  }
  
  Clust_results[[3+S+s]] = dynamicMods
  ARI["CC_True_Class",3+S+s] = adjustedRandIndex(CC_True_Class,dynamicMods)
  ARI["MOC_True_Class",3+S+s] = adjustedRandIndex(MOC_True_Class,dynamicMods)
  ARI["MSOC_True_Class",3+S+s] = adjustedRandIndex(MSOC_True_Class,dynamicMods)
  
  print(paste0(s,"th study MOC"))
  rm(sft,connectivity,k,percentChange,percentChange2,softPower,adjacency,TOM,dissTOM,geneTree,dynamicMods)
}

##Maximum Same Omic correlation for each study##
for (s in 1:S){
  corM3 = Studies[[s]]$Max_Same_Omic_Cor
  sft = pickSoftThreshold.fromSimilarity(similarity=as.matrix(abs(corM3)), powerVector = powers, verbose = 0)
  
  connectivity = sft$fitIndices$mean.k.
  k = length(connectivity)
  percentChange = - (connectivity[-1]-connectivity[-k])/connectivity[-k]
  percentChange2 = - (percentChange[-1]-percentChange[-(k-1)])/percentChange[-(k-1)]
  
  softPower = pickPower(percentChange2)
  adjacency = abs(corM3)^softPower
  TOM = TOMsimilarity(as.matrix(adjacency))
  dissTOM = 1-TOM
  
  geneTree = hclust(as.dist(dissTOM), method = "average")
  minModuleSize = minClusterSize; #minimum number of clusters
  # Module identification using dynamic tree cut:
  dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                              deepSplit = 2, pamRespectsDendro = FALSE,
                              minClusterSize = minModuleSize)
  
  ###Considering each unassigned gene as a cluster####
  names(dynamicMods) = NULL
  unclass_indexes = dynamicMods == 0
  n_unclass = sum(unclass_indexes)
  if (n_unclass != 0){
    dynamicMods[unclass_indexes] = max(dynamicMods) + (1:n_unclass)
  }
  
  Clust_results[[3+2*S+s]] = dynamicMods
  ARI["CC_True_Class",3+2*S+s] = adjustedRandIndex(CC_True_Class,dynamicMods)
  ARI["MOC_True_Class",3+2*S+s] = adjustedRandIndex(MOC_True_Class,dynamicMods)
  ARI["MSOC_True_Class",3+2*S+s] = adjustedRandIndex(MSOC_True_Class,dynamicMods)
  
  print(paste0(s,"th study MSOC"))
  rm(sft,connectivity,k,percentChange,percentChange2,softPower,adjacency,TOM,dissTOM,geneTree,dynamicMods)
}

###One-omic methods### with omic 1 only of the 1st study
single_omic1 = Studies[[1]]$omic[,seq(1,total_individuals,numb_omic)]
M1 = abs(cor(single_omic1)) #Correlation matrix from omic 1
# single_omic2 = omic[,seq(2,3000,3)]
# single_omic3 = omic[,seq(3,3000,3)]


#The first omic#
sft = pickSoftThreshold.fromSimilarity(similarity=M1, powerVector = powers, verbose = 0)

#One-omic with the new cut function#
connectivity = sft$fitIndices$mean.k.
k = length(connectivity)
percentChange = - (connectivity[-1]-connectivity[-k])/connectivity[-k]
percentChange2 = - (percentChange[-1]-percentChange[-(k-1)])/percentChange[-(k-1)]
softPower = pickPower(percentChange2)
adjacency = M1^softPower
TOM = TOMsimilarity(as.matrix(adjacency))
dissTOM = 1-TOM
geneTree = hclust(as.dist(dissTOM), method = "average")
minModuleSize = minClusterSize; #minimum number of clusters
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)

###Considering each unassigned gene as a cluster####
names(dynamicMods) = NULL
unclass_indexes = dynamicMods == 0
n_unclass = sum(unclass_indexes)
if (n_unclass != 0){
  dynamicMods[unclass_indexes] = max(dynamicMods) + (1:n_unclass)
}

Clust_results[[3+3*S+1]] = dynamicMods
ARI["CC_True_Class",3+3*S+1] = adjustedRandIndex(CC_True_Class,dynamicMods)
ARI["MOC_True_Class",3+3*S+1] = adjustedRandIndex(MOC_True_Class,dynamicMods)
ARI["MSOC_True_Class",3+3*S+1] = adjustedRandIndex(MSOC_True_Class,dynamicMods)

##One-omic traditional WGCNA##
rm(softPower,adjacency,TOM,dissTOM,geneTree,dynamicMods)

PowerCut3 = function(powers,trunc.R.sq,connectivity,min.trunc.R.sq,changeCut){
  n = length(connectivity)
  percentChange = - (connectivity[-1]-connectivity[-n])/connectivity[-n]
  interest = trunc.R.sq >= min.trunc.R.sq
  if (sum(interest) == 0){
    selectedpower = powers[which.max(trunc.R.sq)]
  } else {
    potentialpowers = powers[interest]
    for (i in potentialpowers){
      index = which(powers == i)
      if (index[1] != n){
        if (percentChange[index[1]] < changeCut){
          selectedpower = i
          break
        }
        
        if (index[1] != (n-1)){#This might not be necessary, since the percent change is always decreasing when the powers are evenly scaled
          if (percentChange[index[1]] < percentChange[index[1]+1]){
            selectedpower = i
            break
          }
        }
      }
    }
    selectedpower = i
  }
  return(selectedpower)
}
softPower = PowerCut3(sft$fitIndices$Power,sft$fitIndices$truncated.R.sq,sft$fitIndices$mean.k.,0.6,1/3)
adjacency = M1^softPower
TOM = TOMsimilarity(as.matrix(adjacency))
dissTOM = 1-TOM

geneTree = hclust(as.dist(dissTOM), method = "average")
minModuleSize = minClusterSize; #minimum number of clusters
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)

###Considering each unassigned gene as a cluster####
names(dynamicMods) = NULL
unclass_indexes = dynamicMods == 0
n_unclass = sum(unclass_indexes)
if (n_unclass != 0){
  dynamicMods[unclass_indexes] = max(dynamicMods) + (1:n_unclass)
}

Clust_results[[3+3*S+2]] = dynamicMods
ARI["CC_True_Class",3+3*S+2] = adjustedRandIndex(CC_True_Class,dynamicMods)
ARI["MOC_True_Class",3+3*S+2] = adjustedRandIndex(MOC_True_Class,dynamicMods)
ARI["MSOC_True_Class",3+3*S+2] = adjustedRandIndex(MSOC_True_Class,dynamicMods)

##K-means with the true k##
set.seed(12345)
Kmean_result = kmeans(t(scale(single_omic1)),centers = length(table(CC_True_Class)),nstart = 25)
Clust_results[[3+3*S+3]] = Kmean_result$cluster
ARI["CC_True_Class",3+3*S+3] = adjustedRandIndex(CC_True_Class,Kmean_result$cluster)

set.seed(12345)
Kmean_result = kmeans(t(scale(single_omic1)),centers = length(table(MOC_True_Class)),nstart = 25)
Clust_results[[3+3*S+4]] = Kmean_result$cluster
ARI["MOC_True_Class",3+3*S+3] = adjustedRandIndex(MOC_True_Class,Kmean_result$cluster)

set.seed(12345)
Kmean_result = kmeans(t(scale(single_omic1)),centers = length(table(MSOC_True_Class)),iter.max = 20,nstart = 25)
Clust_results[[3+3*S+5]] = Kmean_result$cluster
ARI["MSOC_True_Class",3+3*S+3] = adjustedRandIndex(MSOC_True_Class,Kmean_result$cluster)

##K-means with the gap statistics for choosing the optimal k between 1:130##
set.seed(12345)
gaps = clusGap(t(scale(single_omic1)), kmeans, iter.max = 20,
               K.max=max(length(table(MSOC_True_Class)),length(table(MOC_True_Class)),length(table(CC_True_Class)))+10, 
               B = 100, d.power = 1)
rm(Kmean_result)
opt_k = maxSE(gaps$Tab[,"gap"], gaps$Tab[,"SE.sim"],method = "globalSEmax")
set.seed(12345)
Kmean_result = kmeans(t(scale(single_omic1)),centers = opt_k,iter.max = 20,nstart = 25)

Clust_results[[3+3*S+6]] = Kmean_result$cluster
ARI["CC_True_Class",3+3*S+4] = adjustedRandIndex(CC_True_Class,Kmean_result$cluster)
ARI["MOC_True_Class",3+3*S+4] = adjustedRandIndex(MOC_True_Class,Kmean_result$cluster)
ARI["MSOC_True_Class",3+3*S+4] = adjustedRandIndex(MSOC_True_Class,Kmean_result$cluster)


####Whole covariance matrix methods#####
rm(comb_rho1,comb_rho2,comb_rho3,Sum_Wr1,Sum_Wr2,Sum_Wr3,SumW1,SumW2,SumW3,sft,connectivity,k,percentChange,percentChange2,softPower,adjacency,TOM,dissTOM,geneTree,dynamicMods)
Sum_Wr1 = rep(0,pstar)
Sum_Wr2 = rep(0,pstar)
Sum_Wr3 = rep(0,pstar)
SumW1 = matrix(rep(0,pstar^2), nrow = pstar)
SumW2 = matrix(rep(0,pstar^2), nrow = pstar)
SumW3 = matrix(rep(0,pstar^2), nrow = pstar)
for (s in 1:S){
  #We are rounding the inverses to insure we numerically get positive definite matrix as output
  #This will insure we theoretically get the combined correlation's value in [0,1]
  W1 = round(solve(Studies[[s]]$rho1_cov, tol = 10^(-25)),digits=5)
  W2 = round(solve(Studies[[s]]$rho2_cov, tol = 10^(-25)),digits=5)
  W3 = round(solve(Studies[[s]]$rho3_cov, tol = 10^(-25)),digits=5)
  SumW1 = SumW1 + W1
  SumW2 = SumW2 + W2
  SumW3 = SumW3 + W3
  Sum_Wr1 = Sum_Wr1 + W1%*%Studies[[s]]$rho1
  Sum_Wr2 = Sum_Wr2 + W2%*%Studies[[s]]$rho2
  Sum_Wr3 = Sum_Wr3 + W3%*%Studies[[s]]$rho3
  
  print(s)
}
comb_rho1 = round(solve(SumW1, tol = 10^(-25)),digits=5)%*%Sum_Wr1
comb_rho2 = round(solve(SumW2, tol = 10^(-25)),digits=5)%*%Sum_Wr2
comb_rho3 = round(solve(SumW3, tol = 10^(-25)),digits=5)%*%Sum_Wr3

#Some for epsilon machine are greater than 1
# Old_comb_rho1 = comb_rho1
# Old_comb_rho2 = comb_rho2
# Old_comb_rho3 = comb_rho3
# 
# #Making values greater than 1, one
# comb_rho1[comb_rho1 > 1] = 1
# comb_rho2[comb_rho2 > 1] = 1
# comb_rho3[comb_rho3 > 1] = 1


#Making the combiened CC, MOC and MSOC matrices#
comb_CC = diag(rep(1,G))
comb_MOC = diag(rep(1,G))
comb_MSOC = diag(rep(1,G))
k = 1
for (i in 1:(G-1)){
  for (j in (i+1):G){
    comb_CC[i,j] = comb_rho1[k]
    comb_CC[j,i] = comb_rho1[k]
    comb_MOC[i,j] = comb_rho2[k]
    comb_MOC[j,i] = comb_rho2[k]
    comb_MSOC[i,j] = comb_rho3[k]
    comb_MSOC[j,i] = comb_rho3[k]
    k = k+1
  }
  print(i)
}










Results_evaluator <- function(Meta_cors,Meta_cors_Diag,cors){
  ARI = rep(NA,10+3*S)
  names(ARI) = c("Meta_CC","Meta_MOC","Meta_MSOC",paste(c("Meta_CC","Meta_MOC","Meta_MSOC"),"_Diag"),
                 paste0("CC-",1:5),paste0("MOC-",1:5),paste0("MSOC-",1:5),"WGCNA(1)-1","WGCNA(2)-1", 
                 "K-M(1)-1", "K-M(2)-1")
  Clust_results = list()
  comb_CC = Meta_cors$comb_CC
  comb_MOC = Meta_cors$comb_MOC
  comb_MSOC = Meta_cors$comb_MSOC
  omic = cors[[1]]$omic#Omic data of the 1st study
  
  ##The Meta-CC method##
  sft = pickSoftThreshold.fromSimilarity(similarity=abs(comb_CC), powerVector = powers, verbose = 0)
  
  connectivity = sft$fitIndices$mean.k.
  k = length(connectivity)
  percentChange = - (connectivity[-1]-connectivity[-k])/connectivity[-k]
  percentChange2 = - (percentChange[-1]-percentChange[-(k-1)])/percentChange[-(k-1)]
  
  softPower = which.max(percentChange2) + 2
  adjacency = abs(comb_CC)^softPower
  TOM = TOMsimilarity(as.matrix(adjacency))
  dissTOM = 1-TOM
  
  geneTree = hclust(as.dist(dissTOM), method = "average")
  minModuleSize = minClusterSize; #minimum number of clusters
  # Module identification using dynamic tree cut:
  dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                              deepSplit = 2, pamRespectsDendro = FALSE,
                              minClusterSize = minModuleSize)
  
  ###Considering each unassigned gene as a cluster####
  names(dynamicMods) = NULL
  unclass_indexes = dynamicMods == 0
  n_unclass = sum(unclass_indexes)
  if (n_unclass != 0){
    dynamicMods[unclass_indexes] = max(dynamicMods) + (1:n_unclass)
  }
  
  ARI[1] = adjustedRandIndex(TrueClass,dynamicMods)
  Clust_results[[1]] = dynamicMods
  
  
  ##The Meta-MOC method##
  rm(sft,connectivity,k,percentChange,percentChange2,softPower,adjacency,TOM,dissTOM,geneTree,dynamicMods)
  sft = pickSoftThreshold.fromSimilarity(similarity=abs(comb_MOC), powerVector = powers, verbose = 0)
  
  connectivity = sft$fitIndices$mean.k.
  k = length(connectivity)
  percentChange = - (connectivity[-1]-connectivity[-k])/connectivity[-k]
  percentChange2 = - (percentChange[-1]-percentChange[-(k-1)])/percentChange[-(k-1)]
  
  softPower = which.max(percentChange2) + 2
  adjacency = abs(comb_MOC)^softPower
  TOM = TOMsimilarity(as.matrix(adjacency))
  dissTOM = 1-TOM
  
  geneTree = hclust(as.dist(dissTOM), method = "average")
  minModuleSize = minClusterSize; #minimum number of clusters
  # Module identification using dynamic tree cut:
  dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                              deepSplit = 2, pamRespectsDendro = FALSE,
                              minClusterSize = minModuleSize)
  
  ###Considering each unassigned gene as a cluster####
  names(dynamicMods) = NULL
  unclass_indexes = dynamicMods == 0
  n_unclass = sum(unclass_indexes)
  if (n_unclass != 0){
    dynamicMods[unclass_indexes] = max(dynamicMods) + (1:n_unclass)
  }
  
  ARI[2] = adjustedRandIndex(TrueClass,dynamicMods)
  Clust_results[[2]] = dynamicMods
  
  
  ##The Meta-MSOC method##
  rm(sft,connectivity,k,percentChange,percentChange2,softPower,adjacency,TOM,dissTOM,geneTree,dynamicMods)
  sft = pickSoftThreshold.fromSimilarity(similarity=abs(comb_MSOC), powerVector = powers, verbose = 0)
  
  connectivity = sft$fitIndices$mean.k.
  k = length(connectivity)
  percentChange = - (connectivity[-1]-connectivity[-k])/connectivity[-k]
  percentChange2 = - (percentChange[-1]-percentChange[-(k-1)])/percentChange[-(k-1)]
  
  softPower = which.max(percentChange2) + 2
  adjacency = abs(comb_MSOC)^softPower
  TOM = TOMsimilarity(as.matrix(adjacency))
  dissTOM = 1-TOM
  
  geneTree = hclust(as.dist(dissTOM), method = "average")
  minModuleSize = minClusterSize; #minimum number of clusters
  # Module identification using dynamic tree cut:
  dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                              deepSplit = 2, pamRespectsDendro = FALSE,
                              minClusterSize = minModuleSize)
  
  ###Considering each unassigned gene as a cluster####
  names(dynamicMods) = NULL
  unclass_indexes = dynamicMods == 0
  n_unclass = sum(unclass_indexes)
  if (n_unclass != 0){
    dynamicMods[unclass_indexes] = max(dynamicMods) + (1:n_unclass)
  }
  
  ARI[3] = adjustedRandIndex(TrueClass,dynamicMods)
  Clust_results[[3]] = dynamicMods
  
  ###For results from the diagonal correlation covariance matrix####
  rm(comb_CC,comb_MOC,comb_MSOC)
  comb_CC = Meta_cors_Diag$comb_CC
  comb_MOC = Meta_cors_Diag$comb_MOC
  comb_MSOC = Meta_cors_Diag$comb_MSOC
  
  ##The Meta-CC method##
  rm(sft,connectivity,k,percentChange,percentChange2,softPower,adjacency,TOM,dissTOM,geneTree,dynamicMods)
  sft = pickSoftThreshold.fromSimilarity(similarity=abs(comb_CC), powerVector = powers, verbose = 0)
  
  connectivity = sft$fitIndices$mean.k.
  k = length(connectivity)
  percentChange = - (connectivity[-1]-connectivity[-k])/connectivity[-k]
  percentChange2 = - (percentChange[-1]-percentChange[-(k-1)])/percentChange[-(k-1)]
  
  softPower = which.max(percentChange2) + 2
  adjacency = abs(comb_CC)^softPower
  TOM = TOMsimilarity(as.matrix(adjacency))
  dissTOM = 1-TOM
  
  geneTree = hclust(as.dist(dissTOM), method = "average")
  minModuleSize = minClusterSize; #minimum number of clusters
  # Module identification using dynamic tree cut:
  dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                              deepSplit = 2, pamRespectsDendro = FALSE,
                              minClusterSize = minModuleSize)
  
  ###Considering each unassigned gene as a cluster####
  names(dynamicMods) = NULL
  unclass_indexes = dynamicMods == 0
  n_unclass = sum(unclass_indexes)
  if (n_unclass != 0){
    dynamicMods[unclass_indexes] = max(dynamicMods) + (1:n_unclass)
  }
  
  ARI[4] = adjustedRandIndex(TrueClass,dynamicMods)
  Clust_results[[4]] = dynamicMods
  
  
  ##The Meta-MOC method##
  rm(sft,connectivity,k,percentChange,percentChange2,softPower,adjacency,TOM,dissTOM,geneTree,dynamicMods)
  sft = pickSoftThreshold.fromSimilarity(similarity=abs(comb_MOC), powerVector = powers, verbose = 0)
  
  connectivity = sft$fitIndices$mean.k.
  k = length(connectivity)
  percentChange = - (connectivity[-1]-connectivity[-k])/connectivity[-k]
  percentChange2 = - (percentChange[-1]-percentChange[-(k-1)])/percentChange[-(k-1)]
  
  softPower = which.max(percentChange2) + 2
  adjacency = abs(comb_MOC)^softPower
  TOM = TOMsimilarity(as.matrix(adjacency))
  dissTOM = 1-TOM
  
  geneTree = hclust(as.dist(dissTOM), method = "average")
  minModuleSize = minClusterSize; #minimum number of clusters
  # Module identification using dynamic tree cut:
  dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                              deepSplit = 2, pamRespectsDendro = FALSE,
                              minClusterSize = minModuleSize)
  
  ###Considering each unassigned gene as a cluster####
  names(dynamicMods) = NULL
  unclass_indexes = dynamicMods == 0
  n_unclass = sum(unclass_indexes)
  if (n_unclass != 0){
    dynamicMods[unclass_indexes] = max(dynamicMods) + (1:n_unclass)
  }
  
  ARI[5] = adjustedRandIndex(TrueClass,dynamicMods)
  Clust_results[[5]] = dynamicMods
  
  
  ##The Meta-MSOC method##
  rm(sft,connectivity,k,percentChange,percentChange2,softPower,adjacency,TOM,dissTOM,geneTree,dynamicMods)
  sft = pickSoftThreshold.fromSimilarity(similarity=abs(comb_MSOC), powerVector = powers, verbose = 0)
  
  connectivity = sft$fitIndices$mean.k.
  k = length(connectivity)
  percentChange = - (connectivity[-1]-connectivity[-k])/connectivity[-k]
  percentChange2 = - (percentChange[-1]-percentChange[-(k-1)])/percentChange[-(k-1)]
  
  softPower = which.max(percentChange2) + 2
  adjacency = abs(comb_MSOC)^softPower
  TOM = TOMsimilarity(as.matrix(adjacency))
  dissTOM = 1-TOM
  
  geneTree = hclust(as.dist(dissTOM), method = "average")
  minModuleSize = minClusterSize; #minimum number of clusters
  # Module identification using dynamic tree cut:
  dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                              deepSplit = 2, pamRespectsDendro = FALSE,
                              minClusterSize = minModuleSize)
  
  ###Considering each unassigned gene as a cluster####
  names(dynamicMods) = NULL
  unclass_indexes = dynamicMods == 0
  n_unclass = sum(unclass_indexes)
  if (n_unclass != 0){
    dynamicMods[unclass_indexes] = max(dynamicMods) + (1:n_unclass)
  }
  
  ARI[6] = adjustedRandIndex(TrueClass,dynamicMods)
  Clust_results[[6]] = dynamicMods
  
  ##CC for the 1st study##
  rm(sft,connectivity,k,percentChange,percentChange2,softPower,adjacency,TOM,dissTOM,geneTree,dynamicMods)
  for (s in 1:S){
    corM = cors[[s]]$cancorM #Considering the sth study
    sft = pickSoftThreshold.fromSimilarity(similarity=as.matrix(corM), powerVector = powers, verbose = 0)
    
    connectivity = sft$fitIndices$mean.k.
    k = length(connectivity)
    percentChange = - (connectivity[-1]-connectivity[-k])/connectivity[-k]
    percentChange2 = - (percentChange[-1]-percentChange[-(k-1)])/percentChange[-(k-1)]
    
    softPower = which.max(percentChange2) + 2
    adjacency = corM^softPower
    TOM = TOMsimilarity(as.matrix(adjacency))
    dissTOM = 1-TOM
    
    geneTree = hclust(as.dist(dissTOM), method = "average")
    minModuleSize = minClusterSize; #minimum number of clusters
    # Module identification using dynamic tree cut:
    dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                                deepSplit = 2, pamRespectsDendro = FALSE,
                                minClusterSize = minModuleSize)
    
    ###Considering each unassigned gene as a cluster####
    names(dynamicMods) = NULL
    unclass_indexes = dynamicMods == 0
    n_unclass = sum(unclass_indexes)
    if (n_unclass != 0){
      dynamicMods[unclass_indexes] = max(dynamicMods) + (1:n_unclass)
    }
    
    ARI[6+s] = adjustedRandIndex(TrueClass,dynamicMods)
    Clust_results[[6+s]] = dynamicMods
    print(paste0(s,"th study CC"))
    rm(corM,sft,connectivity,k,percentChange,percentChange2,softPower,adjacency,TOM,dissTOM,geneTree,dynamicMods)
  }
  
  
  ##Maximum omic correlation for the 1st study##
  
  for (s in 1:S){
    corM2 = cors[[s]]$Max_Omic_Cor
    sft = pickSoftThreshold.fromSimilarity(similarity=as.matrix(abs(corM2)), powerVector = powers, verbose = 0)
    
    connectivity = sft$fitIndices$mean.k.
    k = length(connectivity)
    percentChange = - (connectivity[-1]-connectivity[-k])/connectivity[-k]
    percentChange2 = - (percentChange[-1]-percentChange[-(k-1)])/percentChange[-(k-1)]
    
    softPower = which.max(percentChange2) + 2
    adjacency = abs(corM2)^softPower
    TOM = TOMsimilarity(as.matrix(adjacency))
    dissTOM = 1-TOM
    
    geneTree = hclust(as.dist(dissTOM), method = "average")
    minModuleSize = minClusterSize; #minimum number of clusters
    # Module identification using dynamic tree cut:
    dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                                deepSplit = 2, pamRespectsDendro = FALSE,
                                minClusterSize = minModuleSize)
    
    ###Considering each unassigned gene as a cluster####
    names(dynamicMods) = NULL
    unclass_indexes = dynamicMods == 0
    n_unclass = sum(unclass_indexes)
    if (n_unclass != 0){
      dynamicMods[unclass_indexes] = max(dynamicMods) + (1:n_unclass)
    }
    
    ARI[6+S+s] = adjustedRandIndex(TrueClass,dynamicMods)
    Clust_results[[6+S+s]] = dynamicMods
    print(paste0(s,"th study MOC"))
    rm(sft,connectivity,k,percentChange,percentChange2,softPower,adjacency,TOM,dissTOM,geneTree,dynamicMods)
  }
  
  
  ##Maximum Same Omic correlation for the 1st study##
  for (s in 1:S){
    corM3 = cors[[s]]$Max_Same_Omic_Cor
    sft = pickSoftThreshold.fromSimilarity(similarity=as.matrix(abs(corM3)), powerVector = powers, verbose = 0)
    
    connectivity = sft$fitIndices$mean.k.
    k = length(connectivity)
    percentChange = - (connectivity[-1]-connectivity[-k])/connectivity[-k]
    percentChange2 = - (percentChange[-1]-percentChange[-(k-1)])/percentChange[-(k-1)]
    
    softPower = which.max(percentChange2) + 2
    adjacency = abs(corM3)^softPower
    TOM = TOMsimilarity(as.matrix(adjacency))
    dissTOM = 1-TOM
    
    geneTree = hclust(as.dist(dissTOM), method = "average")
    minModuleSize = minClusterSize; #minimum number of clusters
    # Module identification using dynamic tree cut:
    dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                                deepSplit = 2, pamRespectsDendro = FALSE,
                                minClusterSize = minModuleSize)
    
    ###Considering each unassigned gene as a cluster####
    names(dynamicMods) = NULL
    unclass_indexes = dynamicMods == 0
    n_unclass = sum(unclass_indexes)
    if (n_unclass != 0){
      dynamicMods[unclass_indexes] = max(dynamicMods) + (1:n_unclass)
    }
    
    ARI[6+2*S+s] = adjustedRandIndex(TrueClass,dynamicMods)
    Clust_results[[6+2*S+s]] = dynamicMods
    print(paste0(s,"th study MSOC"))
    rm(sft,connectivity,k,percentChange,percentChange2,softPower,adjacency,TOM,dissTOM,geneTree,dynamicMods)
  }
  
  ###One-omic methods### with omic 1 only of the 1st study
  single_omic1 = omic[,seq(1,total_individuals,3)]
  M1 = abs(cor(single_omic1)) #Correlation matrix from omic 1
  # single_omic2 = omic[,seq(2,3000,3)]
  # single_omic3 = omic[,seq(3,3000,3)]
  
  
  #The first omic#
  sft = pickSoftThreshold.fromSimilarity(similarity=M1, powerVector = powers, verbose = 0)
  
  #One-omic with the new cut function#
  connectivity = sft$fitIndices$mean.k.
  k = length(connectivity)
  percentChange = - (connectivity[-1]-connectivity[-k])/connectivity[-k]
  percentChange2 = - (percentChange[-1]-percentChange[-(k-1)])/percentChange[-(k-1)]
  softPower = which.max(percentChange2) + 2
  adjacency = M1^softPower
  TOM = TOMsimilarity(as.matrix(adjacency))
  dissTOM = 1-TOM
  geneTree = hclust(as.dist(dissTOM), method = "average")
  minModuleSize = minClusterSize; #minimum number of clusters
  # Module identification using dynamic tree cut:
  dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                              deepSplit = 2, pamRespectsDendro = FALSE,
                              minClusterSize = minModuleSize)
  
  ###Considering each unassigned gene as a cluster####
  names(dynamicMods) = NULL
  unclass_indexes = dynamicMods == 0
  n_unclass = sum(unclass_indexes)
  if (n_unclass != 0){
    dynamicMods[unclass_indexes] = max(dynamicMods) + (1:n_unclass)
  }
  
  ARI[6+3*S+1] = adjustedRandIndex(TrueClass,dynamicMods)
  Clust_results[[6+3*S+1]] = dynamicMods
  
  
  ##One-omic traditional WGCNA##
  rm(softPower,adjacency,TOM,dissTOM,geneTree,dynamicMods)
  softPower = PowerCut3(sft$fitIndices$Power,sft$fitIndices$truncated.R.sq,sft$fitIndices$mean.k.,0.6,1/3)
  adjacency = M1^softPower
  TOM = TOMsimilarity(as.matrix(adjacency))
  dissTOM = 1-TOM
  
  geneTree = hclust(as.dist(dissTOM), method = "average")
  minModuleSize = minClusterSize; #minimum number of clusters
  # Module identification using dynamic tree cut:
  dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                              deepSplit = 2, pamRespectsDendro = FALSE,
                              minClusterSize = minModuleSize)
  
  ###Considering each unassigned gene as a cluster####
  names(dynamicMods) = NULL
  unclass_indexes = dynamicMods == 0
  n_unclass = sum(unclass_indexes)
  if (n_unclass != 0){
    dynamicMods[unclass_indexes] = max(dynamicMods) + (1:n_unclass)
  }
  
  ARI[6+3*S+2] = adjustedRandIndex(TrueClass,dynamicMods)
  Clust_results[[6+3*S+2]] = dynamicMods
  
  ##K-means with the true k##
  set.seed(12345)
  Kmean_result = kmeans(t(scale(single_omic1)),centers = length(table(TrueClass)),nstart = 25)
  ARI[6+3*S+3] = adjustedRandIndex(TrueClass,Kmean_result$cluster)
  Clust_results[[6+3*S+3]] = Kmean_result$cluster
  
  
  ##K-means with the gap statistics for choosing the optimal k between 1:130##
  set.seed(12345)
  gaps = clusGap(t(scale(single_omic1)), kmeans, K.max=length(table(TrueClass))+10, B = 100, d.power = 1)
  rm(Kmean_result)
  opt_k = maxSE(gaps$Tab[,"gap"], gaps$Tab[,"SE.sim"],method = "globalSEmax")
  set.seed(12345)
  Kmean_result = kmeans(t(scale(single_omic1)),centers = opt_k,nstart = 25)
  ARI[6+3*S+4] = adjustedRandIndex(TrueClass,Kmean_result$cluster)
  Clust_results[[6+3*S+4]] = Kmean_result$cluster
  
  names(Clust_results) = c("Meta_CC","Meta_MOC","Meta_MSOC",paste(c("Meta_CC","Meta_MOC","Meta_MSOC"),"_Diag"),
                           paste0("CC-",1:5),paste0("MOC-",1:5),paste0("MSOC-",1:5),"WGCNA(1)-1","WGCNA(2)-1", 
                           "K-M(1)-1", "K-M(2)-1")
  
  list(ARI=ARI,Clusters=Clust_results)
}

