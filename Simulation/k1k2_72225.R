# Running third: 7/23/25 8:15 PM
# corrected dimension selection
##k1=k2=3,...,10, 
##-->consider k1,k2 known since this always has to be estimated,
##   we want to isolate effect of KR vs MASE.
##m=101, n=1000
##trueK = max(k1,k2),k1+k2,floor(k1*k2/2),k1*k2 #note that early versions have floor((k1+k2)/2), which may be <max(k1,k2).
##sigma = 0.01,0.05,0.1,0.25,0.5,1
##Z/U, KR/MASE
##known k1,k2, un/known K


library(irlba)
library(mclust)
library(MGLM)
library(ggplot2)
library(GGally)
library(pracma)
library(MASS)
library(ggpubr)
library(RColorBrewer)
library(reshape2)

pacman::p_load(irlba,mclust,pracma,MASS,reshape2)
require(doParallel)
cluster <- makeCluster(10)
registerDoParallel(cluster)

class2mat <- function(tauvec){
  tauvec<- as.numeric(tauvec)
  k <- length(unique(tauvec))
  taumat <- matrix(rep(0,length(tauvec)*k),nrow=length(tauvec))
  for(i in 1:length(tauvec)){
    taumat[i,tauvec[i]]<- 1
  }
  return(taumat)
}

procrustes2 <- function(X, Y, type = "I") {
  if(type == "C"){
    X <- X/norm(X, type = "F")*sqrt(nrow(X))
    Y <- Y/norm(Y, type = "F")*sqrt(nrow(Y))
  }
  if(type == "D"){
    tX <- rowSums(X^2)
    tX[tX <= 1e-15] <- 1
    tY <- rowSums(Y^2)
    tY[tY <= 1e-15] <- 1
    X <- X/sqrt(tX)
    Y <- Y/sqrt(tY)
  }
  
  tmp <- t(X) %*% Y
  tmp.svd <- svd(tmp)
  W <- tmp.svd$u %*% t(tmp.svd$v)
  newX <- X %*% W
  return(list(newX = newX, error = norm(newX-Y, type="F"), W = W))
}

drawJointClustering <- function(n,Propmat){
  k <- dim(Propmat)
  Z <- t(rmultinom(n,1,c(Propmat)))
  Zarray <- array(c(t(Z)),dim=c(k[1],k[2],n))
  Z1 <- apply(Zarray,c(3,1),sum)
  Z2 <- apply(Zarray,c(3,2),sum)
  return(list(Z=Z,Z1=Z1,Z2=Z2))
}

getElbows <- function(dat, n = 3, threshold = FALSE, plot = FALSE, main="") {
  ## Given a decreasingly sorted vector, return the given number of elbows
  ##
  ## Args:
  ##   dat: a input vector (e.g. a vector of standard deviations), or a input feature matrix.
  ##   n: the number of returned elbows.
  ##   threshold: either FALSE or a number. If threshold is a number, then all
  ##   the elements in d that are not larger than the threshold will be ignored.
  ##   plot: logical. When T, it depicts a scree plot with highlighted elbows.
  ##
  ## Return:
  ##   q: a vector of length n.
  ##
  ## Reference:
  ##   Zhu, Mu and Ghodsi, Ali (2006), "Automatic dimensionality selection from
  ##   the scree plot via the use of profile likelihood", Computational
  ##   Statistics & Data Analysis, Volume 51 Issue 2, pp 918-930, November, 2006. 
  
  #  if (is.unsorted(-d))
  
  
  if (is.matrix(dat)) {
    d <- sort(apply(dat,2,sd), decreasing=TRUE)
  } else {
    d <- sort(dat,decreasing=TRUE)
  }
  
  if (!is.logical(threshold))
    d <- d[d > threshold]
  
  p <- length(d)
  if (p == 0)
    stop(paste("d must have elements that are larger than the threshold ",
               threshold), "!", sep="")
  
  lq <- rep(0.0, p)                     # log likelihood, function of q
  for (q in 1:p) {
    mu1 <- mean(d[1:q])
    mu2 <- mean(d[-(1:q)])              # = NaN when q = p
    sigma2 <- (sum((d[1:q] - mu1)^2) + sum((d[-(1:q)] - mu2)^2)) /
      (p - 1 - (q < p))
    lq[q] <- sum( dnorm(  d[1:q ], mu1, sqrt(sigma2), log=TRUE) ) +
      sum( dnorm(d[-(1:q)], mu2, sqrt(sigma2), log=TRUE) )
  }
  
  q <- which.max(lq)
  if(is.na(q)){
    q <- p
  }
  if (n > 1 && q < (p-1)) {
    q <- c(q, q + getElbows(d[(q+1):p], n-1, plot=FALSE))
  }
  
  if (plot==TRUE) {
    if (is.matrix(dat)) {
      sdv <- d # apply(dat,2,sd)
      plot(sdv,type="b",xlab="dim",ylab="stdev",main=main)
      points(q,sdv[q],col=2,pch=19)
    } else {
      plot(dat, type="b",main=main)
      points(q,dat[q],col=2,pch=19)
    }
  }
  
  return(q)
}

compareWithTruth <- function(embeddings,num_clusts,truth){
  eH <- hc(embeddings)
  CR <- c(hclass(eH,num_clusts))
  return(adjustedRandIndex(truth,CR))
}

getProbMatrix <- function(k1,k2,trueK){
  #trueK must be >= max(k1,k2)
  pvec <- sqrt(runif(k1*k2))
  pvec <- sort(pvec,decreasing=TRUE)[sample(1:trueK)]
  pvec <- pvec/sum(pvec)
  if(k1<=k2){
    rowk <- k1
    colk <- k2
  }else{
    rowk <- k2
    colk <- k1
  }
  P <- matrix(rep(0,k1*k2),nrow=rowk)
  diag(P[,1:rowk]) <- pvec[1:rowk]
  if(colk>rowk){
    for(i in (rowk+1):colk){
      P[sample(1:rowk,1),i] <- pvec[i]
    }
  }
  if(trueK>colk){
    P[sample(which(P==0),trueK-colk)] <- pvec[(colk+1):trueK]
  }
  if(k1>k2){
    P <- t(P)
  }
  return(P)
}

kvec <- 3:10
m <- 101
sigmavec <- c(0.01,0.05,0.1,0.25,0.5,1)

num_k1 <- length(kvec)
num_K <- 4
num_sigma <- length(sigmavec)
dimvec <- c(num_k1,num_k1,4,num_sigma)
nsetups <- num_k1/2*(num_k1+1)*4*num_sigma

ntrials <- 100
mult <- 1.96/sqrt(ntrials)
#order is mean(KR)-mult*sd, mean(KR), mean(KR)+mult*sd, 
#         mean(MASE)-mult*sd, mean(MASE), mean(MASE)+mult*sd
Zresults.known <- array(0,dim=c(dimvec,6))
Uresults.known <- array(0,dim=c(dimvec,6))
Zresults.unknown <- array(0,dim=c(dimvec,6))
Uresults.unknown <- array(0,dim=c(dimvec,6))
Zkestresults <- array(0,dim=c(dimvec,6))
Ukestresults <- array(0,dim=c(dimvec,6))

n <- 1000
s<- 1

verybeginning <- Sys.time()

for(k1idx in 1:num_k1){
  k1 <- kvec[k1idx]
  for(k2idx in k1idx:num_k1){
    k2 <- kvec[k2idx]
    trueKvec <- c(max(k1,k2),k1+k2,floor(k1*k2/2),k1*k2)
    for(Kidx in 1:4){
      trueK <- trueKvec[Kidx]
      for(sigmaidx in 1:num_sigma){
        sigma <- sigmavec[sigmaidx]
        #Order is Zkr, Zmase, Ukr, Umase
        trialresultsall <- matrix(0,nrow=ntrials,ncol=12)
        stime <- Sys.time()
        trialresultsall <- foreach (i = 1:ntrials,.combine=rbind,.packages=c("MASS","MGLM","mclust")) %dopar% {
          P <- getProbMatrix(k1,k2,trueK)
          jC <- drawJointClustering(n,P)
          truelabels <- c(jC$Z%*%(1:dim(jC$Z)[2]))
          Y1 <- jC$Z1%*%mvrnorm(k1,mu=rep(0,m),Sigma=diag(m))+
            mvrnorm(n,mu=rep(0,m),Sigma=sigma*diag(m))
          Y2 <- jC$Z2%*%mvrnorm(k2,mu=rep(0,m),Sigma=diag(m))+
            mvrnorm(n,mu=rep(0,m),Sigma=sigma*diag(m))
          svd1 <- svd(Y1)
          svd2 <- svd(Y2)
          U1 <- svd1$u[,1:k1]
          U2 <- svd2$u[,1:k2]
          Z1 <- class2mat(kmeans(U1,centers=k1)$cluster)
          Z2 <- class2mat(kmeans(U2,centers=k2)$cluster)
          Zkr <- kr(Z1,Z2)
          Zmase <- cbind(Z1,Z2)
          Ukr <- kr(as.matrix(U1),as.matrix(U2))
          Umase <- cbind(U1,U2)
          svd.Z.kr <- svd(Zkr)
          svd.Z.mase <- svd(Zmase)
          svd.U.kr <- svd(Ukr)
          svd.U.mase <- svd(Umase)
          #known cluster numbers
          known_results <- c(compareWithTruth(svd.Z.kr$u[,1:trueK],trueK,truelabels),
                             compareWithTruth(svd.Z.mase$u[,1:min(dim(Zmase)[2],trueK)],trueK,truelabels),
                             compareWithTruth(svd.U.kr$u[,1:trueK],trueK,truelabels),
                             compareWithTruth(svd.U.mase$u[,1:min(dim(Umase)[2],trueK)],trueK,truelabels))
          #unknown cluster numbers
          d1 <- getElbows(svd.Z.kr$d,n=2)
          d1 <- d1[length(d1)]
          d2 <- getElbows(svd.Z.mase$d,n=2)
          d2 <- d2[length(d2)]
          d3 <- getElbows(svd.U.kr$d,n=2)
          d3 <- d3[length(d3)]
          d4 <- getElbows(svd.U.mase$d,n=2)
          d4 <- d4[length(d4)]
          estK.Z.kr <- min(d1,sum(svd.Z.kr$d>0))
          estK.Z.mase <- min(d2,sum(svd.Z.mase$d>0))
          estK.U.kr <- min(d3,sum(svd.U.kr$d>0))
          estK.U.mase <- min(d4,sum(svd.U.mase$d>0))
          unknown_results <- c(compareWithTruth(svd.Z.kr$u[,1:estK.Z.kr],estK.Z.kr,truelabels),
                               compareWithTruth(svd.Z.mase$u[,1:estK.Z.mase],estK.Z.mase,truelabels),
                               compareWithTruth(svd.U.kr$u[,1:estK.U.kr],estK.U.kr,truelabels),
                               compareWithTruth(svd.U.mase$u[,1:estK.U.mase],estK.U.mase,truelabels))
          estk_results <- c(abs(estK.Z.kr-trueK),abs(estK.Z.mase-trueK),
                            abs(estK.U.kr-trueK),abs(estK.U.mase-trueK))
          c(known_results,unknown_results,estk_results)
        }
        means <- apply(trialresultsall,2,mean)
        sds <- apply(trialresultsall,2,sd)
        Zresults.known[k1idx,k2idx,Kidx,sigmaidx,] <- 
          c(max(c(means[1]-mult*sds[1],0)),means[1],min(c(1,means[1]+mult*sds[1])),
            max(c(means[2]-mult*sds[2],0)),means[2],min(c(1,means[2]+mult*sds[2])))
        Zresults.known[k2idx,k1idx,Kidx,sigmaidx,] <- Zresults.known[k1idx,k2idx,Kidx,sigmaidx,]
        Uresults.known[k1idx,k2idx,Kidx,sigmaidx,] <- 
          c(max(c(means[3]-mult*sds[3],0)),means[3],min(c(1,means[3]+mult*sds[3])),
            max(c(means[4]-mult*sds[4],0)),means[4],min(c(1,means[4]+mult*sds[4])))
        Uresults.known[k2idx,k1idx,Kidx,sigmaidx,] <- Uresults.known[k1idx,k2idx,Kidx,sigmaidx,]
        Zresults.unknown[k1idx,k2idx,Kidx,sigmaidx,] <- 
          c(max(c(means[5]-mult*sds[5],0)),means[5],min(c(1,means[5]+mult*sds[5])),
            max(c(means[6]-mult*sds[6],0)),means[6],min(c(1,means[6]+mult*sds[6])))
        Zresults.unknown[k2idx,k1idx,Kidx,sigmaidx,] <- Zresults.unknown[k1idx,k2idx,Kidx,sigmaidx,]
        Uresults.unknown[k1idx,k2idx,Kidx,sigmaidx,] <- 
          c(max(c(means[7]-mult*sds[7],0)),means[7],min(c(1,means[7]+mult*sds[7])),
            max(c(means[8]-mult*sds[8],0)),means[8],min(c(1,means[8]+mult*sds[8])))
        Uresults.unknown[k2idx,k1idx,Kidx,sigmaidx,] <- Uresults.unknown[k1idx,k2idx,Kidx,sigmaidx,]
        
        Zkestresults[k1idx,k2idx,Kidx,sigmaidx,] <-
          c(max(c(means[9]-mult*sds[9],0)),means[9],means[9]+mult*sds[9],
            max(c(means[10]-mult*sds[10],0)),means[10],means[10]+mult*sds[10])
        Zkestresults[k2idx,k1idx,Kidx,sigmaidx,] <- Zkestresults[k1idx,k2idx,Kidx,sigmaidx,]
        Ukestresults[k1idx,k2idx,Kidx,sigmaidx,] <-
          c(max(c(means[11]-mult*sds[11],0)),means[11],means[11]+mult*sds[11],
            max(c(means[12]-mult*sds[12],0)),means[12],means[12]+mult*sds[12])
        Ukestresults[k2idx,k1idx,Kidx,sigmaidx,] <- Ukestresults[k1idx,k2idx,Kidx,sigmaidx,]
        cat('Setup: ', s, ' of ', nsetups,'\n')
        cat('k1: ', k1idx, ' of ', num_k1,'\n')
        cat('k2: ', k2idx+1-k1idx, ' of ', num_k1+1-k1idx,'\n')
        cat('trueK: ', Kidx, ' of ', 4,'\n')
        cat('sigma: ', sigmaidx, ' of ', num_sigma,'\n')
        etime <- Sys.time()
        cat('Elapsed time: ',difftime(etime,stime,units="secs")[[1]],' seconds.\n')
        cat('Estimated time remaining: ',difftime(etime,verybeginning,units="secs")[[1]]*(nsetups-s)/s,' seconds.\n')
        s<- s+1
      }
    }
  }
}
stopCluster(cl=cluster)
save(kvec,m,sigmavec,ntrials,mult,Zresults.known,Uresults.known,
     Zresults.unknown,Uresults.unknown,Zkestresults,Ukestresults,
     file="ExperimentResults7-22-25_varyingk1k2.RData")







