library('Matrix')
library('Rcpp')
library('BayesLogit')
sourceCpp('Gibbs.cpp')
############### EM inference function
EMinfer <- function(pvals){
  gene.num <- length(pvals)
  
  alpha0 <- 1
  alpha1 <- 0.2
  pi0 <- 0.6
  pi1 <- 0.4
  # initial EM for estimating parameters of p-value distribution
  for(i in 1:2000){
    tmp1 <- pi1*alpha1*pvals^(alpha1-1)
    tmp0 <- pi0*alpha0*(pvals)^(alpha0-1)
    proZ <- tmp1/(tmp1+tmp0)
    
    pi0 <- sum(1-proZ)/length(proZ)
    pi1 <- sum(proZ) / length(proZ)
    
    alpha0 <- - sum(1 - proZ)/ sum( (1-proZ)* log(pvals) )
    alpha1 <- - sum(proZ)/ sum( proZ* log(pvals) )
    
  }
  print(sprintf("inital estimation via EM : alpha0 %f alpha1 %f pi0 %f pi1 %f",alpha0,alpha1,pi0,pi1))
  return(proZ)
}

############################################## simulation function
simulation <- function(simulated.alpha0=1,simulated.alpha1=0.2,simulated.gamma=-2,simulated.beta,pvals,network){
  Neighbor.index  <- network$index
  Neighbor.weight <- network$weight
  Neighbor.count  <- network$count
  Neighbor.sumW   <- network$sumW
  
  simulated.Z <- rbinom(gene.num,1,prob = 1/(1+exp(-simulated.gamma)))
  simulated.Neighbor.statistics <- matrix(0.0, nrow=gene.num, ncol = K)
  K <- length(Neighbor.weight)
  z.num <- sum(simulated.Z)
  
  for(k in 1:K){
    for(i in 1:length(pvals)){
      if(Neighbor.count[k,i] > 0){
        tmpZ  <- simulated.Z[Neighbor.index[[k]][[i]]]
        tmpW  <- Neighbor.weight[[k]][[i]]
        simulated.Neighbor.statistics[i,k] <- ( sum(tmpW[tmpZ==1])-sum(tmpW[tmpZ==0] ))
      }else{
        simulated.Neighbor.statistics[i,k] <- 0
      }
    }
  }
  

  
  for(l in 1:20){
    for(g in 1:length(pvals)){
      
      tmp <- simulated.Neighbor.statistics[g,] %*% simulated.beta + simulated.gamma
      proba <- 1/(1+exp(-tmp))
      before <- simulated.Z[g]
      simulated.Z[g] <- rbinom(1,1,prob = proba)
      after <- simulated.Z[g]

      if(before != after){
      
      for(k in 1:K){
        tmpW <- Neighbor.weight[[k]][[g]]
        idx <- 1
        for(i in Neighbor.index[[k]][[g]]){
          simulated.Neighbor.statistics[i,k] <- simulated.Neighbor.statistics[i,k] + (after-before)*2*tmpW[idx]
          idx <- idx + 1
          }
      }
      }
      
    }
  }
  
  for(i in 1:length(pvals)){
    if(simulated.Z[i] == 1){
      pvals[i] <- rbeta(1,simulated.alpha1,1)
    }else{
      pvals[i] <- rbeta(1,simulated.alpha0,1)
    }
  }
  
  data <- list(pvals = pvals, true.Z = simulated.Z)

  return(data)
}
