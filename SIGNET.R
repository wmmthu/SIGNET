library('Matrix')
library('Rcpp')
library('BayesLogit')
sourceCpp('Gibbs.cpp')

# SIGNET
SIGNET_infer <- function(iter=50000,verbose=FALSE,pvals,network){
  
  Neighbor.index  <- network$index
  Neighbor.weight <- network$weight
  Neighbor.count  <- network$count
  Neighbor.sumW   <- network$sumW
  K <- length(Neighbor.weight)
  
  gene.num <- length(pvals)
  Z <- matrix(0,nrow = gene.num,ncol = (iter+1))
  
  alpha0 <- 1
  alpha1 <- 0.2
  pi0 <- 0.6
  pi1 <- 0.4
  
  gamma  <- 0
  beta  <- seq(0,0,length.out = K)
  
  probI <- 1/K
  g <- 100
  
  para.list <- matrix(0,ncol=(2*K+4),nrow = (iter+1))
  para.list[1,] <- c(alpha0,alpha1,gamma,beta,beta,1)
  
  n0 <- 0
  n1 <- 0
  Neighbor.statistics <- matrix(0.0, nrow=gene.num, ncol = K)
  n.monitor <- 100
  
  # initial EM for estimating parameters of p-value distribution
  for(i in 1:200){
    tmp1 <- pi1*alpha1*pvals^(alpha1-1)
    tmp0 <- pi0*alpha0*(pvals)^(alpha0-1)
    proZ <- tmp1/(tmp1+tmp0)
    
    pi0 <- sum(1-proZ)/length(proZ)
    pi1 <- sum(proZ) / length(proZ)
    
    alpha0 <- - sum(1 - proZ)/ sum( (1-proZ)* log(pvals) )
    alpha1 <- - sum(proZ)/ sum( proZ* log(pvals) )
    
  }
  gamma <- -log(1/pi1 - 1)
  print(sprintf("inital estimation via EM : alpha0 %f alpha1 %f pi0 %f pi1 %f",alpha0,alpha1,pi0,pi1))
  para.list[1,] <- c(alpha0,alpha1,gamma,beta,beta,1)
  init.Z <- rbinom(gene.num,1,proZ)
  Z[,1] <- init.Z
  
  
  for(k in 1:K){
    for(i in 1:length(pvals)){
      if(Neighbor.count[k,i] > 0){
        tmpZ  <- Z[Neighbor.index[[k]][[i]],1]
        tmpW  <- Neighbor.weight[[k]][[i]]
        Neighbor.statistics[i,k] <- ( sum(tmpW[tmpZ==1])-sum(tmpW[tmpZ==0]))
      }else{
        Neighbor.statistics[i,k] <- 0
      }
    }
  }
  
  w <- numeric(gene.num)
  for(t in 2:(iter+1)){
    n0 <- sum(Z[,t-1]==0)
    n1 <- sum(Z[,t-1]==1)
    
    # update parameter
    # alpha0
    para.list[t,1] <- rgamma(1,shape=n0+1, rate = 1 - sum(log(pvals[Z[,t-1]==0])) )
    # alpha1
    para.list[t,2] <- rgamma(1,shape=n1+1, rate = 1 - sum(log(pvals[Z[,t-1]==1])))

    if(para.list[t,1] < para.list[t,2]){
	    para.list[t,c(1,2)] <- para.list[t,c(2,1)]
    }
    
    # beta and gamma 
    beta <- para.list[t-1,3:(K+3)]
    betaI <- c(1,para.list[t-1,(K+4):(2*K+3)])
    
    X <- cbind(1,Neighbor.statistics)
    Xr <- as.matrix(X[,which(betaI==1)])
    br <- beta[which(betaI==1)]
    
    w <- rpg(num=gene.num,h=1,z=Xr%*%br)
    vv <- diag(sum(betaI) ) * para.list[t-1,2*K+4]
    vv[1,1] <- 5
    
    V <- solve(solve(vv) + t(Xr)%*%(Xr*w))
    
    m <- V%*%(t(Xr)%*%(Z[,t-1]-1/2))
    newBeta <- m + t(chol(V))%*%rnorm(sum(betaI))
    para.list[t,which(betaI==1)+2] <- newBeta 

    para.list[t,3] <- min(para.list[t,3],0)

    para.list[t,which(betaI==0)+2] <- rnorm(sum(betaI==0),0,para.list[t-1,2*K+4]/g)
    # sample indicator
    beta <- para.list[t,3:(K+3)]
    Xrbr <- Xr%*%br
    for(i in 1:K){
      zz <- (Z[,t-1] - 1/2) / w - Xrbr
      if(betaI[i+1]==0){
        zz <- zz - beta[i+1]*X[,i+1]
      }
      p1 <- log(probI) - 1/2 * sum(zz^2 * w) + log(dnorm(beta[i+1],0,para.list[t-1,2*K+4]))
      zz <- zz + beta[i+1]*X[,i+1]
      p0 <- log(1 - probI) - 1/2 * sum(zz^2 * w) + log(dnorm(beta[i+1],0,para.list[t-1,2*K+4]/g))
      para.list[t,3+K+i] <- rbinom(1,1,prob=1/(1+exp(p0-p1)) )
    }

    # update tau
    
    curBeta <- para.list[t,4:(K+3)]
    curI    <- para.list[t,(K+4):(2*K+3)]
    para.list[t,2*K+4] <- 1/rgamma(1,shape=K/2 + 2, rate = 1/2*(sum(curBeta[curI==1]^2) + g*sum(curBeta[curI==0]^2))+1)
    
    # udpate Z sequtially
    Z[,t] <- Gibbs(Neighbor.statistics,Neighbor.index,Neighbor.weight,Neighbor.sumW,para.list[t,],pvals,Z[,t-1])
    
    if( t > 10 & (t/n.monitor) == round(t/n.monitor) & verbose){
      cat("run : ",t,"\n")
      cat("para:",para.list[t,],"\n")
      cat("probI:", probI, '\n')
    }
    
    
  }
  
  data <- data.frame(pval=pvals)
  data['localFDR'] <- 1 - apply(Z[,floor(iter/2):iter],1,mean)
  
  
  if(K>1){
    network.included.prob <- apply(para.list[floor(iter/2):iter,(K+4):(2*K+3)],2,mean)
  }else{
    network.included.prob <- mean(para.list[floor(iter/2):iter,(K+4):(2*K+3)])
  }
  result <- list(data=data,network.included.prob=network.included.prob,para.list=para.list)
  
  return(result)
}


SIGNET <- function(p_value_file, network_dir, iters = 20000, remove_HLA=TRUE, edge_threshold=0.01){
  
  # load GWAS p value file
  GWAS.p <- read.table(p_value_file,header=T,stringsAsFactors = F)
  
  # remove genes fall within HLA region
  if(remove_HLA == TRUE){
    HLA.gene.index <- which(GWAS.p[,1]=='chr6' & GWAS.p[,3] > 25000000 & GWAS.p[,2] < 35000000)
    if(length(HLA.gene.index) > 0){
      GWAS.p <- GWAS.p[-HLA.gene.index,]
    }
    Genes  <- GWAS.p[,'gene_symbol']
  }
  
  # record near gene pairs
  near.pair <- matrix(NA,ncol=2,nrow=0)
  near.pair.str <- c()
  for(i in 1:dim(GWAS.p)[1]){
    overlap.start <- GWAS.p[i,'start'] - 1000000
    overlap.end   <- GWAS.p[i,'end'] + 1000000
    
    overlap.index <- which((GWAS.p[,1]==GWAS.p[i,1]) & ( (GWAS.p[,'start']<overlap.end & GWAS.p[,'start']> overlap.start)|(GWAS.p[,'end']<overlap.end & GWAS.p[,'end']> overlap.start) )  )
    overlap.index <- overlap.index[which(overlap.index!=i)]
    near.pair <- rbind(near.pair, matrix(c(seq(i,i,length.out = length(overlap.index)), overlap.index), ncol=2,byrow = F))
    
    near.pair.str <- c(near.pair.str,paste(GWAS.p[i,'gene_symbol'], GWAS.p[overlap.index,'gene_symbol']))
    near.pair.str <- c(near.pair.str,paste(GWAS.p[overlap.index,'gene_symbol'], GWAS.p[i,'gene_symbol']))
  }
  
  # load tissue-specific gene networks
  Network <- list()
  W <- list()
  net.filedir  <- network_dir
  net.filelist <- list.files(net.filedir)
  
  network.names <- c()
  edge.threshold <- edge_threshold
  for(i in 1:length(net.filelist)){
    print(sprintf('loading network %d',i))
    tmp <- strsplit(net.filelist[i],'-')
    network.names <- c(network.names , net.filelist[i])
    Network[[i]] <- read.table(sprintf('%s/%s',net.filedir, net.filelist[i]), stringsAsFactors = F)
    Network[[i]] <- Network[[i]][which(Network[[i]][,3]>edge.threshold),]
    
    network.pair <- paste(Network[[i]][,1],Network[[i]][,2])
    Network[[i]] <- Network[[i]][which(is.na(match(network.pair,near.pair.str))),]
    
    W[[i]] <- Matrix(0,nrow=length(Genes),ncol=length(Genes),sparse=TRUE)
    
    index.1 <- match(Network[[i]][,1],Genes)
    index.2 <- match(Network[[i]][,2],Genes)
    
    both.index <- which(!is.na(index.1) & !is.na(index.2))
    
    index.1 <- index.1[both.index]
    index.2 <- index.2[both.index]
    tmp.index <- cbind(index.1,index.2)
    W[[i]][tmp.index] <- Network[[i]][both.index,3]
    
    tmp.index <- cbind(index.2,index.1)
    W[[i]][tmp.index] <- Network[[i]][both.index,3]
    
  }
  
  K <- length(Network)
  
  gene.num <- length(Genes)
  Neighbor.index  <- list()
  Neighbor.weight <- list()
  Neighbor.count  <- matrix(0,nrow=K,ncol=gene.num)
  Neighbor.sumW  <- matrix(0,nrow=K,ncol=gene.num)
  
  # processing tissue-specific gene networks
  for(k in 1:K){
    print(sprintf('processing network %d',k))
    Neighbor.index[[k]]  <- list()
    Neighbor.weight[[k]] <- list()
    
    for(i in 1:length(Genes)){
      Neighbor.index[[k]][[i]]  <- which(W[[k]][i,]>0)
      Neighbor.weight[[k]][[i]] <- W[[k]][i,Neighbor.index[[k]][[i]]]
    }
    
    Neighbor.count[k,] <- apply(W[[k]]>0,1,sum)
    Neighbor.sumW[k,]  <- apply(W[[k]],1,sum)
  }
  
  # perform inference and output results
  pvalues <- GWAS.p[,'pvalue']
  gene.num <- length(pvalues)
  nets <- list(index=Neighbor.index, weight = Neighbor.weight, count = Neighbor.count, sumW = Neighbor.sumW)
  result <- SIGNET_infer(iter=iters,verbose=TRUE,pvals=pvalues,network=nets)
  result$data[,'Gene'] <- Genes
  result$network.name <- network.names
  return(result)
}
