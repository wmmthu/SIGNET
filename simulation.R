library('Rcpp')
source('SIGNET.R')
source('utils.R')
library('PRROC')


f = '12_multiple_sclerosis.txt'
GWAS.p <- read.table(sprintf('../GWAS gene scores v1/%s',f),header=T,stringsAsFactors = F)
HLA.gene.index <- which(GWAS.p[,1]=='chr6' & GWAS.p[,3] > 25000000 & GWAS.p[,2] < 35000000)
if(length(HLA.gene.index) > 0){
  GWAS.p <- GWAS.p[-HLA.gene.index,]
}

Genes  <- GWAS.p[,'gene_symbol']
gene.num <- length(Genes)


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


simulated.network.num <- 32
K <- simulated.network.num
Network <- list()
W <- list()
# tissue specific regulatory network
net.filedir  <- '../regulatorycircuits/Network_compendium/Tissue-specific_regulatory_networks_FANTOM5-v1/32_high-level_networks/'
net.filelist <- list.files(net.filedir)

network.names <- c()
Neighbor.index  <- list()
Neighbor.weight <- list()
for(i in 1:simulated.network.num){
  print(i)
  tmp <- strsplit(net.filelist[i],'-')
  network.names <- c(network.names , net.filelist[i])
  Network[[i]] <- read.table(sprintf('%s/%s',net.filedir, net.filelist[i]), stringsAsFactors = F)
  
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

weighted <- T
if(!weighted){
for(i in 1:K){
  W[[i]] <- (W[[i]]>0) * 1
}
}


Neighbor.index  <- list()
Neighbor.weight <- list()
Neighbor.count  <- matrix(0,nrow=K,ncol=gene.num)
Neighbor.sumW  <- matrix(0,nrow=K,ncol=gene.num)


for(k in 1:K){
  print(k)
  Neighbor.index[[k]]  <- list()
  Neighbor.weight[[k]] <- list()

  for(i in 1:length(Genes)){
    Neighbor.index[[k]][[i]]  <- which(W[[k]][i,]>0)
    Neighbor.weight[[k]][[i]] <- W[[k]][i,Neighbor.index[[k]][[i]]]
  }
  
  Neighbor.count[k,] <- apply(W[[k]]>0,1,sum)
  Neighbor.sumW[k,]  <- apply(W[[k]],1,sum)
}


pvalues <- GWAS.p[,'pvalue']
gene.num <- length(pvalues)
nets <- list(index=Neighbor.index, weight = Neighbor.weight, count = Neighbor.count, sumW = Neighbor.sumW)

individual_nets <- list()
for(i in 1:K){
  tmp <- list(index=list(Neighbor.index[[i]]), weight = list(Neighbor.weight[[i]]), count = list(Neighbor.count[[i]]) )
  individual_nets[[i]] <- tmp
}



############################################## simulation
alpha1.list <- c(0.2)
gamma.list <- c(-4,-3,-2,-1)
beta.list <- c(1,2,5)
num.list <- c(1,2,3)
replicate.num <- 2
iters <- 10000

for(s.alpha1 in alpha1.list){
  for(s.gamma in gamma.list){
    
    AUCs <- matrix(0,nrow = 0, ncol = 14)
    colnames(AUCs) <- c('alpha1','beta','gamma','num','auc.pval','auc.infer44','auc.infer44.net','prc.pval','prc.infer44','prc.infer44.net','true.num','estimated.alpha0','estimated.alpha1','estimated.gamma')
    beta.truth <- matrix(0,nrow=0,ncol=10)
    beta.pip   <- matrix(0,nrow=0,ncol=10)
    
    for(s.beta in beta.list){
      for(s.num in num.list){
      for(rr in 1:replicate.num){
        
        start <- proc.time()
        
        simulated.beta  <- seq(0,0,length.out = K)
        simulated.beta[sample(1:K,s.num)] <- s.beta
        data <- simulation(simulated.alpha0 = 1, simulated.alpha1 = s.alpha1, simulated.gamma = s.gamma, simulated.beta = simulated.beta, pvals = pvalues, network=nets)
        result <- SIGNET_infer(iter = iters, verbose=F, pvals = data$pvals, network = nets)
        
        
        auc.pval <- roc.curve(1 - data$pvals[data$true.Z==1],1 - data$pvals[data$true.Z==0])$auc
        auc.infer <- roc.curve(1 - result$data$localFDR[data$true.Z==1],1 - result$data$localFDR[data$true.Z==0])$auc
        
        prc.pval <- pr.curve(1 - data$pvals[data$true.Z==1], 1 - data$pvals[data$true.Z==0])$auc.integral
        prc.infer <- pr.curve(1 - result$data$localFDR[data$true.Z==1],1 - result$data$localFDR[data$true.Z==0])$auc.integral
        
        
        if(sum(simulated.beta)>0){
        prc.infer.net <- pr.curve(result$network.included.prob[simulated.beta>0],result$network.included.prob[simulated.beta==0])$auc.integral
        auc.infer.net <- roc.curve(result$network.included.prob[simulated.beta>0],result$network.included.prob[simulated.beta==0])$auc
        }else{
          auc.infer.net <- 0
          prc.infer.net <- 0
        }
        
        estimates <- apply(result$para.list[(iters/2+1):iters,1:3],2,mean)
        
        rr <- c(s.alpha1,s.beta,s.gamma,s.num,auc.pval,auc.infer,auc.infer.net,prc.pval,prc.infer,prc.infer.net,sum(data$true.Z),estimates)
        AUCs <- rbind(AUCs,rr)
        print(rr)
        
        beta.truth <- rbind(beta.truth,simulated.beta)
        beta.pip <- rbind(beta.pip, result$network.included.prob)
      
        print(proc.time() - start)
      }
    }
    }
    
    outf <- sprintf('simulation.alpha1.%f.gamma.%f.bin',s.alpha1,s.gamma)
    simulation.result <- list(AUC=AUCs,beta.truth=beta.truth,beta.pip=beta.pip)
    save(simulation.result,file=outf)
  }
}


