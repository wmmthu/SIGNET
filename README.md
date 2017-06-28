# SIGNET

Code for "Simultaneous identification of disease genes and responsive tissues via Bayesian integration of multiple regulatory networks"

![Schematic diagram of SIGNET](https://github.com/wmmthu/SIGNET/raw/master/signet.jpg)

Requirements:
```
1) R package "Matrix" (https://cran.r-project.org/web/packages/Matrix/)
This package is used to support sparse matrix, which is used when loading gene networks (highly sparse).

2) R package "Rcpp" (https://cran.r-project.org/web/packages/Rcpp/index.html)
This package is used to support efficient implementation of Gibbs sampling

3) R package "BayesLogit" (https://cran.r-project.org/web/packages/BayesLogit/index.html)
This package provides an efficient sampler for Polya-Gamma random variable
```
Usage : 
```
(1) Compute gene-level p-values with PASCAL (http://www2.unil.ch/cbg/index.php?title=Pascal)

(2) Download tissue-specific gene networks (or other alternative context-specific gene networks), e.g. regulatorycircuits (http://regulatorycircuits.org/) and put these networks into one folder. 
Notice that only unzipped network files are allowed to put into this folder.

(3) Use the function "SIGNET" from SIGNET.R to perform inference.

command    : result <- SIGNET(p_value_file, network_dir, iters = 20000, remove_HLA=TRUE, edge_threshold=0.01)  
parameters :  
p_value_file   : gene-level p-values file generated in step (1)  
network_dir    : the folder path for gene networks obtained in step (2)  
iters          : the number of iterations for performing MCMC sampling  
remove_HLA     : whether or not to remove genes located at HLA region  
edge_threshold : the threshold for remove noisy edges in gene networks  

output : a list of result with components as  
result$network.name : the names for gene networks, obtained as file names for these networks  
result$data         : a dataframe , each row is a gene and three columns denote gene names ('Gene'), gene-level p-values ('pval'), gene-level local FDR ('localFDR')  
result$para.list    : a dataframe containing parameters of each MCMC sampling step  
result$network.included.prob : the posterior includsion probabilities for these gene networks  
```

If you have any questions, please contact me (Mengmeng Wu, wmm15@mails.tsinghua.edu.cn).
