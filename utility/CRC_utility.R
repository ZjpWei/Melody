## Sequence depth rarify function for BW model
.rarefy <- function (otu.tab, ss) {  
  set.seed(ss)
  depth = min(rowSums(otu.tab))  
  otu.tab <- as.matrix(otu.tab)  
  ind <- (rowSums(otu.tab) < depth)  
  sam.discard <- rownames(otu.tab)[ind]  
  otu.tab <- otu.tab[!ind, ]  
  rarefy <- function(x, depth) {    
    y <- sample(rep(1:length(x), x), depth)    
    y.tab <- table(y)    
    z <- numeric(length(x))    
    z[as.numeric(names(y.tab))] <- y.tab    
    z  
  }  
  otu.tab.rff <- t(apply(otu.tab, 1, rarefy, depth))  
  rownames(otu.tab.rff) <- rownames(otu.tab)  
  colnames(otu.tab.rff) <- colnames(otu.tab)  
  return(otu.tab.rff)
}

## Bootstrap function for simulation
boostrap <- function(vec){
  boot <- sample(vec, size = length(vec), replace = TRUE)
  return(sort(boot))
}
