

calcNeighbors <- function(fit, lambda, type, level, v) {

  n_cats <- level[v]

  if(type[v]!="c") { #continuous case
    coefs_bin <- as.matrix(coef(fit, s = lambda)[-1, ]) != 0 #nonzero?
    n_neighbors <- colSums(coefs_bin)
  }
  if(type[v]=="c"){ #categorical case
    m_neighbors <- matrix(0, ncol = length(fit$lambda), nrow = n_cats)
    coefs_bin <- vector("list", length=n_cats)
    for(ca in 1:n_cats){
      coefs_bin[[ca]] <- as.matrix(coef(fit, s = lambda)[[ca]][-1,]) != 0 #nonzero?
    }
    n_neighbors <- colSums(Reduce('+', coefs_bin)!=0) # rule: a predictor has a nonzero parameter with 1 category of the y, then we have a neighborhood relation
  }

  return(n_neighbors)

}







