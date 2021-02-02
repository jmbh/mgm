

getSign <- function(l_w_ind,
                    l_w_par,
                    type,
                    set_signdefined,
                    overparameterize,
                    ord)
  
  
  
{
  
  
  # ---------- Compute Aux Variables ----------
  
  pair <- l_w_ind[[1]]
  
  outlist <- list("voteSign" = NULL,
                  "Signs" = NULL)
  
  # ---------- A) For overparameterize = TRUE ----------
  
  if(overparameterize) {
    
    # We have to take care of 3 cases: con-con, cat-cat, con-cat
    
    # I) ----- continuous-continuous -----
    
    if(all(type[pair]!='c')) {
      v_sign <- rep(NA, ord+1)
      for(u in 1:(ord+1)) v_sign[u] <-  sign(l_w_par[[u]])
      int_sign <- sign(mean(v_sign)) # Majority vote; if equal, we get 0 = undefined
      
      outlist$voteSign <- int_sign
      outlist$Signs <- v_sign
      
    } # end if: I)
    
    
    # II) ----- binary-binary -----
    
    if(all(type[pair]=='c')) {
      
      if(ord == 1) {
        ## only need to check parameters for one category, because of symmetry
        # regression1
        sign1 <- 0 #  default in case regression1 led to a zero estimate
        if(l_w_par[[1]][[1]][1] != 0) sign1 <- sign(l_w_par[[1]][[1]][1])
        if(l_w_par[[1]][[1]][2] != 0) sign1 <- - sign(l_w_par[[1]][[1]][2])
        # regression2
        sign2 <- 0 #  default in case regression1 led to a zero estimate
        if(l_w_par[[2]][[1]][1] != 0) sign2 <- sign(l_w_par[[2]][[1]][1])
        if(l_w_par[[2]][[1]][2] != 0) sign2 <- - sign(l_w_par[[2]][[1]][2])
        int_sign <- sign(mean(c(sign1, sign2))) # Majority vote
        
        outlist$voteSign <- int_sign
        outlist$Signs <- c(sign1, sign2)
        
      } else {
        int_sign <- 0 # no sign defined for interactions of order k>2
        
        outlist$Signs <- NA
        outlist$voteSign <- int_sign
        
      }
      
      
    } # end if: II)
    
    
    # III) ----- continuous-binary -----
    
    if(any(type[pair] %in% 'c') & any(type[pair] %in% c('p', 'g')) ) {
      
      if(ord == 1) {
        sign1 <- sign2 <- 0 # set default in case one direction has zero estimates
        # need to know which list entry in l_w_par corresponds to which regression: cont <- binary or cont -> binary; I do that by the fixed dimensionality of the parameter vector/matrix
        if(is.null(dim(l_w_par[[1]])))  { #is.null -> continuous, else binary
          sign1 <- sign(as.numeric(l_w_par[[1]][[2]]))
        } else {
          if(l_w_par[[1]][1, 1] != 0) sign1 <- - sign(l_w_par[[1]][1,1]) # positive value for state 1 means negative 'pairwise relationship'
          if(l_w_par[[1]][2, 1] != 0) sign1 <- sign(l_w_par[[1]][2,1])
        }
        if(is.null(dim(l_w_par[[2]])))  {
          sign2 <- sign(l_w_par[[2]][2])
        } else {
          if(l_w_par[[2]][1,1] != 0) sign2 <- - sign(l_w_par[[2]][1,1]) # positive value for state 1 means negative 'pairwise relationship'
          if(l_w_par[[2]][2,1] != 0) sign2 <- sign(l_w_par[[2]][2,1])
        }
        int_sign <- sign(mean(c(sign1,sign2))) # Majority vote
        
        outlist$voteSign <- int_sign
        outlist$Signs <- c(sign1, sign2)
        
      } else {
        int_sign <- 0 # no sign defined for interactions of order k>2
        
        outlist$Signs <- NA
        outlist$voteSign <- int_sign
      }
      
    } # end if: III)
    
    
  } # end if: overparameterize?
  
  # ---------- B) For overparameterize = FALSE ----------

  if(!overparameterize) {
    
    
    # I) ----- continuous-continuous -----
    
    if(all(type[pair]!='c')) {
      v_sign <- rep(NA, ord+1)
      for(u in 1:(ord+1)) v_sign[u] <-  sign(l_w_par[[u]])
      int_sign <- sign(mean(v_sign)) # Majority vote; if equal, we get 0 = undefined
      
      outlist$voteSign <- int_sign
      outlist$Signs <- v_sign
      
    } # end if: I)
    
    
    # II) ----- binary-binary -----
    
    if(all(type[pair]=='c')) {
      if(ord == 1) {
        sign1 <- sign(as.numeric(l_w_par[[1]][[2]]))
        sign2 <- sign(as.numeric(l_w_par[[2]][[2]]))
        int_sign <- sign(mean(c(sign1, sign2)))
        
        outlist$voteSign <- int_sign
        outlist$Signs <- c(sign1, sign2)
        
      } else {
        int_sign <- 0 # no sign defined for interactions of order k>2
        
        outlist$Signs <- NA
        outlist$voteSign <- int_sign
      }
    } # end if: II)
    
    
    # III) ----- continuous-binary -----
    
    if(any(type[pair] %in% 'c') & any(type[pair] %in% c('p', 'g')) ) {
      
      if(ord == 1) {
        
        if(length(l_w_par[[1]]) == 1) {
          sign1 <- sign(as.numeric(l_w_par[[1]])) # if interaction A is cont <- binary
        } else {
          sign1 <- sign(as.numeric(l_w_par[[1]][[2]])) # if interaction A is binary <- cont
        }
        
        # same for second interaction
        if(length(l_w_par[[2]]) == 1) {
          sign2 <- sign(as.numeric(l_w_par[[2]]))
        } else {
          sign2 <- sign(as.numeric(l_w_par[[2]][[2]]))
        }
        
        int_sign <- sign(mean(c(sign1, sign2)))
        
        outlist$voteSign <- int_sign
        outlist$Signs <- c(sign1, sign2)
        
      } else {
        int_sign <- 0 # no sign defined for interactions of order k>2
        outlist$Signs <- NA
        outlist$voteSign <- int_sign
      }
      
    } # end if: III)
    
    
  } # end if: overparameterize?
  
  # ---------- Return Sign ----------
  
  return(outlist)
  
} # eoF
