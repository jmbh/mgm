
mgmfit <- function(
  data, #data matrix, col=variables
  type, #data type for col 1:ncol; c=categorical, g=gaussian, p=poisson, e=exponential
  lev, #number of categories of categorical variables, continuous variables have level=1
  lambda.sel="EBIC", #method for penalization parameter (lambda) -selection 
  folds=10, #folds in case CV is used for lambda selection
  gam=.25, #tuning parameter for EBIC, in case EBIC is used for lambda selection
  d=1, #maximal degree of the true graph
  rule.reg="AND", #parameter-aggregation of categorical variables
  rule.cat="OR",  #either "OR"- conditional independence; or matrix that specifies costumized rules
  pbar = TRUE # shows a progress bar if TRUE
)

{
  
  
  dev_glmnet <- function (object, ...) 
  {
    dev = object$dev
    nulldev = object$nulldev
    (1 - dev) * nulldev
  }
  
  # step 1: sanity checks & info from data
  stopifnot(ncol(data)==length(type)) # type vector has to match data
  stopifnot(ncol(data)==length(lev)) # level vector has to match data
  if(sum(apply(data, 2, function(x) class(x)!="numeric"))>0) stop("Only numeric values permitted!")
  if(sum(apply(cbind(data[,type=="p"],rep(1,nrow(data))), 2, function(x) sum(x-round(x))))>0) stop("Only integers permitted for Poisson random variables!")
  
  n <- nrow(data)
  nNode <- ncol(data)
  
  # step 2: prepare data
  #data[,type!="c" & type!="p"] <- scale(data[,type!="c" & type!="p"]) #standardize continuous variables
  data <- as.data.frame(data) #necessary for formula input
  colnames(data) <- paste("V",1:nNode, sep="") #necessary for formula input
  
  
  #compare entered and empirical levels  
  emp_lev <- numeric(nNode) + 1
  for(z in 1:nNode)
  {
    if(type[z]=="c") {
      emp_lev[z] <- length(unique(data[,z]))
      data[,z] <- as.factor(data[,z]) #turn categoricals into factors (for formula function)  
    }
  }
  
  
  # warning: entered lev = emp lev ?
  if(sum(emp_lev[type=="c"]!=lev[type=="c"])) warning("Entered levels are not equal to empirical levels. Empirical levels are used.")
  
  
  #indexing-dummy that helps us to put parameters involving categorical variables in the right place
  dummy_par.sort <- logical() # this will tell us, in which places of each row we fill parameters; this is because whenn gauss <- cat; we get m-1 parameter
  dummy_levels <- numeric()
  for(i in 1:nNode)
  {
    lev.e <- 1
    tar.e <- TRUE
    if(emp_lev[i]>1)
    {
      tar.e <- c(FALSE,rep(TRUE,emp_lev[i]-1))
      lev.e <- emp_lev[i]-1
    }
    dummy_par.sort <- c(dummy_par.sort, tar.e)
    dummy_levels <- c(dummy_levels, lev.e)
  }
  
  dummy_matrix.int <- cbind(dummy_levels, 1:length(dummy_levels))
  dummy.ind <- unlist(apply(dummy_matrix.int,1,function(x) { rep(x[2],x[1])}))
  
  dummy_matrix <- cbind(emp_lev, 1:length(emp_lev))
  ind <- as.numeric(unlist(apply(dummy_matrix,1,function(x) { rep(x[2],x[1])})))
  
  
  # step 3: create storage for parameters
  model.par.matrix <- matrix(0, sum(emp_lev), sum(emp_lev))
  m_lambdas <- matrix(0,nNode,2) #storing lambda threshold and actual lambda
  
  #progress bar
  if(pbar==TRUE) {
    pb <- txtProgressBar(min = 0, max=nNode, initial=0, char="-", style = 3)
  }
  
  # step 4: estimation
  for(v in seq_len(nNode))
  {
    # step 4.1: compute design matrix (adding interactions as a function of d)
    if(d>(nNode-1)) {
      stop("Order of interactions can be maximal the number of predictors!")
    } else if (d==1){ form <- as.formula(paste(colnames(data)[v],"~ (.)"))
    } else { form <- as.formula(paste(colnames(data)[v],"~ (.)^",d)) }
    
    
    X <- model.matrix(form, data=data)[,-1]
    
    #define link function
    if(type[v]=="c") {
      fam <- "multinomial"
    } else if(type[v]=="g" | type[v]=="e") { #should be inverse link for "e", but currently not avail. for glmnet
      fam <- "gaussian"
    } else if(type[v]=="p") {
      fam <- "poisson"
    }
    
    
    # step 4.2: select alpha & call glmnet
    
    #lambda selection with EBIC
    if(lambda.sel=="EBIC") {
      
      fit <- glmnet(X, data[,v], family=fam, alpha=1)
      
      #glmnet doesnt give us the pseudo LL, therefore we have to calculate it
      
      #calculate LL_Null depending on cat/cont
      if(type[v]=="g") {
        
        mean.i <- coef(fit, s=1)[1] #mean given by intercept model
        LL_null <- sum(dnorm(data[,v],mean.i,1, log=TRUE))
        
      } else if(type[v]=="e") {
        
        mean.i <- coef(fit, s=1)[1] #mean given by intercept model
        LL_null <- sum(dnorm(data[,v],mean.i,1, log=TRUE))
        
      } else if(type[v]=="p") {
        
        mean.i <- coef(fit, s=1)[1] #mean given by intercept model
        LL_null <- sum(dpois(data[,v],exp(mean.i), log=TRUE))
        
      } else if(type[v]=="c") {
        
        n_cats <- emp_lev[v]
        
        #dummy matrices to compute LL
        ind_dum <- cbind(matrix(0,n,n_cats), as.numeric(data[,v]))
        ind_mat <- t(apply(ind_dum, 1, function(x) {x[x[4]] <- 1; return(x[1:n_cats])} ))
        int_mat <- matrix(0,n,n_cats)
        for(ca in 1:n_cats) { int_mat[,ca] <- coef(fit, s=1)[[ca]][1] }
        
        LL_null <- 1/n * (sum(ind_mat * int_mat) - n*log(sum(exp(int_mat[1,]))) ) #LL multinomial from glmnet paper
      }
      
      # calc LL_sat
      LL_sat <- 1/2 * fit$nulldev + LL_null
      
      # calc LL for all lambdas
      dev <- dev_glmnet(fit)
      
      LL <- - 1/2 * dev + LL_sat
      
      n_lambdas <- length(fit$lambda)
      
      # calculation of nonzero neighborhoods
      if(type[v]!="c") { #continuous case
        coefs_bin <- as.matrix(coef(fit)[-1,]) != 0 #nonzero?
        n_neighbors <- colSums(coefs_bin)
      }
      if(type[v]=="c"){ #categorical case
        m_neighbors <- matrix(0,ncol=n_lambdas, nrow=n_cats)
        coefs_bin <- vector("list", length=n_cats)
        for(ca in 1:n_cats){
          coefs_bin[[ca]] <- as.matrix(coef(fit)[[ca]][-1,]) != 0 #nonzero?
        }
        n_neighbors <- colSums(Reduce('+', coefs_bin)!=0) #rule: a predictor has a nonzero parameter with 1 category of the y, then we have a neighborhood relation
      }
      
      # calc all EBICs
      EBIC_lambda <- -2*LL + n_neighbors * log(n) + 2*gam*n_neighbors*log(ncol(X)) 
      lambda_select <- fit$lambda[which(EBIC_lambda==min(EBIC_lambda))]
      coefs <- coef(fit, s=lambda_select) #lambda value with highest EBIC
      
      
      # lambda selection with CV
    } else {
      
      fit <- cv.glmnet(X, data[,v], family=fam, alpha=1, nfolds=folds, type.measure = "deviance")
      lambda_select <-  fit$lambda.min
      coefs <- coef(fit, s=lambda_select)
      
    } # end of estimation; 
    
    
    #list to matrix; cut out intercepts
    coefsm <- matrix(do.call(rbind,lapply(coefs, as.numeric)),nrow=emp_lev[v])[,-1]
    
    # step 4.3: save lambda + save & apply tau threshold
    m_lambdas[v,2] <- bound <- sqrt(d) * sqrt(sum(coefsm^2)) * sqrt(log(nNode)/n)
    m_lambdas[v,1] <- lambda_select
    coefsm[abs(coefsm)<bound]<-0 # apply tau threshold
    
    
    # step 4.4: write into model.par.matrix
    
    #select corresponding row in model par matrix & fill in
    #get correct dummy
    dummy_par.sort.v <- dummy_par.sort[ind!=v]
    
    #select corresponding row(s) in model par matrix & fill in
    # continuous
    if(emp_lev[v]==1) {
      exp.n.c <- length(model.par.matrix[ind==v,ind!=v][dummy_par.sort.v]) #number of coefficients
      model.par.matrix[ind==v,ind!=v][dummy_par.sort.v] <- coefsm[1:(exp.n.c)]
      
    } else { # categorical
      
      for(L in 1:emp_lev[v])
      {
        exp.n.c <- length(model.par.matrix[ind==v,ind!=v][,dummy_par.sort.v][L,])
        model.par.matrix[ind==v,ind!=v][,dummy_par.sort.v][L,] <- coefsm[L,1:(exp.n.c)]
      }
    }
    
    #progress bar
    if(pbar==TRUE) {
      setTxtProgressBar(pb, v)
    }
    
    
  } # end variable-loop
  
  # step 5: derivates on model parameter matrix
  
  
  # 5.1: aggregate on within categories
  
  f_agg_cats <- function(model.par.matrix, rule) {
    
    #select only colums where paramater are actually estimated (glmnet estimates k-1 not k parameters)
    m.p.m <- model.par.matrix[,dummy_par.sort]
    
    # averaging over  columns
    m.p.m.1 <-  t(apply(m.p.m, 1, function(x) {
      
      out <- numeric(0)
      for(i in 1:nNode)
      {
        out.n <- mean(abs(x[dummy.ind==i])) #without abs, this keeps the sign; but because of the glmnet parameterization in categoricals it burries nonzero coefficients in the binary case
        if(rule=="AND") {
          out.n <- out.n * (sum(x[dummy.ind==i]==0)<1) #the second term = 0 when not all coefficients are nonzero
        }
        out <- rbind(out, out.n)
      }
      out <- matrix(out, nrow=1)
    }))
    
    # averaging over rows
    m.p.m.2 <-  apply(m.p.m.1, 2, function(x) {
      out <- numeric()
      for(i in 1:nNode)
      {
        
        out.n <- mean(abs(x[ind==i]))
        if(rule=="AND") {
          out.n <- out.n * (sum(x[ind==i]==0)<1) #the second term = 0 when not all coefficients are nonzero
        }
        out <- rbind(out, out.n)
      }
      out <- matrix(out, ncol=1)
    })
    
  } #end of function
  
  #costumized rule
  if(rule.cat!="OR")
  {
    #check on matrix
    if(sum(dim(rule.cat) == dim(model.par.matrix))!=2) stop("rule.cat must have the same dimension as the parameter matrix!")
    
    #apply rule
    ind_costcat <- ((model.par.matrix + rule.cat)!=1)*1
    diag(ind_costcat)<-0
    ind_costcat_agg <- f_agg_cats(ind_costcat, "AND")
    m.p.m.2 <- f_agg_cats(model.par.matrix, "OR") * ind_costcat_agg
    
  #standard rule: OR which leads to edges that indicate conditional independence
  } else {
    m.p.m.2 <- f_agg_cats(model.par.matrix, "OR")  
  }
  
  
  ### 5.3: aggregate across two regressions
  
  if(rule.reg=="AND") {
    m.p.m.2_nonzero <- m.p.m.2!=0
    m.p.m.2_nonzero <- m.p.m.2_nonzero * t(m.p.m.2_nonzero)
    m.p.m.2 <- m.p.m.2 * m.p.m.2_nonzero
  }
  
  #make matrices symmetric (taking the average)
  wadj <- (m.p.m.2 + t(m.p.m.2))/2 #adjacency matrix
  mpar.matrix.sym <- (model.par.matrix+t(model.par.matrix)) / 2
  
  #create list mapping: parameters <-> variables as input for qgraph "group"-argument
  indvar.map <- indvar.map.label <-  vector("list", length=nNode)
  ind_map <- 1
  for(m in 1:nNode){
    
    #create indices list for qgraph
    indvar.map[[m]] <- ind_map:(ind_map+lev[m]-1)
    
    #create labels for qgraph
    if(lev[m]==1)
    {
      indvar.map.label[[m]] <- m
    } else {
      indvar.map.label[[m]] <- paste(m, 1:lev[m], sep = ".")
    }
    ind_map <- ind_map + lev[m]
  }
  
  indvar.map.label_all <- do.call(c, indvar.map.label)
  
  #dichotomize
  adj <- (wadj!=0)*1
  
  # step 6: output
  output_list <- list("adj"=adj, "wadj"=wadj, "wpar.matrix" = model.par.matrix, 
                      "wpar.matrix.sym"=mpar.matrix.sym, "indvar.map"=indvar.map, 
                      "indvar.map.labels"=indvar.map.label_all, "lambda"=m_lambdas[,1])
  

class(output_list) <- "mgm"
  
  return(output_list)
  
} 











