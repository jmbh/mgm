# jonashalsbeck@gmail.com; August 2019


# Function to plot LHS legend

plotLegend <- function(labels, 
                       TM, 
                       margins, 
                       ylim, 
                       plot_y, 
                       cex.label, 
                       p = p) {
  
  
  # Generate label vector
  if(is.null(labels)) {
    label_vec <- paste0(TM[, 1], " - ", TM[, 2])
  } else {
    tar_mat_label <- TM[ ,1:2]
    tar_mat_label <- apply(tar_mat_label, 1:2, as.character)
    for(i in 1:p) tar_mat_label[tar_mat_label == i] <- labels[i]
    label_vec <- paste0(tar_mat_label[, 1], " - ", tar_mat_label[, 2])
  }
  
  # Plot  
  par(mar=margins)
  plot.new()
  plot.window(xlim = c(-1, 1), ylim = ylim)
  text(0, plot_y, label_vec, cex = cex.label)
  
} # eoF


# Function to plot data

plotData <- function(TM, 
                     axis.ticks,
                     ylim,
                     plot_y,
                     lwd.qtl, 
                     cex.bg, 
                     cex.mean, 
                     margins,
                     bgcol = "white") {
  
  # Some settings  
  xlim <- range(axis.ticks)
  
  # Setup canvas
  par(mar=margins)
  plot.new()
  plot.window(xlim = xlim, ylim = ylim)
  
  # background color 
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = bgcol, border = FALSE)
  # rect(xleft = xlim[1], ybottom = ylim[1], xright = xlim[2], ytop = ylim[2], col = bgcol, border = FALSE)
  
  if(is.null(axis.ticks)) axis.ticks <- round(seq(xlim[1], xlim[2], length = 5), 2)
  axis(3, axis.ticks, lwd=0)
  
  abline(h = plot_y, col = "grey")
  if(0 %in% axis.ticks) abline(v = 0, lty=2, col = "black") # zero line
  
  # Plot quantiles
  segments(x0 = TM[, 4], 
           y0 = plot_y, 
           x1 = TM[, 5],
           y1 = plot_y, lwd = lwd.qtl)
  
  # Plot prop>0
  points(TM[, 3], plot_y, pch=19, col=bgcol, cex = cex.bg)
  points(TM[, 3], plot_y, pch=21, col="black", cex = cex.bg)
  text(TM[, 3], plot_y, TM[, 6], cex = cex.mean)
  
  
} # eoF



plotRes <- function(object,
                    quantiles = c(.05, .95),
                    labels = NULL,
                    decreasing = TRUE, 
                    cut = NULL,
                    cex.label = .75,
                    lwd.qtl = 2, 
                    cex.mean = .55, 
                    cex.bg = 2.7, 
                    axis.ticks = c(-.5, -.25, 0, .25, .5, .75, 1), 
                    axis.ticks.mod = NULL,
                    layout.width.labels = .2, 
                    layout.gap.pw.mod = .15, 
                    table = FALSE)

{
  
  
  # ----------------------------------------------------------------------
  # ---------- Aux variables ---------------------------------------------
  # ----------------------------------------------------------------------
  
  # Fill in defaults
  
  # Get vars
  mod <- object$call$object$call$moderators
  type <- object$call$object$call$type
  p <- length(type)
  nB <- object$call$nB
  
  # ----------------------------------------------------------------------
  # ---------- Input Checks ----------------------------------------------
  # ----------------------------------------------------------------------
  
  # Input checks: basic
  if(!("core" %in% class(object))) stop("plotRes() currently only supports resampled mgm() objects.")
  if(!("resample" %in% class(object))) stop("PplotRes() only takes resample objects as input (see ?resample).")
  # if(missing(quantiles)) stop("No quantiles specified.")
  if(object$call$object$call$k > 3) stop("Currently only implemented with mgms with k <= 3.")
  
  # Check whether there is only a single moderator
  if(class(mod)[1] == "numeric") if(length(mod) > 1) stop("plotRes() is currently only implemented for MNNs with a single moderator.")
  
  if(class(mod)[1] == "matrix") {
    
    # find whether one variable is in each 3-way interaction; if yes, that's the one we treat as the moderator
    n_mod <- nrow(mod)
    m_check <- matrix(NA, n_mod, p)
    
    for(i in 1:p) for(m in 1:n_mod) m_check[m, i] <- i %in% mod[m, ] 
    v_check <- apply(m_check, 2, function(x) all(x == TRUE))
    
    mod_fix <- which(v_check)
    if(is.null(mod_fix)) stop("plotRes() is currently only implemented for MNNs with a single moderator.")
    
    # pretend that we have all possible moderation effect of mod; this will result in unmodeled moderation effects to be set to zero (as they should be) 
    mod <- mod_fix
  }
  
  
  # ----------------------------------------------------------------------
  # ---------- A.1) mgm objects (unmoderated) ----------------------------
  # ----------------------------------------------------------------------
  
  if(is.null(mod)) {
    
    # ---------- Preprocessing ----------  
    
    # Get basic info
    dims <- dim(object$bootParameters)
    p <- dims[1]
    nB <- dims[3]
    n_pars <- p*(p-1) / 2
    
    # Collapse into edge x property matrix
    tar_mat <- matrix(NA, nrow=n_pars, ncol = 6)
    colnames(tar_mat) <- c("Variable A", "Variable B", "Mean", "qtl_low", "qtl_high", "propLtZ")
    
    counter <- 1
    for(row in 1:p) {
      for(col in row:p) {
        if(row!=col){
          
          # Variable ids
          tar_mat[counter, 1] <- row
          tar_mat[counter, 2] <- col
          
          # Quantiles
          qtls <- quantile(object$bootParameters[row, col, ], probs = quantiles)
          tar_mat[counter, 3] <- mean(object$bootParameters[row, col, ])
          tar_mat[counter, 4] <- qtls[1]
          tar_mat[counter, 5] <- qtls[2]
          tar_mat[counter, 6] <- mean(abs(object$bootParameters[row, col, ]) > 0) # proportion estimates > 0
          
          # update counter
          counter <- counter + 1
        }
      }
    }
    
    tar_mat <- round(tar_mat, 2)
    
    # Order
    tar_mat <- tar_mat[order(tar_mat[,3], decreasing = decreasing), ]
    
    # Subset (cut)
    if(is.null(cut)) {
      TM <- tar_mat
    } else {
      TM <- tar_mat[cut, ]
    }
    
    if(!table) {
    # ---------- Plotting ----------  
    
    # Compute aux variables for plotting
    n_rows <- nrow(TM)
    
    ylim <- c(0, 1)
    plot_y <- seq(ylim[2], ylim[1], length = n_rows)
    margins <- c(0, .5, 3, .5)
    
    # ----- Setup layout ----
    
    lmat <- matrix(1:2, nrow=1)
    lo <- layout(lmat, widths = c(layout.width.labels, 1))
    
    # ----- Part A: Legend ----
    
    plotLegend(labels = labels, 
               TM = TM, 
               margins = margins, 
               ylim = ylim, 
               plot_y = plot_y, 
               cex.label = cex.label, p=p)
    
    
    # ----- Part B: Data ----
    
    plotData(TM = TM, 
             axis.ticks = axis.ticks, 
             ylim = ylim, 
             plot_y = plot_y, 
             lwd.qtl = lwd.qtl, 
             cex.bg = cex.bg, 
             cex.mean = cex.mean, 
             margins=margins)
    
    
    # Return table insteaed
    } else {
      
      return(TM)

    }
    
  } # end if: moderation?
  
  
  
  
  # ----------------------------------------------------------------------
  # ---------- A.1) mgm objects (moderated; 1 moderator) -----------------
  # ----------------------------------------------------------------------
  
  if(!is.null(mod)) {
    
    # ---------- Preprocessing ----------  
    
    # ----- Get estimates out of model object -----
    
    m_ind_allpw <- t(combn(1:p, 2)) # list all possible 2-way interactions
    n_pw <- nrow(m_ind_allpw) # how many?
    m_pw <- m_mod <- matrix(0, nrow=n_pw, ncol=nB) # Create storage for nB pairwise and moderation effects
    
    # loop through bootstrapped objects
    for(b in 1:nB) {
      
      mod_b <- object$models[[b]]
      mod_b_pw <- matrix(mod_b$interactions$indicator[[1]], ncol=2)
      mod_b_m <- matrix(mod_b$interactions$indicator[[2]], ncol=3)
      
      # Loop pairwise & get pairwise
      
      if(nrow(mod_b_pw) > 0) for(i in 1:nrow(mod_b_pw)) {
        ind_pw <- which(apply(m_ind_allpw, 1, function(x) sum(x %in% mod_b_pw[i, ])) == 2)
        m_pw[ind_pw, b] <- mod_b$interactions$weightsAgg[[1]][[i]][1]
        
        # Add sign if applicable
        if(mod_b$interactions$signs[[1]][[i]] == -1) m_pw[ind_pw, b] <- m_pw[ind_pw, b] * -1
          
        # # add sign if applicable
        # if(all(type[mod_b_pw] == "g")) {
        #   v_est <- unlist(mod_b$interactions$weights[[1]][[i]])
        #   v_sign <- sign(v_est)
        #   sign_sel <- as.numeric(names(which.max(table(v_sign)))) # ugly!
        #   m_pw[ind_pw, b] <- m_pw[ind_pw, b] * sign_sel
        # }
        
      } # end: pw int
      
      # Loop moderation & get moderation
      if(nrow(mod_b_m) > 0) for(i in 1:nrow(mod_b_m)) {
        ind_mod <- which(apply(m_ind_allpw, 1, function(x) sum(x %in% mod_b_m[i, ][-which(mod_b_m[i, ]==mod)])) == 2)
        m_mod[ind_mod, b] <- mod_b$interactions$weightsAgg[[2]][[i]][1] #* mod_b$interactions$signs[[2]][[i]] # add sign
        
        # Add sign if applicable
        if(mod_b$interactions$signs[[2]][[i]] == -1) m_mod[ind_mod, b]  <- m_mod[ind_mod, b]  * -1
        
        # # add sign if applicable
        # if(all(type[mod_b_m] == "g")) {
        #   v_est <- unlist(mod_b$interactions$weights[[2]][[i]])
        #   v_sign <- sign(v_est)
        #   sign_sel <- as.numeric(names(which.max(table(v_sign)))) # ugly!
        #   m_mod[ind_mod, b] <- m_mod[ind_mod, b] * sign_sel
        # }
        
      } # end: mod 
      
    } # end: bootstrap it
    
    
    # ----- Compute quantile / median / prop>0 -----
    
    tar_mat_pw <- tar_mat_mod <- matrix(NA, nrow=n_pw, ncol = 6)
    colnames(tar_mat_pw) <- colnames(tar_mat_pw) <- c("Variable A", "Variable B", "Mean", "qtl_low", "qtl_high", "propLtZ")
    
    # Pairwise
    tar_mat_pw[, 1:2] <- m_ind_allpw
    tar_mat_pw[, 3] <- rowMeans(m_pw)
    tar_mat_pw[, 4:5] <- t(apply(m_pw, 1, function(x) quantile(x, probs = quantiles)))
    tar_mat_pw[, 6] <- apply(m_pw, 1, function(x) mean(x!=0) )
    tar_mat_pw <- round(tar_mat_pw, 2)
    
    # Moderation
    tar_mat_mod[, 1:2] <- m_ind_allpw
    tar_mat_mod[, 3] <- rowMeans(m_mod)
    tar_mat_mod[, 4:5] <- t(apply(m_mod, 1, function(x) quantile(x, probs = quantiles)))
    tar_mat_mod[, 6] <- apply(m_mod, 1, function(x) mean(x!=0) )
    tar_mat_mod <- round(tar_mat_mod, 2)
    
    # Order
    ord <- order(tar_mat_pw[,3], decreasing = decreasing)
    tar_mat_pw <- tar_mat_pw[ord, ]
    tar_mat_mod <- tar_mat_mod[ord, ]
    
    # Subset (cut argument)
    if(is.null(cut)) {
      TM_pw <- tar_mat_pw
      TM_mod <- tar_mat_mod
    } else {
      TM_pw <- tar_mat_pw[cut, ]
      TM_mod <- tar_mat_mod[cut, ]
    }
    
    if(!table) {
      
    # ---------- Plotting ---------- 
    
    n_rows <- nrow(TM_pw)
    ylim <- c(0, 1)
    plot_y <- seq(ylim[2], ylim[1], length = n_rows)
    margins <- c(0, .5, 3, .5)
    
    # ----- Setup layout -----
    
    lmat <- rbind(c(1,2,7,3), 
                  c(4,5,8,6))
    
    lo <- layout(lmat, 
                 widths = c(layout.width.labels, 1, layout.gap.pw.mod, 1), 
                 heights = c(.07, .9))
    
    # ----- 0) Plot Top Legend -----
    
    par(mar=rep(0, 4))
    
    plot.new() # fill 1
    
    plot.new() 
    plot.window(xlim=c(0, 1), ylim=c(0,1))
    text(0.5, 0.5, "Pairwise effects", cex=1.5, adj = 0.5)
    
    plot.new() 
    plot.window(xlim=c(0, 1), ylim=c(0,1))
    text(0.5, 0.5, "Moderation effects", cex=1.5, adj = 0.5)
    
    
    
    # ----- A) Plot LHS Legend -----
    
    plotLegend(labels = labels, 
               TM = TM_pw, 
               margins = margins, 
               ylim = ylim, 
               plot_y = plot_y, 
               cex.label = cex.label, p=p)
    
    # ----- B) Plot Pairwise effects -----
    
    
    plotData(TM = TM_pw, 
             axis.ticks = axis.ticks, 
             ylim = ylim, 
             plot_y = plot_y, 
             lwd.qtl = lwd.qtl, 
             cex.bg = cex.bg, 
             cex.mean = cex.mean, 
             margins = margins)
    
    # ----- C) Plot Moderation effects -----
    
    if(is.null(axis.ticks.mod)) axis.ticks.mod <- axis.ticks
    
    plotData(TM = TM_mod, 
             axis.ticks = axis.ticks.mod, 
             ylim = ylim, 
             plot_y = plot_y, 
             lwd.qtl = lwd.qtl, 
             cex.bg = cex.bg, 
             cex.mean = cex.mean, 
             bgcol = "white", 
             margins = margins)
    
    # Return table instead
    } else {
      
      # make table
      colnames(TM_mod)[3:6] <- c("Mod_Mean", "Mod_qtl_low", "Mod_qtl_high", "Mod_propLtZ") 
      out_table <- cbind(TM_pw, TM_mod)
      return(out_table)
      
    }
    
  } # end if: moderation?
  
  # ----------------------------------------------------------------------
  # ---------- B) mvar objects -------------------------------------------
  # ----------------------------------------------------------------------
  
  # To do ...
  
  
  
} # eoF

