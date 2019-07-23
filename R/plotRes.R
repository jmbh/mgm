

plotRes <- function(object,
                    quantiles,
                    labels = NULL,
                    decreasing = TRUE, 
                    cut = NULL,
                    cex.label = .75,
                    lwd.qtl = 2, 
                    cex.mean = .5, 
                    cex.bg = 3.5, 
                    axis.ticks = c(-.5, -.25, 0, .25, .5, .75, 1))
  
{
  
  # ---------- Input Checks ----------
  
  browser()
  
  if(!("core" %in% class(object))) stop("plotRes() currently only supports resampled mgm() objects.")
  if(!("resample" %in% class(object))) stop("PplotRes() only takes resample objects as input (see ?resample).")
  if(missing(quantiles)) stop("No quantiles specified.")
  if(object$call$k > 3) stop("Currently only implemented with mgms with k <= 3.")
  

  # ---------- A) mgm objects ----------
  
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
  
  
  # Order
  tar_mat <- tar_mat[order(tar_mat[,3], decreasing = decreasing), ]
  
  # Subset (cut)
  if(is.null(cut)) {
    TM <- tar_mat
  } else {
    TM <- tar_mat[cut, ]
  }
  
  # ---------- Plotting ----------  
  
  # Compute aux variables for plotting
  n_rows <- nrow(TM)

  ylim <- c(0, 1)
  plot_y <- seq(ylim[2], ylim[1], length = n_rows)
  mar <- c(0, .5, 3, .5)
  
  # ----- Setup layout ----
  
  lmat <- matrix(1:2, nrow=1)
  lo <- layout(lmat, widths = c(.2, 1))
  
  # ----- Part A: Legend ----
  
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
  par(mar=mar)
  plot.new()
  plot.window(xlim = c(-1, 1), ylim = ylim)
  text(0, plot_y, label_vec, cex = cex.label)
  
  
  # ----- Part B: Data ----
  
  # Some settings  
  xlim <- range(axis.ticks)
  
  # Setup canvas
  par(mar=mar)
  plot.new()
  plot.window(xlim = xlim, ylim = ylim)
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
  points(TM[, 3], plot_y, pch=20, col="white", cex = cex.bg)
  text(TM[, 3], plot_y, TM[, 6], cex = cex.mean)
  
  
  
  # ---------- B) mvar objects ----------
  
  
  
  
  
} # eoF

