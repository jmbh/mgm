# ******************************************************************************
# From https://github.com/jmbh/mgmDocumentation/blob/master/examples_tvmgm.R
# An example of applying mgm::tvmgm to Gaussian data of dimension 67x150
# ******************************************************************************

# jonashaslbeck@gmail.com; January 2019

# ----------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------
# ------------------- Code and Examples: Installation ------------------------------------------------
# ----------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------

library(mgm) # 1.2-5
library(qgraph)

# !!! Make sure to set the working directory to the path of the present R-file !!!

figDir <- "figures"
codeDir <- ""
fileDir <- "files"

# ----------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------
# ------------------- Fruit Fly Application: "Time-varying Mixed Graphical Models' -------------------
# ----------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------

head(fruitfly_data$data[,1:5])

# -------------------- Estimating Bandwidth parameter --------------------

set.seed(1)
# Note: This can take a while ~1h on a MacBook Pro
p <- ncol(fruitfly_data$data)
bw_tvmgm <- bwSelect(data = fruitfly_data$data,
                     type = rep('g', p),
                     level = rep(1, p),
                     bwSeq = c(0.1, 0.2, 0.3, 0.4),
                     bwFolds = 5,
                     bwFoldsize = 5,
                     modeltype = 'mgm', k = 2,
                     threshold = 'none', ruleReg = 'OR',
                     timepoints = fruitfly_data$timevector)

# saveRDS(bw_tvmgm, file = paste0(codeDir, "bw_tvmgm.RDS"))
bw_tvmgm <- readRDS(file = paste0(fileDir, "bw_tvmgm.RDS"))

round(bw_tvmgm$meanError, 3)
which.min(bw_tvmgm$meanError)


# -------------------- Estimating time-varying MGM --------------------

set.seed(1)
fit_tvmgm <- tvmgm(data = fruitfly_data$data,
                   type = rep("g", p),
                   level = rep(1, p),
                   timepoints = fruitfly_data$timevector,
                   estpoints = seq(0, 1, length=20),
                   k = 2,
                   bandwidth = 0.3,
                   threshold = "none",
                   ruleReg = "OR")

# saveRDS(fit_tvmgm, file = paste0(codeDir, 'fit_tvmgm.RDS'))
fit_tvmgm <- readRDS(file = paste0(fileDir, 'fit_tvmgm.RDS'))

# get wadj
wadj <- fit_tvmgm$pairwise$wadj
adj <- wadj
adj[adj!=0]<-1

# number of edges across estimation points
n_edges <- apply(adj, 3, sum)/2



# -------------------- Make Predictions from time-varying MGM --------------------

pred_tvmgm <- predict(object = fit_tvmgm,
                      data = fruitfly_data$data,
                      tvMethod = "weighted")


# -------------------- Visualizing Mixed Graphical Models --------------------

n <- nrow( fruitfly_data$data)
estpoints <- seq(0, n, length=20)

# Selected estimation points
E_select <- c(2, 6, 13)
round(E_select / 20 * 67, 2) # estimation points on true time scale

## z) Some work on the data

# normalized Time vector
tv <- fruitfly_data$timevector
tv <- tv - min(tv)
tv <- tv / max(tv)


# Color shade for Nodes
wDegree <- list()
for(i in 1:20)  wDegree[[i]] <- colSums(adj[,,i]) + 1; wDegree[[i]][wDegree[[i]]>9]<-9
n_color <- max(unlist(wDegree))
node_cols <- RColorBrewer::brewer.pal(n_color, 'Blues')


scale_size <- 1
pdf(paste0(figDir,'Fig_tvmgm_fruitfly_example.pdf'), width = 10*scale_size, height = 7*scale_size)

# a) Setup Layout
lom <- matrix(c(1,1,1,
                2,3,4), ncol=3, byrow = TRUE)
lo <- layout(lom, heights = c(1, 1))
# layout.show(lo)


# b) Top: Edges/number of observations used over time

# compute proportions of time
fruitfly_data$stages
csms <- c(0,cumsum(tv)[c(32, 42, 60, 67)])
csms <- csms / max(csms)

# Setup plot area
plot.new()
par(mar = c(3,3,2,3))
plot.window(xlim=c(-.02, 1), ylim=c(0,270))

# Set up rectangle-stages
library(scales)
col_stagesAlph <- alpha(RColorBrewer::brewer.pal(4, 'Set1'), .3)
col_stages <- RColorBrewer::brewer.pal(4, 'Set1')
for(i in 1:4) rect(csms[i], 0, csms[i+1], 250, col = col_stagesAlph[i])
dash_size <- 5
title(xlab = 'Estimation Points', cex.lab = 1.5, line=1, col.lab='blue')
segments(estpoints/67, -dash_size, estpoints/67, +dash_size, col='blue')

# Add time line
arrows(0, 265, 1, 265, length = .1)
mtext('Time                                ', 3, srt = 90, col='black')
mtext('                       / Measurements across time', 3, srt = 90, col='red')
# Add measurements
segments(tv, 265-dash_size, tv, 265+dash_size, col='red')


# add axis left: (number of edges)
y_ticks <- c(0, 50, 150, 200, 250)
tick_horiz <- .01
for(i in 1:5) {
  segments(0, y_ticks[i] , -tick_horiz,  y_ticks[i])
  text(-tick_horiz-.015, y_ticks[i], y_ticks[i])
}
mtext('Number of Edges', 2, srt = 90)


# add axis right: Proportion used sample size
y_ticks_labels <- c(0, 25, 50, 67)
y_ticks <- y_ticks_labels * 250 / 67
tick_horiz <- .01
for(i in 1:5) {
  segments(1, y_ticks[i] , 1+tick_horiz,  y_ticks[i])
  text(1+tick_horiz+.015, y_ticks[i], y_ticks_labels[i])
}
mtext('Local N', 4, srt = 90, line=1.3, col='black')


# Add Description of Stages
shift_left <- .015
cex_stages <- 1.5
text(csms[2]-shift_left, 80, 'Embryo', srt=90, cex = cex_stages, col = col_stages[1])
text(csms[3]-shift_left, 80, 'Larva', srt=90, cex = cex_stages, col = col_stages[2])
text(csms[4]-shift_left, 190, 'Pupa', srt=90, cex = cex_stages, col = col_stages[3])
text(csms[5]-shift_left, 190, 'Adult', srt=90, cex = cex_stages, col = col_stages[4])

# Plot data: number of edges
points(estpoints/67, n_edges, pch=20)
lines(estpoints/67, n_edges, pch=20)

# Plot data: proportion used Sample size
rel_y <- 250/67
points(estpoints/67, fit_tvmgm$Ne*rel_y, pch=21)
lines(estpoints/67, fit_tvmgm$Ne*rel_y, pch=20, lty=2)

# legend
legend(.6, 210, c('Number of Edges', 'Local n'), lty=1:2, pch=20:21, bty='n', cex=1.4)


# c) Bottom: Graphs at three different time points

# preliminary plotting
for(i in E_select) qgraph(adj[,,i],
                          layout='spring',
                          repulsion=1.05,
                          labels = F,
                          # edge.color = fit_tvmgm_nTH$pairwise$edgecolor[,,i],
                          color = node_cols[wDegree[[i]]],
                          mar = c(6, 6, 6, 6))


dev.off()
