
# -------------------------------------------------------------------------
# --------------- Loading packages & Data ---------------------------------
# -------------------------------------------------------------------------
# install.packages('huge')
# install_github("jmbh/mgm")

# install.packages('moments')
# install.packages('outliers')

# library(devtools)
# # library(mgm)
# library(qgraph)
# library('huge')
# library('moments')
# library(outliers)


# ---------- Continuous Moderation ----------


data_NCT <- readRDS(file = paste0("data_NCT.RDS"))
#data_NCT$names <- c("CRP", "Sad Mood", "Anhedonia", "Sleep", "Fatigue", "Appetite", "Motor", "Concen", "Feels Bas", "Death")

# Fit initial model


moderators1b <- rbind(c(1,2,3),
                      c(1,2,4),
                      c(1,2,5),
                      c(1,2,6),
                      c(1,2,7),
                      c(1,2,8),
                      c(1,2,9),
                      c(1,2,10),
                      c(1,3,4),
                      c(1,3,5),
                      c(1,3,6),
                      c(1,3,7),
                      c(1,3,8),
                      c(1,3,9),
                      c(1,3,10),
                      c(1,4,5),
                      c(1,4,6),
                      c(1,4,7),
                      c(1,4,8),
                      c(1,4,9),
                      c(1,4,10),
                      c(1,5,6),
                      c(1,5,7),
                      c(1,5,8),
                      c(1,5,9),
                      c(1,5,10),
                      c(1,6,7),
                      c(1,6,8),
                      c(1,6,9),
                      c(1,6,10),
                      c(1,7,8),
                      c(1,7,9),
                      c(1,7,10),
                      c(1,8,9),
                      c(1,8,10),
                      c(1,9,10))

head(data_NCT$data)

modelCRP_cov <- mgm(data = data_NCT$data, 
                    type = data_NCT$type, 
                    level = data_NCT$level, 
                    lambdaSel = "EBIC", 
                    binarySign = TRUE, 
                    lambdaGam = 0, 
                    lambdaSeq = 0, 
                    moderators = moderators1b)

modelCRP_cov$interactions$indicator

modelCRP_cov$nodemodels[[1]]$model



modelCRP_cov_full <- mgm(data = data_NCT$data, 
                         type = data_NCT$type, 
                         level = data_NCT$level, 
                         lambdaSel = "EBIC", 
                         binarySign = TRUE, 
                         lambdaGam = 0, 
                         lambdaSeq = 0, 
                         moderators = 1)

modelCRP_cov_full$interactions$indicator






