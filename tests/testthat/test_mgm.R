

test_that("Pairwise mgm examples return expected estimates", {
  
  # --- Fit pairwise with EBIC ---
  fit_k2 <- mgm(data = autism_data$data,
                type = autism_data$type,
                level = autism_data$lev,
                lambdaSel = "EBIC", 
                lambdaGam = 0.25,
                k = 2, 
                pbar = FALSE, 
                signInfo = FALSE) 
  
  # check one pairwise interactiion
  expect_equal(round(fit_k2$pairwise$wadj[1, 4], 3), 0.152)
  # and one intercept
  expect_equal(round(fit_k2$intercepts[[1]][[1]], 3), 0.664)
  
  
  # --- Fit pairwise with CV ---
  set.seed(1)
  fit_k2 <- mgm(data = autism_data$data,
                type = autism_data$type,
                level = autism_data$lev,
                lambdaSel = "CV", 
                lambdaFolds = 10,
                k = 2, 
                pbar = FALSE, 
                signInfo = FALSE) 
  
  # check one pairwise interactiion
  expect_equal(round(fit_k2$pairwise$wadj[1, 4], 3), 0.173)
  # and one intercept
  expect_equal(round(fit_k2$intercepts[[1]][[1]], 3), 0.72)
  

})