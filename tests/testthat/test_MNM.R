

test_that("MNM with k-specification", {
  
  data(autism_data)
  
  # Fit MGM with pairwise & three-way interactions
  set.seed(1)
  fit_k3 <- mgm(data = autism_data$data,
                type = autism_data$type,
                level = autism_data$lev,
                lambdaSel = "CV",
                lambdaFolds = 10,
                k = 3, 
                pbar = FALSE, 
                signInfo = FALSE)
  
  # Number of 3-way interactions = 11 
  expect_equal(nrow(fit_k3$interactions$indicator[[2]]), 10)
  
  
  # Plot Factor Graph (only to catch errors ...)
  FactorGraph(object = fit_k3, 
              PairwiseAsEdge = FALSE, 
              labels = autism_data$colnames)
  

})