

print.mgm <- function(x,
                      ...)



{

  model_classes <- c('Mixed Graphical Model (MGM)',
                     'mixed Vector Autoregressive (mVAR) model',
                     'Time-varying Mixed Graphical Model (tv-MGM)',
                     'Time-varying mixed Vector Autoregressive (tv-mVAR) model')


  # ---------- print for fit objects ----------

  if(!('predicted' %in% class(x)) & !('bwSelect' %in% class(x))) {

    if(class(x) == 'mgm') {
      cat('mgm fit-object',
          '\n\nModel class: ', model_classes[1],
          '\nOrder: ' , x$call$k,
          '\nNodes: ' , length(x$call$type))
    }


    if(class(x) == 'mvar') {
      cat('mgm fit-object',
          '\n\nModel class: ', model_classes[2],
          '\nLags: ' , x$call$lags,
          '\nNodes: ' , length(x$call$type))
    }



    if(class(x) == 'tvmgm') {
      cat('mgm fit-object',
          '\n\nModel class: ', model_classes[3],
          '\nOrder: ' , x$call$k,
          '\nNodes: ' , length(x$call$type),
          '\nEstimation points: ' , length(x$call$estpoints))
    }


    if(class(x) == 'tvmvar') {
      cat('mgm fit-object',
          '\n\nModel class: ', model_classes[4],
          '\nLags: ' , x$call$lags,
          '\nNodes: ' , length(x$call$type),
          '\nEstimation points: ' , length(x$call$estpoints))
    }

  }


  # ---------- print for prediction object ----------

  if('predicted' %in% class(x)) {

    if('mgm' %in% class(x)) mc <- model_classes[1]
    if('mvar' %in% class(x)) mc <- model_classes[2]
    if('tvmgm' %in% class(x)) mc <- model_classes[3]
    if('tvmvar' %in% class(x)) mc <- model_classes[4]


    cat('mgm prediction-object',
        '\n\nModel class: ', mc,
        '\nError Types:', paste(names(x$call$errorCon), names(x$call$errorCat)))

  }


  # ---------- print for bwSelect object ----------


  if('bwSelect' %in% class(x)) {

    if('mgm' %in% class(x)) mc <- model_classes[3]
    if('mvar' %in% class(x)) mc <- model_classes[4]

    cat('mgm bandwidth selection-object',
        '\n\nModel class: ', mc,
        '\nBandwith path: ', paste(x$call$bwSeq),
        '\nNumber of Folds: ', paste(x$call$bwFolds),
        '\nFoldsize: ', paste(x$call$bwFoldsize),
        '\nOptimal Bandwidth: ', x$call$bwSeq[which.min(x$meanError)])

  }




} # eoF
