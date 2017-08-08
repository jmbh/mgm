
# Gives summary of resampling object

print.resample <- function(x, # output from resample()
                             ...)
  
  
{
  
  obj_class <- class(x$call$object)
  
  model_classes <- c('Mixed Graphical Model (MGM)',
                     'mixed Vector Autoregressive (mVAR) model',
                     'Time-varying Mixed Graphical Model (tv-MGM)',
                     'Time-varying mixed Vector Autoregressive (tv-mVAR) model')
  
  if('core' %in% obj_class) obj_class_display <- model_classes[1]
  if('mvar' %in% obj_class) obj_class_display <- model_classes[2]
  if('tvmgm' %in% obj_class) obj_class_display <- model_classes[3]
  if('tvmvar' %in% obj_class) obj_class_display <- model_classes[4]

  
  if('tvmvar' %in% obj_class | 'tvmgm' %in% obj_class) {
    
    cat('resample-object',
        '\n\nModel class: ', obj_class_display,
        '\nBootstrap samples: ' , x$call$nB,
        '\nBlocks: ' , x$call$blocks)
    
  } else {
    
    cat('resample-object',
        '\n\nModel class: ', obj_class_display,
        '\nBootstrap samples: ' , x$call$nB)
    
  }
  
    
    
} # eoF


