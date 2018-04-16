


print.int <- function(x,
                      ...)
 
  
  
{
  
  # Sign Defined?  
  if(x$sign == 0) sign <- "Not Defined" else sign <- x$sign 
  
  # Create Console Output
  cat('Interaction:', paste(x$variables, collapse = "-"),
        '\nWeight: ', x$edgeweight,
        '\nSign: ' , sign)

}
  
  