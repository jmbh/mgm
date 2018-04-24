


print.int <- function(x,
                      ...)
 
  
  
{
  
  # Sign Defined?  
  if(x$sign == 0) sign <- "NA" else sign <- x$sign 
  
  # Signtext
  if(x$sign == 0) signtext <- " (Not defined)"
  if(x$sign == 1) signtext <- " (Positive)"
  if(x$sign == -1) signtext <- " (Negative)"
  
  # Create Console Output
  cat('Interaction:', paste(x$variables, collapse = "-"),
        '\nWeight: ', x$edgeweight,
        '\nSign: ' , sign, signtext)

}
  
  