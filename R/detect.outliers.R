# Author: Komal S Rathi
# Date: 08/15/2019
# Function: using Tukey's method for outlier detection
detect.outliers <- function(x){
  iqr <- 1.5*IQR(x)
  
  # overexpression cutoff
  perc.75 <- quantile(x, probs = 0.75)[[1]]
  over <- perc.75 + iqr
  
  # underexpression cutoff
  perc.25 <- quantile(x, probs = 0.25)[[1]]
  under <- perc.25-iqr
  
  x <- ifelse(x < under, -1, 
              ifelse(x > over, 1, 0))
  return(x)
}