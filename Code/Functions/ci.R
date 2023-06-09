ci <- function(x){
  x <- na.omit(x)
  CI <- 1.96*(sd(x)/sqrt(length(x)))
  return(CI)
}

library("compiler")
ci <- cmpfun(ci)