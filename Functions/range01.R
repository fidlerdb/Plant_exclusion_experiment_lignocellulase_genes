library(compiler)
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
cmpfun(range01)