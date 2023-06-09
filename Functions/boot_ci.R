
# Create a function to produce bootstrapped 95% confidence intervals
# ()
boot_ci <- function(x, B = NULL){
  # Deal with a NULL bootstrap value
  if(is.null(B)){B <- 1000}
  
  # Find how many rows we have 
  n = length(x)
  
  # Create a matrix of samples to populate with recored values
  boot.samples = matrix(sample(x, size = n*B, replace = TRUE)
                        , B, n)
  # Calculate the mans of those values
  boot.statistics = apply(boot.samples, 1, mean)
  
  # Take the 2.5%-ile and the 97.5%-ile and call these u_ci and l_ci
  l_ci <- quantile(boot.statistics, 0.025)
  u_ci <- quantile(boot.statistics, 0.975)
  
  # Give these back to the user
  return(data.frame(l_ci = l_ci, u_ci = u_ci))
}

# Use byte compilation to speed this up
require(compiler)
boot_ci <- cmpfun(boot_ci)

# boot.mean(x = df[df$Species == 'Sorangium cellulosum' 
#                   & df$Treatment == 'b-New',]$Abundance#, B = 10000
#           )
