# Retreives the exposure-response parameters based on the
# MBMA logistic regression analysis across vedotin ADCs for 
# different MMAE exposure metrics.
# * B0 (intercept) is the first vector value
# * B1 (slope) is the second vector value 

get_betas <- function(pkm = "CMAXMMAE") { 
  
  switch(pkm,
         CMAXMMAE = c(-4.895, 0.433),
         AUCMMAE  = c(-4.008, 1.144),
         CAVGMMAE = c(-4.252, 0.719))

}