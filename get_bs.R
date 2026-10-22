# Simulates a single subject based on adult data from the NHANES database. 
# A male or female subject is generated from a 50% probability. From there,
# weight (BW) and height (HT) are simulated based on sex-specific distributions 
# and correlation of BW and HT. Common body size metrics of ideal body weight (IBW),
# adjusted IBW (AIBW), and body surface area (BSA) are calculated.
#
# Units for each body size metrics are reported in:
# * Body weight (BW)                  : kilograms
# * Height (HT)                       : inches 
# * Ideal body weight (IBW)           : kilograms
# * Adjusted ideal body weight (AIBW) : kilograms
# * body surface area (BSA)           : meters^2

get_bs <- function() { 
  
  ismale <- rbinom(n = 1, size = 1, prob = 0.5)
    
         if (ismale == 1) { corr <- matrix(c(1, 0.420, 0.420, 1), 2, 2) %>% chol()
                            etas <- matrix(rnorm(2), ncol = 2) %*% corr
                            bw   <- 86.4*exp(etas[1]*0.232)
                            ht   <- 68.6*exp(etas[2]*0.0445)
                            ibw  <- 50.0 + 2.3*(ht - 60)
                              
  } else if (ismale == 0) { corr <- matrix(c(1, 0.295, 0.295, 1), 2, 2) %>% chol()
                            etas <- matrix(rnorm(2), ncol = 2) %*% corr
                            bw   <- 75.3*exp(etas[1]*0.263)
                            ht   <- 63.2*exp(etas[2]*0.0440)
                            ibw  <- 45.5 + 2.3*(ht - 60) 
                              
  } else stop("'ismale' must be 0 or 1")
    
  c("bw"     = bw, 
    "bsa"    = (bw * ht * 2.54 / 3600)^0.5, 
    "aibw"   = ibw + 0.4 * (bw - ibw), 
    "ismale" = ismale, 
    "ht"     = ht, 
    "ibw"    = ibw)
    
}
