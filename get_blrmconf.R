get_blrmconf <- function(pop, 
                         doses = NULL, doseu = NULL, ii = NULL, addl = NULL, dltwindow = NULL, dar = NULL, mod = NULL, pkm = "CMAXMMAE", 
                         tdlt = 0.25, noprior = FALSE) {
  
  if (nrow(pop) > 1) { 
  
    betas <- get_betas(pkm)
    tpk   <- get_tpk(doses = doses, doseu = doseu, ii = ii, addl = addl, dltwindow = dltwindow, dar = dar, mod = mod, pkm = pkm)
    
    if (noprior) { prior_intercept <- default_prior_intercept(binomial())
                   prior           <- default_prior_coef(binomial())
          } else { prior_intercept <- normal(location = betas[1] + betas[2]*mean(pop$PK), scale = 0.1)
                   prior           <- normal(location = betas[2],                         scale = 0.1) } 
    
    stan_glm(formula         = DLT ~ PK, 
             data            = pop, 
             family          = binomial(link = "logit"), 
             prior_intercept = prior_intercept,
             prior           = prior, 
             refresh         = 0) %>% 
      posterior_epred(newdata = data.frame(PK = tpk)) %>% 
      {. > tdlt + 0.08} %>% 
      colMeans() 
  
  } else NA
  
}
