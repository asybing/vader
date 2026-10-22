get_sub <- function(id, coh, dose, doseu = NULL, ii = NULL, addl = NULL, dltwindow = NULL, dar = NULL, mod = NULL, pkm = "CMAXMMAE") { 
  
  # Simulate PK ----------------------------------------------------------------
  bs  <- get_bs()
  bsu <- switch(doseu, mg = 1, mgkg = bs["bw"], mgkgaibw = bs["aibw"], mgm2 = bs["bsa"]) 
  pk  <- mod %>% param(DAR0 = dar) %>%
                 ev(amt = dose * bsu, cmt = 1, ii = ii, addl = addl) %>% 
                 mrgsim_df(end = dltwindow * 24, delta = 1) %>% 
                 summarise(CMAXMMAE = max(MMAE),
                           AUCMMAE  = last(AUCMMAE) / 1000,
                           CAVGMMAE = last(AUCMMAE) / last(time)) %>%
                 ungroup() %>%
                 pull(!!sym(pkm)) 
  
  # Simulate DLT ---------------------------------------------------------------      
  betas <- get_betas(pkm)
  pdlt  <- exp(betas[1] + betas[2]*pk)/(1 + exp(betas[1] + betas[2]*pk))
  dlt   <- rbinom(n = 1, size = 1, prob = pdlt)
  
  data.frame(ID     = id,
             COH    = coh, 
             DOSE   = dose, 
             ISMALE = bs[["ismale"]],
             BW     = bs[["bw"]], 
             HT     = bs[["ht"]], 
             AIBW   = bs[["aibw"]], 
             PK     = pk, 
             PDLT   = pdlt, 
             DLT    = dlt)
    
}
