get_tmtd <- function(doses = NULL, doseu = NULL, ii = NULL, addl = NULL, dltwindow = NULL, dar = NULL, mod = NULL, pkm = "CMAXMMAE", 
                     tdlt = 0.25) {
  
  betas <- get_betas(pkm)
  tpk   <- get_tpk(doses = doses, doseu = doseu, ii = ii, addl = addl, dltwindow = dltwindow, dar = dar, mod = mod, pkm = pkm)
  
  max(doses[tpk < (log(tdlt/(1 - tdlt)) - betas[1])/betas[2]]) 

}
