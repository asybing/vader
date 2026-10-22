get_tpk <- function(doses = NULL, doseu = NULL, dar = NULL, ii = NULL, addl = NULL, dltwindow = NULL, mod = NULL, pkm = "CMAXMMAE") { 
  
  tbs <- switch(doseu, mg = 1, mgkg = 80, mgkgaibw = 60) # typical body size
  
  mod %>% zero_re() %>% 
          param(DAR0 = dar) %>%
          ev(data.frame(ID = seq_along(doses), 
                       amt = doses * tbs, 
                      time = 0,
                       cmt = 1, 
                        ii = ii, 
                      addl = addl) %>% as.ev()) %>% 
          mrgsim_df(end = dltwindow * 24, delta = 1) %>% 
          group_by(ID) %>%
          summarise(CMAXMMAE = max(MMAE),
                    AUCMMAE  = last(AUCMMAE) / 1000,
                    CAVGMMAE = last(AUCMMAE) / last(time)) %>%
          ungroup() %>%
          pull(!!sym(pkm))
  
}
