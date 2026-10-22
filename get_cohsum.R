get_cohsum <- function(pop, ncoh,
                       doses = NULL, doseu = NULL, ii = NULL, addl = NULL, dltwindow = NULL, dar = NULL, mod = NULL, pkm = "CMAXMMAE", 
                       tdlt = 0.25, detype = "3+3", noprior = FALSE) { 

  pop %>% mutate(COH = factor(COH, levels = 1:ncoh)) %>% 
          group_by(COH, DOSE, .drop = F) %>% 
          summarize(N      = n(), 
                    NDLT   = sum(DLT), 
                    RDLT   = NDLT/N %>% round(3), 
                   .groups = "drop_last") %>% 
          ungroup() %>% 
          mutate(CONF = switch(detype, "3+3" = NA,
                                        boin = 1 - pbeta(tdlt, shape1 = NDLT + 1, shape2 = N - NDLT + 1) %>% round(3),
                                        blrm = pop %>% get_blrmconf(doses = doses, doseu = doseu, ii = ii, addl = addl, dltwindow = dltwindow, dar = dar, mod = mod, pkm = pkm, tdlt = tdlt, noprior = noprior)),
                UNTOL = switch(detype, "3+3" = (N >= 3 & N <= 6 & NDLT >= 2) | (N > 6 & NDLT >= 3),
                                        boin =  N >= 3      & CONF > 0.95,
                                        blrm = !is.na(CONF) & CONF > 0.25))
  
}
