### ----------------------------------------------------------------------------
### simde : Single Dose Escalation (DE) Clinical Trial Simulator (sim) for Vedotin ADCs 
### ----------------------------------------------------------------------------

# Workflow ---------------------------------------------------------------------
# simde()            :  
# |-check_inputs()   : 1.) Conduct input and data checks before proceeding
# |-get_tmtd()       : 2.) Calculate the 'true' MTD based on the assumed 
# | |-get_betas()    :     exposure-response relationship and the typical 
# | |-get_tpk()      :     exposure based on the poppk model provided
# |-get_sub()        : 3.) enroll a simulated subject with a BW, HT, BSA, AIBW
# | |-get_bs()       :     exposure, and if they experienced a DLT.
# | |-get_betas()    :
# |-get_cohsum()     : 4.) Summarise the cohort level data including which N, NDLT,  
# | |-get_blrmconf() :     RDLT, and which dose levels are toxic (based on 3+3, 
# | | |-get_betas()  :     BOIN, or BLRM)
# | | |-get_tpk()    : 
# |-get_udss()       : 5.) Decide whether to go up, down, stay, or stop (UDSS)

# Inputs -----------------------------------------------------------------------
#
# Inputs for a single dose escalation trial
# Drug Parameters 
# * doses: A vector of numerics for the dose levels to investigate. Default of NULL.
# * doseu: A string of "mg", "mgkg", "mgkgaibw", or "mgm2" for dose units. Default of NULL.
# * ii: A numeric in hours of the frequency of dosing during the DLT window. Default of NULL.
# * addl: An integer # of additional doses to be administered in the DLT window. Default of NULL.
# * dltwindow: A numeric in days of the length of the DLT window. Default of NULL.
# * dar: A numeric of the drug-antibody ratio (DAR) for the vedotin ADC. Default of NULL.
# * mod: An mrgsolve model object. Default of NULL.
# * pkm: A string of "CMAXMMAE", "AUCMMAE", or "CAVGMMAE" for the PK metric to base the exposure-response off of. Default of "CMAXMMAE".
#
# Clinical Trial Parameters
# * tdlt: A numeric target DLT rate. Default of 0.25.
# * maxnc: An integer maximum number of subjects per cohort. Default of 12.
# * maxnt: An integer maximum number of subjects per trial. Default of 90.
# * detype: A string of "3+3", "boin", or "blrm" for the dose escalation strategy. Default of "3+3".
#

simde <- function(doses = NULL, doseu = NULL, ii = NULL, addl = NULL, dltwindow = NULL, dar = NULL, mod = NULL, pkm = "CMAXMMAE", 
                  tdlt = 0.25, maxnc = 12, maxnt = 90, detype = "3+3", noprior = FALSE) {

  # Check inputs ---------------------------------------------------------------
  check_inputs(doses = doses, doseu = doseu, ii = ii, addl = addl, dltwindow = dltwindow, dar = dar, mod = mod, pkm = pkm, 
               tdlt = tdlt, maxnc = maxnc, maxnt = maxnt, detype = detype, noprior = noprior)
  
  # Initialize trial variables -------------------------------------------------
  tmtd   <- get_tmtd(doses = doses, doseu = doseu, ii = ii, addl = addl, dltwindow = dltwindow, dar = dar, mod = mod, pkm = pkm, tdlt = tdlt)   
  topcoh <- ncoh <- length(doses)    
  n.udss <- pop <- data.frame()
  id     <- coh <- 1; dose <- doses[1] # initialize first subject
  
  # Enroll patients until MTD or MAD is met ------------------------------------
  repeat { 
    
    sub    <- get_sub(id = id, coh = coh, dose = dose, doseu = doseu, ii = ii, addl = addl, dltwindow = dltwindow, dar = dar, mod = mod, pkm = pkm)
    pop    <- pop %>% rbind(sub)
    cohsum <- pop %>% get_cohsum(ncoh = ncoh, doses = doses, doseu = doseu, ii = ii, addl = addl, dltwindow = dltwindow, dar = dar, mod = mod, pkm = pkm, tdlt = tdlt, detype = detype, noprior = noprior)
    udss   <- pop %>% get_udss(cohsum = cohsum, tdlt = tdlt, topcoh = topcoh, maxnc = maxnc, maxnt = maxnt, detype = detype)
    n.udss <- n.udss %>% rbind(udss)
   
    if (udss$D == "stop") break else { 
      
      if (detype %in% c("3+3", "boin")) {
        
        untol <- cohsum %>% filter(DOSE == dose) %>% pull(UNTOL) # decide if this is an untolerable dose (ud)
        if (untol) { doses  <- doses[-(coh:topcoh)]
                     topcoh <- length(doses) }  
      
      }
      
      # Define next subject ----------------------------------------------------
      id   <- id + 1
      coh  <- switch(udss$D, up = coh + 1, down = coh - 1, stay = coh)
      dose <- doses[coh]
    
    }

  }
  
  # Declare MAD / MTD ----------------------------------------------------------
  if (cohsum$UNTOL[1]) { maxd  <- doses[1] 
                         maxdl <- "MAD" 
                } else { maxd  <- cohsum %>% filter(!UNTOL & N > 0) %>% pull(DOSE) %>% max()
                         maxdl <- if (any(cohsum$UNTOL)) "MTD" else "MAD" }

  # Output ---------------------------------------------------------------------
  list(pop    = inner_join(pop, n.udss, by = "ID"), 
       cohsum = cohsum, 
       n      = pop %>% nrow(), 
       ndlt   = pop %>% summarise(NDLT = sum(DLT))         %>% pull(NDLT),
       ntox   = pop %>% summarise(NTOX = sum(DOSE > tmtd)) %>% pull(NTOX), 
       maxd   = maxd, 
       maxdl  = maxdl, 
       tmtd   = tmtd)

}
