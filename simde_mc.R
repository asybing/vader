### ----------------------------------------------------------------------------
### simde_mc : Monte Carlo (MC) Dose Escalation (DE) Clinical Trial Simulator (sim) for Vedotin ADCs 
### ----------------------------------------------------------------------------

# Workflow ---------------------------------------------------------------------
# simde_mc()           : 1.) Simulate N monte carlo (MC) clinical trial simulations
# |-check_inputs()     : 2.) Conduct input and data checks before proceeding
# |-simde()            :  
# | |-check_inputs()   : 
# | |-get_tmtd()       : 3.) Calculate the 'true' MTD based on the assumed 
# | | |-get_betas()    :     exposure-response relationship and the typical 
# | | |-get_tpk()      :     exposure based on the poppk model provided
# | |-get_sub()        : 4.) enroll a simulated subject with a BW, HT, BSA, AIBW
# | | |-get_bs()       :     exposure, and if they experienced a DLT.
# | | |-get_betas()    :
# | |-get_cohsum()     : 5.) Summarise the cohort level data including which N, NDLT,  
# | | |-get_blrmconf() :     RDLT, and which dose levels are toxic (based on 3+3, 
# | | | |-get_betas()  :     BOIN, or BLRM)
# | | | |-get_tpk()    : 
# | |-get_udss()       : 6.) Decide whether to go up, down, stay, or stop (UDSS)
# |-nplot()            : 7.) Generate summaries of the MC simulations (N enrolled, 
# |-mtdplot()          :     N at toxic doses, N DLTs, and how often the 'true' MTD 
#                      :     was identified).
# 
# Inputs -----------------------------------------------------------------------
#
# ## Inputs for the monte carlo simulations
# * reps: An integer of the # of monte carlo simulations to run. Default of 500.
# * patch.title: A string of the title for the final overall plot. Default of "".
#
# ## Inputs for a single dose escalation trial
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
# * noprior: a TRUE or FALSE used only for BLRM, if a prior should be used or not.

simde_mc <- function(reps = 500,
                     doses = NULL, doseu = NULL, ii = NULL, addl = NULL, dltwindow = NULL, dar = NULL, mod = NULL, pkm = "CMAXMMAE",
                     tdlt = 0.25, maxnc = 12, maxnt = 90, detype = "3+3", noprior = FALSE,
                     patch.title = "") {
  
  # Check inputs ---------------------------------------------------------------
  check_inputs(reps = reps, 
               doses = doses, doseu = doseu, ii = ii, addl = addl, dltwindow = dltwindow, dar = dar, mod = mod, pkm = pkm, 
               tdlt = tdlt, maxnc = maxnc, maxnt = maxnt, detype = detype, noprior = noprior,
               patch.title = patch.title)
  
  # Run dose escalations N=reps times ------------------------------------------
  sims <- map(1:reps, ~ simde(doses = doses, doseu = doseu, ii = ii, addl = addl, dltwindow = dltwindow, dar = dar, mod = mod, pkm = pkm,
                              tdlt = tdlt, maxnc = maxnc, maxnt = maxnt, detype = detype, noprior = noprior), 
                              .progress = T)
  
  # Extract data into data.frames ----------------------------------------------
  pops    <- map_df(sims, ~ .x$pop,    .id = "REP")
  cohsums <- map_df(sims, ~ .x$cohsum, .id = "REP")
  popsums <- map_df(sims, ~ data.frame(N     = .x$n,
                                       NDLT  = .x$ndlt,
                                       NTOX  = .x$ntox,
                                       MAXDL = .x$maxdl,
                                       MAXD  = .x$maxd, 
                                       TMTD  = .x$tmtd), .id = "REP")
  
  # Summarise the clinical trial simulations -----------------------------------
  sum <- popsums %>% summarise(across(c(N, NDLT, NTOX, MAXD), 
                               ~ paste0(round(  median(.x              ), 1), " (", 
                                        round(quantile(.x, probs = 0.05), 1), "-", 
                                        round(quantile(.x, probs = 0.95), 1), ")"), .names = "{.col}.MED.P05.P95"),
                     PERC.MAXD = paste0(round(    mean(MAXD == TMTD)*100, 1), "%")) %>% 
                     select(sort(names(.)))
  
  # Plot outputs ---------------------------------------------------------------
     n.plot <-   nplot(popsums, N,    nreps = reps, xlab = "N enrolled", patch.title = patch.title)
  ndlt.plot <-   nplot(popsums, NDLT, nreps = reps, xlab = "N with DLT")
  ntox.plot <-   nplot(popsums, NTOX, nreps = reps, xlab = "N above 'true' MTD")
   mtd.plot <- mtdplot(popsums,       nreps = reps, doses = doses, tmtd = unique(popsums$TMTD)) 
      patch <- n.plot + ndlt.plot + ntox.plot + mtd.plot + plot_layout(nrow = 1, axes = "collect", guides = "collect") 

  # Output ----------------------------------------------------------------------
  list(pops = pops, cohsums = cohsums, popsums = popsums, sum = sum, patch = patch, plots = list(n.plot, ndlt.plot, ntox.plot, mtd.plot),  
       reps = reps, doses = doses, doseu = doseu, ii = ii, addl = addl, dltwindow = dltwindow, dar = dar, mod = mod, pkm = pkm,
       tdlt = tdlt, maxnc = maxnc, maxnt = maxnt, detype = detype)
  
}

