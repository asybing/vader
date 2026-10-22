get_udss <- function(pop, cohsum, topcoh, tdlt, maxnc, maxnt, detype = "3+3") {
  
  id     <- pop %>% slice_max(ID, n = 1) %>% pull(ID)
  coh    <- pop %>% slice_max(ID, n = 1) %>% pull(COH)
  csum   <- cohsum %>% filter(COH == coh) 
  closum <- cohsum %>% filter(COH == coh - 1) 
  cupsum <- cohsum %>% filter(COH == coh + 1) 
  n      <- csum$N
  ndlt   <- csum$NDLT
  total  <- pop %>% nrow()
  rdlt   <- csum$RDLT
  
  ### 3+3 ----------------------------------------------------------------------
  if (detype == "3+3") { 
    
    # Always enroll at least 3 subjects -------------------------------
            if (n < 3) { data.frame(ID = id, D = "stay", REASON = "Need at least 3 subjects per cohort") 
    
    # Special case rules: first dose and highest dose -----------------
     } else if (coh == 1 & ndlt >= 2) { data.frame(ID = id, D = "stop", REASON = "first dose too toxic, stop per 3+3") 
     } else if (coh == topcoh & n <= 6) { 
            if (ndlt <= 1) {            data.frame(ID = id, D = "stay", REASON = "backfill at top cohort")
     } else if (ndlt >= 2) {            data.frame(ID = id, D = "down", REASON = "stop per 3+3, backfill if needed") 
     
    # backfill at highest dose if needed -----------------------------
    }} else if (coh == topcoh) { 
  
            if (ndlt >= 3)  {           data.frame(ID = id, D = "down", REASON = "down per 3+3, backfill if needed") 
     } else if (n >= maxnc) {           data.frame(ID = id, D = "stop", REASON = paste("N =", maxnc, "per cohort reached")) 
     } else if (ndlt <= 2)  {           data.frame(ID = id, D = "stay", REASON = "backfill up to 12 subjects")
       
    # if there are 3 subjects ----------------------------------------         
    }} else if (n == 3) {
            if (ndlt == 0) {            data.frame(ID = id, D = "up",   REASON = "go up per 3+3") 
     } else if (ndlt == 1) {            data.frame(ID = id, D = "stay", REASON = "stay per 3+3") 
     } else if (ndlt >= 2) {            data.frame(ID = id, D = "down", REASON = "down per 3+3, backfill if needed") 
    
    # if there are 3 < n <= 6 subjects -------------------------------
    }} else if (n <= 6) { 
            if (ndlt <= 1 & n == 6) {   data.frame(ID = id, D = "up",   REASON = "go up per 3+3") 
     } else if (ndlt <= 1 & n <  6) {   data.frame(ID = id, D = "stay", REASON = "stay per 3+3") 
     } else if (ndlt >= 2)          {   data.frame(ID = id, D = "down", REASON = "down per 3+3, backfill if needed") 
       
    }} 
    
  ### BOIN ---------------------------------------------------------------------
  } else if (detype == "boin") { 
    
    phi1 <- 0.6*tdlt; phi2 <- 1.4*tdlt
    le   <- log((1 - phi1)/(1 - tdlt))/log(tdlt*(1 - phi1)/phi1/(1 - tdlt)) 
    ld   <- log((1 - tdlt)/(1 - phi2))/log(phi2*(1 - tdlt)/tdlt/(1 - phi2))
    
    # Always enroll at least 3 subjects -----------------------------
           if (n < 3) {                data.frame(ID = id, D = "stay", REASON = "Need at least 3 subjects per cohort") 
             
    # Special case rules: N=3, first dose, and highest dose ---------
    } else if (n == 3 & ndlt == 1)   { data.frame(ID = id, D = "stay", REASON = "For N=3, 3+3 rule overrides BOIN rule") 
    } else if (coh == 1 & rdlt > ld) { data.frame(ID = id, D = "stop", REASON = "Early stop due to de-escalation decision at 1st dose per BOIN.")
    } else if (coh == topcoh) { 
      
             if (n    >= maxnc) {      data.frame(ID = id, D = "stop", REASON = paste(maxnc, "subjects at top dose")) 
      } else if (rdlt <  ld)    {      data.frame(ID = id, D = "stay", REASON = paste("At top dose, keep enrolling until", maxnc, "subjects"))  
      } else if (rdlt >  ld)    {      data.frame(ID = id, D = "down", REASON = "De-escalate per BOIN")  
      
    # Early stopping rules for cohort N and total N -----------------
   }} else if (n     >= maxnc) {       data.frame(ID = id, D = "stop", REASON = paste("N =", maxnc, "per cohort reached")) 
    } else if (total >= maxnt) {       data.frame(ID = id, D = "stop", REASON = paste("N =", maxnt, "for the trial reached")) 
      
    # BOIN rules ----------------------------------------------------
    } else if (rdlt > ld) {            data.frame(ID = id, D = "down", REASON = "Go down per BOIN") 
    } else if (rdlt < le) {            data.frame(ID = id, D = "up",   REASON = "Go up per BOIN") 
    } else                             data.frame(ID = id, D = "stay", REASON = "Stay per BOIN") 
    
  ### BLRM ---------------------------------------------------------------------  
  } else if (detype == "blrm") {
    
    # Always enroll at least 3 subjects -----------------------------
           if (n < 3) { data.frame(ID = id, D = "stay", REASON = "Need at least 3 subjects per cohort") 
      
    # Special case rules: first cohort use 3+3 rules ----------------
    } else if (coh == 1 & n <= 6) {
      
              if (n == 3) {
              if (ndlt == 0) {          data.frame(ID = id, D = "up",   REASON = "Use 3+3 rules for first dose level. go up per 3+3") 
       } else if (ndlt == 1) {          data.frame(ID = id, D = "stay", REASON = "Use 3+3 rules for first dose level. stay per 3+3") 
       } else if (ndlt >= 2) {          data.frame(ID = id, D = "stop", REASON = "Use 3+3 rules for first dose level. stop per 3+3") 
        
      }} else {
              if (ndlt <= 1 & n == 6) { data.frame(ID = id, D = "up",   REASON = "Use 3+3 rules for first dose level. go up per 3+3") 
       } else if (ndlt <= 1 & n <  6) { data.frame(ID = id, D = "stay", REASON = "Use 3+3 rules for first dose level. stay per 3+3") 
       } else if (ndlt >= 2)          { data.frame(ID = id, D = "stop", REASON = "Use 3+3 rules for first dose level. stop per 3+3") 
      }}
      
    # Special case rules: top cohort ---------------------------------
    } else if (coh == topcoh & !csum$UNTOL) { if (       n >= maxnc) data.frame(ID = id, D = "stop", REASON = "Max number of subjects reached at top cohort.") else data.frame(ID = id, D = "stay", REASON = "Keep enrolling until maxnc is reached.") 
    
    # BLRM rules -----------------------------------------------------
    } else if (                 csum$UNTOL) { if (closum$N >= maxnc) data.frame(ID = id, D = "stop", REASON = "Dose is untolerable and N for DL-1 is maxed.")  else data.frame(ID = id, D = "down", REASON = "Dose level is untolerable. go down.") 
    } else if (               cupsum$UNTOL) { if (       n >= maxnc) data.frame(ID = id, D = "stop", REASON = "Max N reached, next dose level is untolerable") else data.frame(ID = id, D = "stay", REASON = "Next dose level is untolerable")   
    } else                                                           data.frame(ID = id, D = "up",   REASON = "dose is tolerable and next dose level is tolerable") 
  
  } else stop("detype must be '3+3', 'boin', or 'blrm'")
  
}

