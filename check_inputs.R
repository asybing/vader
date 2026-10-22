check_inputs <- function(reps = 500, 
                         doses = NULL, doseu = NULL, ii = NULL, addl = NULL, dltwindow = NULL, dar = NULL, mod = NULL, pkm = "CMAXMMAE", 
                         tdlt = 0.25, maxnc = 12, maxnt = 90, detype = "3+3", noprior = FALSE,
                         patch.title = "") { 
                         
  # Check for required packages ----------------------------------------------- 
  if(!requireNamespace("tidyverse", quietly = T)) stop("Please load tidyverse before proceeding")
  if(!requireNamespace("mrgsolve",  quietly = T)) stop("Please load the mrgsolve package before proceeding")
  if(!requireNamespace("rstanarm",  quietly = T)) stop("Please load the rstanarm package before proceeding")
  
  # Inputs for Monte Carlo Simulations Check ----------------------------------
  if (!is.numeric(reps) | reps %% 1 != 0 | length(reps) != 1 | reps <= 0) {
    stop("Error: 'reps' must be a single positive integer.")
  }
  
  if (!is.character(patch.title) | length(patch.title) != 1) {
    stop("Error: 'patch.title' must be a single string.")
  }
  
  # Drug Parameter Check ------------------------------------------------------
  if (!is.numeric(doses) | any(doses <= 0) | is.unsorted(doses, strictly = T) | any(is.na(doses))) {
    stop("Error: 'doses' must be a numeric vector with positive values in ascending order.")
  }

  valid_doseu <- c("mg", "mgkg", "mgkgaibw", "mgm2")
  if (!is.character(doseu) | length(doseu) != 1 | !(doseu %in% valid_doseu)) {
    stop("Error: 'doseu' (dose units) must be a single string. Allowed values: ", paste(valid_doseu, collapse = ", "))
  }

  if (!is.numeric(ii) | length(ii) != 1 | ii <= 0) {
    stop("Error: 'ii' (dosing interval) must be a single positive numeric value.")
  }

  if (!is.numeric(addl) | addl %% 1 != 0 | length(addl) != 1 | addl < 0) {
    stop("Error: 'addl' (additional doses) must be a single non-negative integer.")
  }

  if (!is.numeric(dltwindow) | length(dltwindow) != 1 | dltwindow <= 0) {
    stop("Error: 'dltwindow' must be a single positive numeric value.")
  }
 
  if (!is.numeric(dar) | length(dar) != 1 | dar <= 0) {
    stop("Error: 'dar' (drug-antibody ratio) must be a single positive numeric value.")
  }

  if (!inherits(mod, "mrgmod")) {
    stop("Error: 'mod' (PopPK model) must be an 'mrgsolve' model object (class 'mrgmod').")
  }
  
  if(!("MMAE" %in% mod$capture)) {
    stop(paste("ensure 'MMAE' is in the $CAPTURE block of the model"))
  }
  
  valid_pkm <- c("CMAXMMAE", "AUCMMAE", "CAVGMMAE")
  if (!is.character(pkm) | length(pkm) != 1 | !(pkm %in% valid_pkm)) {
    stop("Error: 'pkm' (PK metric) must be a single string. Allowed values: ", paste(valid_pkm, collapse = ", "))
  }

  # Clinical Trial Parameter Check --------------------------------------------
  if (!is.numeric(tdlt) | length(tdlt) != 1 | tdlt < 0 | tdlt > 1) {
    stop("Error: 'tdlt' (target DLT rate) must be a single numeric value between 0 and 1.")
  }

  if (!is.numeric(maxnc) | maxnc %% 1 != 0 | length(maxnc) != 1 | maxnc <= 0) {
    stop("Error: 'maxnc' (cohort enrollment limit) must be a single positive integer.")
  }

  if (!is.numeric(maxnt) | maxnt %% 1 != 0 | length(maxnt) != 1 | maxnt <= 0) {
    stop("Error: 'maxnt' (trial enrollment limit) must be a single positive integer.")
  }

  if (!is.null(maxnc) & !is.null(maxnt) & maxnc > maxnt) {
    stop("Error: 'maxnc' (cohort enrollment limit) cannot be greater than 'maxnt' (trial enrollment limit).")
  }

  valid_detype <- c("3+3", "boin", "blrm")
  if (!is.character(detype) | length(detype) != 1 | !(detype %in% valid_detype)) {
    stop("Error: 'detype' (dose escalation type) must be a single string. Allowed values: ", paste(valid_detype, collapse = ", "))
  }
  
  if (!is.logical(noprior)) { 
    stop("Error: 'noprior' must be a TRUE or FALSE statement.")
  }
  
}
