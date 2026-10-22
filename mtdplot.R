mtdplot <- function(popsums = NULL, nreps = NULL, doses = NULL, tmtd = NULL) { 
  
  tmtdi <- which(doses == tmtd)
  
  popsums %>% mutate(MAXD  = factor(MAXD,  levels = doses),
                     MAXDL = factor(MAXDL, levels = c("MAD", "MTD"))) %>%  
              count(MAXDL, MAXD, .drop = F) %>% 
              mutate(PERC = 100 * n / nreps) %>% 
              ggplot(aes(MAXD, PERC)) + 
              geom_vline(aes(xintercept = tmtdi, color = "'True' MTD")) +
              geom_col(aes(fill = MAXDL), position = "stack") + 
              scale_x_discrete(name = "Observed MTD/MAD") + 
              scale_y_continuous(name = paste0("% of reps (N=", nreps, ")"), limits = c(0, 100), breaks = seq(0, 100, 25), labels = paste0(seq(0, 100, 25), "%")) 
  
}
