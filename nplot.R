nplot <- function(popsums = NULL, n = NULL, nreps = NULL, xlab = "", patch.title = "") { 
  
  medn <- popsums %>% pull(!!ensym(n)) %>% median()
  maxn <- popsums$N %>% max()
  
  popsums %>% count(!!ensym(n)) %>% 
              mutate(PERC = 100 * n / nreps) %>% 
              ggplot(aes(!!ensym(n), PERC)) + 
              geom_col(fill = NA, color = darkblue) + 
              geom_vline(aes(xintercept = medn, color = "Median"), linetype = "dotted") + 
              ggtitle(patch.title) + 
              scale_x_continuous(name = xlab) + 
              scale_y_continuous(name = paste0("% of reps (N=", nreps, ")"), limits = c(0, 100), breaks = seq(0, 100, 25), labels = paste0(seq(0, 100, 25), "%"))
  
}