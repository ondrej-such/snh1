library(tidyr)
library(dplyr)
library(xtable)

df <- read.csv("data/triples.csv")

df1 <- df |> mutate(omit = pmax(omit12, omit13, omit23)) |> 
            summarize(` $\\boldsymbol{e}_3$ at least as accurate as $\\boldsymbol{wlw}$` = mean(normal >= wlw2),
                      `one of $\\boldsymbol{s}_1, \\boldsymbol{s}_2, \\boldsymbol{s}_3$ at least as accurate as $\\boldsymbol{wlw}$` = mean(omit >= wlw2),
                      `$\\boldsymbol{wlw}$ more accurate than $\\boldsymbol{e}_3, \\boldsymbol{s}_1, \\boldsymbol{s}_2, \\boldsymbol{s}_3$` = mean(wlw2 > normal & wlw2 > omit)) |>
      pivot_longer(everything(), names_to = 'Outcome', values_to = "Portion")
            
                      
                      
                      
           
# colnames(df1) <- c("dataset", "normal", "radial", "Wu-Lin-Weng")
sink("tab-step2.tex")
print(xtable(df1, label = "tab:step2", digits = 2,
             caption = "Comparison of Wu-Lin-Weng's method $\\boldsymbol{wlw}$ and  parameterless coupling methods on three class subsets. "), 
    sanitize.text.function = identity,
    include.rownames = FALSE)
sink()
