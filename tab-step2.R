library(tidyr)
library(dplyr)
library(xtable)

df <- read.csv("data/triples.csv")

df1 <- df |> mutate(omit = pmax(omit12, omit13, omit23)) |> 
            summarize(`normal $\\geq$ Wu-Lin-Weng` = mean(normal >= wlw2),
                      `$\\max(\\boldsymbol{s}_1, \\boldsymbol{s}_2, \\boldsymbol{s}_3) \\geq$ Wu-Lin-Weng` = mean(omit >= wlw2),
                      `Wu-Lin-Weng > max(normal, $\\boldsymbol{s}_1, \\boldsymbol{s}_2, \\boldsymbol{s}_3$)` = mean(wlw2 > normal & wlw2 > omit)) |>
      pivot_longer(everything(), names_to = 'outcome', values_to = "percentage")
            
                      
                      
                      
           
# colnames(df1) <- c("dataset", "normal", "radial", "Wu-Lin-Weng")
sink("tab-step2.tex")
print(xtable(df1, label = "tab:step2", digits = 2,
             caption = "Comparison of Wu-Lin-Weng's method and  parameterless coupling methods on three class subset"), 
    sanitize.text.function = identity,
    include.rownames = FALSE)
sink()
