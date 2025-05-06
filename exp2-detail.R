library(tidyr)
library(dplyr)
library(ggplot2)

library(xtable)
source('lda.R')


dfs <- read_wlws(800, "letter", 0)

# W <- gen_W2(20)

#bcp <- get_bcp(dfs, 5, 12, 16)
bcp <- get_bcp(dfs, 5, 8, 12)

# scores <- eval_scores(bcp, W, score = "acc")$scores

M_proj <- matrix( c(cos(pi/2), sin(pi/2), 
                    cos(7*pi/6), sin(7* pi/6),
                    cos(11 * pi / 6), sin(11 * pi / 6)),
                    ncol = 2, byrow = T)

colnames(M_proj) <- c("x", "y")

# df_pos <- data.frame(W %*% M_proj)

# df_pos$score = apply(scores, 2, mean)

# plt <- ggplot(df_pos, aes(x = x, y = y, color = score)) + geom_point(shape = 21)

N <- 140
xy <- expand.grid(x = seq(-4,8, length.out = N), y = seq(-6,6, length.out = N))

CB <- cbind(M_proj[2,] - M_proj[1,], M_proj[3, ] - M_proj[1,])

print("here 1")
f <- function(w) {
        delta <- w - M_proj[1,]
            v = solve(CB, delta)
                (1-v[1] - v[2])*c(1,0,0) + v[1] * c(0,1,0) + v[2] * c(0,0,1)
}

W1 <- apply(xy, 1, f) |> t()
colnames(W1) <- c("a1", "a2", "a3")
W2 <- subset(W1, (W1[,1] >= 0) & (W1[,2] >=0 ) & (W1[,3] >=0))

print("here 2")
scores1 <- eval_scores(bcp, W1, score = "log")$scores
#scores1 <- eval_scores(bcp, W1, score = "brier")$scores
acc1 <- apply(eval_scores(bcp, W1, score = "acc")$scores,2, mean)
scores2 <- apply(eval_scores(bcp, W2, score = "log")$scores, 2, mean)
#scores2 <- apply(eval_scores(bcp, W2, score = "brier")$scores, 2, mean)
acc2 <- apply(eval_scores(bcp, W2, score = "acc")$scores, 2, mean)

df_pos1 <- data.frame(W1 %*% M_proj)

print("here 3")
df_pos1$score = apply(scores1, 2, mean)

imax = which.max(df_pos1$score)

imax2 = which.max(scores2)

res <- data.frame(optimum = c("Best in \u0394", "Best in D"), 
            `Brier score` = c(scores2[imax2], df_pos1$score[imax]),
            `accuracy` = c(acc2[imax2], acc1[imax]))

res <- cbind(res, rbind(W2[imax2,], W1[imax,]))

sink("tab-step5.tex")
print(xtable(res, label  = "tab:step5", digits = 2,
    caption = "Optima in the simplex and rectangle $D$"), include.rownames = F)
sink()


plt <- ggplot() +   geom_tile(data = df_pos1, mapping= aes(x = x, y = y,  fill = score)) + 
    coord_fixed() +
    annotate("point", x = df_pos1$x[imax], y = df_pos1$y[imax], shape = 21, fill = "black") + 
    geom_polygon(data = M_proj, aes(x = x, y= y), fill = NA, color = "black") + 
   scale_color_gradient(low = "green", high = "red") + xlab("") + ylab("") + 
    scale_fill_gradient(low = "green", high = "red") + theme_minimal(base_size = 9) + theme(
        axis.ticks = element_blank(),    # Remove tick marks
            axis.text = element_blank()      # Remove tick labels
              )

ggsave(
  filename = "graphs/exp2-detail.pdf",
  plot = plt,
  width = 5.5,   # Width in inches
  height = 3.5,  # Height in inches
  units = "in"   # Specify units as inches
)

