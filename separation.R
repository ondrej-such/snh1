library(CVXR)
source("msvm.R")
source("lda.R")

check_sep <- function(f = "dna") {
    print(sprintf("%s ", f ))
    map(0:19, function(run) {
        dfs <- read_wlws(800, f , run)
        dft <- dfs$train
        X <- as.matrix(dplyr::select(dft, -class_id))
        K <- max(dft$class_id)

        w <- Variable(ncol(dft) - 1, cols = 1)
        b <- Variable(1)

        df <- map (1:(K-1), function(i) {
            map ((i+1):K, function(j) {
                X1 <- X[dft$class_id == i,]
                X2 <- X[dft$class_id == j,]
                objective <- Minimize(sum_squares(w))
                constraints <- list (X1 %*% w + b <= -1,
                                     X2 %*% w + b >= 1)
                problem  <- Problem(objective, 
                        constraints = constraints)
                solvers <- c("ECOS")
                map (solvers, function(solver) {
                    result <- solve(problem, solver)
                    data.frame(i = i, j = j, status = result$status, solver = solver)
                }) |> list_rbind()
            }) |> list_rbind()
        }) |>  list_rbind()
        df$dataset = f
        df$K = K
        df
    }) |> list_rbind()
}

