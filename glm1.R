source("msvm.R")
library(glmnet)

files <- c( "dna", 
            "letter", 
            "mnist", 
            "satimage", 
            "segment", 
            "usps", 
            "waveform")

my_pmap <- function(l, f = rbind.data.frame) {
    stopifnot(is.list(l))
    stopifnot(is.list(l[[1]]))
    lapply (1:length(l[[1]]), function(i) {
        ml <- lapply(1:length(l), function(j) l[[j]][[i]])
        do.call(f, ml)
    })
}

glm1 <- function(files = files, workers = 12) {
    plan(multicore, workers = workers)

    map (files, function(file) {
        map (c(300, 800), function(n) {
            future_map (0:19, ~ {
                print(.)
                dfs <- read_wlws(n, dataset = file, .)
                model_glm1(dfs)
            }) |> my_pmap()
        }) |> my_pmap()
    }) |> my_pmap()
}

e <- new.env(parent = emptyenv()) 
e[["wlw2"]] = wu2_ld
e[["normal"]] = normal_ld
e[["radial"]] = stratified_ld


model_glm1 <- function(dfs, alpha = 0) {
    df <- dfs$train
    dft <- dfs$test
    n <- dfs$n
    dataset <- dfs$dataset
    run <- dfs$run

    stopifnot(min(df$class_id) == 1)
    K <- max(df$class_id)

    N <- nrow(dft)

    r <- array(0, dim = c(K, K, N))
    print('done')

    min_lambda <- Inf
    max_lambda <- -Inf

    binary = data.frame(i = vector("integer", 0),
                        j = vector("integer", 0),
                        lambda = vector("numeric", 0)
                        )
                        
    for (i in 1:(K-1)) {
        for (j in (i+1):K) {
            df_ij <- filter(df, class_id %in% c(i,j))
            x <- as.matrix(select(df_ij, -class_id))

            N1 <- sum(df_ij$class_id == i)
            N2 <- sum(df_ij$class_id == j)

            r1 = 1 / (2 * N2 + 2)
            r2 = (2 * N1 + 1) / (2 * N1 + 2)

            v1 <- if_else(df_ij$class_id == i, 1 - r2, r2)
            v2 <- if_else(df_ij$class_id == j, r1, 1 - r1)

            y <- as.matrix(cbind(v2, v1))
            # y <- as.numeric(df$class_id == i) * (1 - 1/ nrow(df)) + 0.5 /nrow(df)
            cv1  <- cv.glmnet(x, y, family = binomial(), alpha = alpha, standardize = F)
            lambda <- cv1$lambda.min
            if (lambda == min(cv1$lambda)) {
                print(sprintf("lower too large %d %d", i, j))
            }
            if (lambda == max(cv1$lambda)) {
                print(sprintf("upper too small %d %d", i, j))
            }
            print(lambda)
            model <- glmnet(x, y, family = binomial(), alpha = 0, lambda = lambda, nlambda = 1, standardize = F)
            print("model")
            yt <- as.matrix(select(dft, -class_id))
            r[i,j, ] = predict(model, yt, type = "link")

            binary[nrow(binary) + 1, ] <- data.frame(i = i, j = j, lambda = lambda)
        } # end loop j
    } # end loop i 
    multi = map (ls(e), function(m) {
        p <- sapply(1:N, function(k) {
                fn <- e[[m]]
                vecp <- fn(r[,,k])
                which.max(vecp)
                })
        data.frame(n = n, method = m, dataset = dataset, run = run, lambda = lambda,  correct = sum(p == dft$class_id))
    } )|> list_rbind()
    list(multi = multi, binary = binary)
}
