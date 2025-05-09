source("msvm.R")
library(glmnet)
library(CalibratR)

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
    r <- lapply (1:length(l[[1]]), function(i) {
        ml <- lapply(1:length(l), function(j) l[[j]][[i]])
        do.call(f, ml)
    })
    names(r) <- names(l[[1]])
    r
}

save_glm1 <- function(... ) {
    l <- glm1(...)
    dir <- "data"
    prefix <- sprintf("%s/glm1-", dir)
    where <- sprintf("%s%s.csv", prefix, c("binary", "multi"))
    print(where)
    write_csv(l$binary, file = where[[1]])
    write_csv(l$multi, file = where[[2]])
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


omit1 <- function(M_logits, idx) {
    inverse_logit <- function(x) {
          return(1 / (1 + exp(-x)))
    }
    M <- matrix(0, nrow = 3, ncol = 3)

    for (i in 1:3) {
        for (j in 1:3) {
            if (i != j)
                M[i,j] = exp(M_logits[i,j])
        }
    }

    p1 <- 1 
    if (idx != 3) {
        p2 <- p1 * M[2, 1] 
    } 
    if (idx != 2) {
        p3 <- p1 * M[3, 1] 
    } else {
        p3 <- p2 * M[3, 1] 
    }

    if (idx == 3) {
        p2 <- p3 * M[2, 3] 
    }
    p <- c(p1, p2, p3)
    p / (sum(p))
}

e <- new.env(parent = emptyenv()) 
e[["wlw2"]] = wu2_ld
e[["normal"]] = normal_ld
e[["radial"]] = stratified_ld

e3 <- new.env(parent = emptyenv())
e3[["omit1"]] = function(x) omit1(x, 1)
e3[["omit2"]] = function(x) omit1(x, 2)
e3[["omit3"]] = function(x) omit1(x, 3)


model_binary <- function(dfs, alpha = 0) {
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


    binary = data.frame(n = vector("integer", 0),
                        dataset = vector("character", 0),
                        run = vector("integer", 0),
                        i = vector("integer", 0),
                        j = vector("integer", 0),
                        alpha = vector("numeric", 0),
                        lambda = vector("numeric", 0),
                        lambda_min = vector("numeric", 0),
                        lambda_max = vector("numeric", 0),
                        bacc = vector("numeric", 0),
                        ece = vector("numeric", 0)
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
            # print(lambda)
            model <- glmnet(x, y, family = binomial(), alpha = 0, lambda = lambda, nlambda = 1, standardize = F)
            yt <- as.matrix(select(dft, -class_id))
            r[i,j, ] = predict(model, yt, type = "link")
            r[j,i, ] = -r[i,j,]

            dfb_ij <- filter(dft, class_id %in% c(i,j))
            xb <- as.matrix(select(dfb_ij, -class_id))
            bacc <- mean((dfb_ij$class_id == i)==  (predict(model, xb, type = "link")> 0))
            ece <- getECE(dfb_ij$class_id == i, predict(model, xb, type = "response"))


            binary[nrow(binary) + 1, ] <- data.frame(
                n = n, 
                dataset = dataset,
                run = run,
                i = i, 
                j = j, 
                alpha =alpha,
                lambda = lambda, 
                lambda_min = min(cv1$lambda),
                lambda_max = max(cv1$lambda), 
                bacc = bacc, 
                ece = ece)
        } # end loop j
    } # end loop i 
    return(list(r = r, binary = binary))
}

model_triple <- function(dfs, alpha  = 0) {
    v <- model_binary(dfs, alpha)
    K <- max(v$binary$j)

    res <- data.frame(
            i = vector("integer", 0),
            j = vector("integer", 0),
            k = vector("integer", 0),
            method = vector("character", 0),
            acc = vector("numeric", 0))


    for (i in 1:(K-2)) {
        for (j in (i+1):(K-1)) {
            for (k in (j+1):K) {
                triple <- c(i,j,k)
                truth <- dfs$test$class_id
                r <- v$r[triple, triple, truth %in% triple]
                N <- dim(r)[3]
                map (ls(e), function(m) {
                    p <- sapply(1:N, function(k) {
                        fn <- e[[m]]
                        vecp <- fn(r[,,k])
                        which.max(vecp)
                    })
                    t2 <- truth[truth %in% triple]
                    q <- sapply(1:length(t2), function(k) {
                        which(t2[k] == triple) 
                    })
                    res[nrow(res) + 1,] <<- list(i = i, j = j , k = k, 
                          method = m, 
                          acc = mean(p == q))
                })
                map (ls(e3), function(m) {
                    p <- sapply(1:N, function(k) {
                        fn <- e3[[m]]
                        vecp <- fn(r[,,k])
                        which.max(vecp)
                    })
                    t2 <- truth[truth %in% triple]
                    q <- sapply(1:length(t2), function(k) {
                        which(t2[k] == triple) 
                    })
                    res[nrow(res) + 1,] <<- list(i = i, j = j , k = k, 
                          method = m, 
                          acc = mean(p == q))
                })
            } # loop k
        } # loop j
    } #loop i 
    res
}
 

model_glm1 <- function(dfs, alpha = 0) {
    dft <- dfs$test
    v <- model_binary(dfs, alpha )
    r <- v$r

    multi = map (ls(e), function(m) {
        p <- sapply(1:N, function(k) {
                fn <- e[[m]]
                vecp <- fn(r[,,k])
                which.max(vecp)
                })
        data.frame(n = n, K = K, method = m, dataset = dataset, run = run, lambda = lambda,  correct = sum(p == dft$class_id))
    } )|> list_rbind()
    list(multi = multi, binary = v$binary)
}
