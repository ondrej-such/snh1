source("msvm.R")
# library(glmnet)
library(CalibratR)
suppressPackageStartupMessages(library(MASS))

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
    M <- M_logits

    if (idx != 1) {
        perm <- if (idx == 2) c(2,1,3) else c(3,2,1)
        M <- M[perm, perm]
    }

        if (M[1,2] >= 0) {
            # p1 >= p2
        if (M[1,3] >=0) {
            # p1 is the largets
            p <- c(1, exp(M[2,1]), exp(M[3,1]))
        } else {
            # p3 is the largest
            p <- c(exp(M[1,3]), exp(M[1,3] + M[2,1]), 1)
        }
    } else {
        # p2 >= p1
        if (M[1,3] >= 0) {
            # p2 is the largest
            p <- c(exp(M[1,2]), 1, exp(M[1,2] + M[3,1]))
        } else {
            # p3 is larger than p1, but unsure yet what is the largest
            if (M[1,3] + M[2,1] >= 0) {
                # p2 is the largest
                p <- c(exp(M[1,2]), 1, exp(M[1,2] + M[3,1]))
            } else {
                # p3 is the largest
                p <- c(exp(M[1,3]), exp(M[1,3] + M[2,1]), 1)
            }
        }
    }
    stopifnot(max(p) <= 1)
    r <- p / sum(p)
    if (idx != 1) r[perm] else r
}


omit2 <- function(M_logits, idx) {
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
        p3 <- p2 * M[3, 2] 
    }

    if (idx == 3) {
        p2 <- p3 * M[2, 3] 
    }
    p <- c(p1, p2, p3)
    res <- p / (sum(p))
    if (sum(is.na(res) > 0)) {
            print(M_logits)
            print(idx)
            print(p)
            print(res)
            stopifnot(F)
    }
    res
}

e <- new.env(parent = emptyenv()) 
e[["wlw2"]] = wu2_ld
e[["normal"]] = normal_ld
e[["radial"]] = stratified_ld

e3 <- new.env(parent = emptyenv())
e3[["omit23"]] = function(x) omit1(x, 1)
e3[["omit13"]] = function(x) omit1(x, 2)
e3[["omit12"]] = function(x) omit1(x, 3)

# Function to identify columns constant within groups
remove_constant_within_groups <- function(data, group_col, tol = 0) {
  # Extract the grouping variable
  groups <- data[[group_col]]
  # Get numeric predictors (exclude the group column)
  predictors <- data[, !names(data) %in% group_col, drop = FALSE]
          
  # Check each column for constancy within groups
  constant_cols <- sapply(predictors, function(col) {
  # Split column by group and check if variance is 0 (or near 0) in any group
  # by_group <- tapply(col, groups, function(x) var(x, na.rm = TRUE))
  # any(by_group == 0 | is.na(by_group))  # NA variance occurs if all values are identical
  by_group <- tapply(col, groups, function(x) var(x, na.rm = TRUE))
  sum(by_group <= tol)
                                })
  # Return names of non-constant columns
  names(predictors)[!constant_cols]
}


lda_binary <- function(dfs, subclasses = 1:max(dfs$train$class_id)) {
    df <- dfs$train
    stopifnot(min(df$class_id) == 1)
    pairs <- expand.grid(i = subclasses, j = subclasses) |> 
                dplyr::filter(i < j)
    dft <- dfs$test
    N <- nrow(dft)

    res1 <- map(1:nrow(pairs), function (r) {
        i <- pairs$i[r]
        j <- pairs$j[r]
        # print(c(i,j))

        df_ij <- filter(df, class_id %in% c(i,j))
        df_ij$class_id = factor(df_ij$class_id)

        # Identify non-constant columns
        group_col <- "class_id"
        tol <- 0
        model <- NULL

        while (is.null(model)) {
            non_constant_cols <- remove_constant_within_groups(df_ij, group_col, tol)

            # Subset the data to keep only non-constant columns + group column
            cleaned_data <- df_ij[, c(group_col, non_constant_cols)]


            m1 <- NULL
            m1 <- tryCatch(lda(class_id == i ~ ., 
                                data = cleaned_data), 
            error = function(e) {
                tol <<- if (tol == 0) 1e-8 else 2 * tol
                return(NULL)
            }, 
            finally = if (tol > 0) print(sprintf("Tol increased to %g for %d %d", tol, i,j)))
            # print(is.null(m1))
            model <- m1
        }
        dft <- dfs$test[, c(group_col, non_constant_cols)]
        pred <- predict(model, dft)
        dfb_ij <- filter(dft, class_id %in% c(i,j))
        pred_ij <- predict(model, dfb_ij)

        bacc <- mean((dfb_ij$class_id == i) ==  pred_ij$class)
        ece <- getECE(dfb_ij$class_id == i, pred_ij$posterior[,2])
        pr = predict(model, dft )
        prp <- pr$posterior
        n1 <- sum(cleaned_data$class_id == i)
        n2 <- sum(cleaned_data$class_id == j)
        m1 <- mean(predict(model, cleaned_data)$x[cleaned_data$class_id == i])
        m2 <- mean(predict(model, cleaned_data)$x[cleaned_data$class_id == j])
        v1 <- var(predict(model, cleaned_data)$x[cleaned_data$class_id == i])
        v2 <- var(predict(model, cleaned_data)$x[cleaned_data$class_id == j])
        # print(c(i,j))
        # print(c(m1,m2, (m1*n1 + m2* n2)/ (n1 + n2)))
        r <- log(prp[,2]) - log(prp[,1])
        r1 <- (m1  - m2) * pr$x + (m2^2 - m1^2)/2 + log(n1/n2)
        if (1e-9 < abs(max(r -  r1))) {
            print(c(i,j, abs(max(r - r1))))

        }
        # print(abs(max(prp[,2] - 1 / (1 + exp( - 2* m1 * pr$x)))))
        
        list(i = i, 
             j = j, 
             r = r1,
             # r = r,
             col_used = length(non_constant_cols),
             bacc = bacc,
             ece = ece)
        })

    n <- dfs$n
    dataset <- dfs$dataset
    run <- dfs$run
    K <- length(subclasses)
    R <- array(0, dim = c(K, K, N))

    binary = map(res1, function(li) {
        logis_r <- li$r
        row <- which(li$i == subclasses)
        col <- which(li$j == subclasses)
        R[row, col,] <<- logis_r
        R[col, row, ] <<- -logis_r
        data.frame(n = n, dataset = dataset, run = run, 
                i = li$i, j = li$j, bacc = li$bacc, ece = li$ece, col_used = li$col_used)
    }) |> list_rbind()

    return(list(r = R, binary = binary, subclasses = subclasses))
}

lda_pred3 <- function(dfs, i, j, k) {
    triple <- c(i,j,k)
    v <- lda_binary(dfs, subclasses = triple)
    K <- max(v$binary$j)
    truth <- dfs$test$class_id
    r <- v$r
    idx <- truth %in% triple
    dft <- dfs$test
    N <- nrow(dft)
    t2 <- truth[idx]
    q <- sapply(1:length(t2), function(k) {
            which(t2[k] == triple) 
        })
    map (c(ls(e), ls(e3)), function(m) {
        p <- sapply(1:N, function(k) {
            fn <- e[[m]]
            if (is.null(fn)) fn<- e3[[m]]
            fn(r[,,k])
        }) |> t()
        colnames(p) <- c("p1", "p2", "p3")
        pred <- apply(p[idx,],1, which.max)
        data.frame(p[idx,], 
                    id = 1:sum(idx), orig_id = (1:N)[idx], 
                    truth = q, pred = pred, method = m, 
                    orig_truth = truth[idx])
    }) |>  list_rbind()
    # print(dim(df1))
    # df2 <- map (ls(e3), function(m) {
        # p <- sapply(1:N, function(k) {
            # fn <- e3[[m]]
            # fn(r[,,k])
        # }) |> t()
        # colnames(p) <- c("p1", "p2", "p3")
        # data.frame(p[idx,], truth = q, method = m)
        # }) |> list_rbind()
# 
    # rbind(df1, df2)
}

lda_triples <- function(dfs) {
    v <- lda_binary(dfs)
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

mix_triple <- function(wlws, i, j, k, alpha = 0) {
    v <- lda_binary(dfs, alpha)
    K <- max(v$binary$j)
    triple <- c(i,j,k)
    truth <- dfs$test$class_id
}
 
lda_multi <- function(dfs) {
    dft <- dfs$test
    v <- lda_binary(dfs)
    r <- v$r
    N <- nrow(dft)
    K <- max(dft$class_id)
    run = dfs$run
    n = dfs$n

    # print("binary done")
    multi = map (ls(e), function(m) {
        # print(m)
        p <- sapply(1:N, function(k) {
                fn <- e[[m]]
                vecp <- fn(r[,,k])
                which.max(vecp)
                })
        data.frame(n = n, K = K, method = m, correct = sum(p == dft$class_id))
    } )|> list_rbind()
    # Now, let's do Hinton's oracle
    mr <- sapply(1:N, function(i) 
        sum(r[dft$class_id[i],,i] > 0))
    # print(mr)
    multi <- rbind(multi,
        data.frame(n = n, K = K, method = "oracle", correct = sum(mr == K-1)))

    multi$dataset = dfs$dataset
    multi$run = run
    list(multi = multi, binary = v$binary)
}

write_multi <- function(runs = 20, workers = 11) {
    plan(multicore, workers = workers)
    df <- map(files, function(f) {
        future_map(0:(runs - 1), function (r) {
            dfs <- read_wlws(800, f, r)
            lda_multi(dfs)$multi
        }) |> list_rbind()
    }) |> list_rbind()
    write.csv(df, file = "data/multi.csv")
}
