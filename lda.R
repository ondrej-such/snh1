source("msvm.R")
# library(glmnet)
library(CalibratR)
suppressPackageStartupMessages(library(MASS))
library(dplyr)
suppressPackageStartupMessages(library(tidyr))
library(stringr)

files <- c( "dna", 
            "letter", 
            "mnist", 
            "satimage", 
            "segment", 
            "usps", 
            "waveform")

dataset_names <- function() {
    fs <- dir("unzips/n800", pattern = "*.t-\\d*")
    m1 <- str_match(fs, "(.*).scale.t-(\\d+)")
    unique(m1[,2])
}

make_folds <- function(grp, f = 5) {
    gs <- unique(grp)
    K <- length(gs)
    res <- vector("integer", length(grp))
    for (g in gs) {
        N <- sum(grp == g)
        N5 <- N %/% f
        ubds <- N5 * seq(1, f)
        R5 <- N %% f
        v <- as.integer(rep(0, f))
        if (R5 > 0) {
            pos <- sample(1:f, size = R5, replace = F)
            v[pos] <- 1
        }
        cumv <- cumsum(v)
        cubds <- ubds + cumv
        idx <- grp == g
        perm <- sample(1:N, replace = F, size = N)  - 1
        folds <- sapply(perm, function(x) sum(x >= cubds)) + 1
        res[idx] <- folds
    }
    res
}


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

log_omit1 <- function(M_logits, idx) {
    M <- M_logits
    stopifnot(nrow(M_logits) ==  3)
    stopifnot(nrow(M_logits) ==  3)
    stopifnot(idx >= 1)
    stopifnot(idx <= 3)

    if (idx != 1) {
        perm <- if (idx == 2) c(2,1,3) else c(3,2,1)
        M <- M[perm, perm]
    }

    if (M[1,2] >= 0) {
            # p1 >= p2
        if (M[1,3] >=0) {
            # p1 is the largets
            p <- c(0, M[2,1], M[3,1])
        } else {
# p3 is the largest
            p <- c(M[1,3], M[1,3] + M[2,1], 0)
        }
    } else {
        # p2 >= p1
        if (M[1,3] >= 0) {
            # p2 is the largest
            p <- c(M[1,2], 0, M[1,2] + M[3,1])
        } else {
            # p3 is larger than p1, but unsure yet what is the largest
            if (M[1,3] + M[2,1] >= 0) {
                # p2 is the largest
                p <- c(M[1,2], 0, M[1,2] + M[3,1])
            } else {
                # p3 is the largest
                p <- c(M[1,3], M[1,3] + M[2,1], 0)
            }
        }
    }
    stopifnot(max(p) <= 0)
    if (idx != 1) p[perm] else p
}

omit1 <- function(M_logits, idx) {
    lp <- log_omit1(M_logits, idx)
    p <- exp(lp)
    p / sum(p)
}

half1 <- function(M_logits, idx) {
    s = rep(0, 3)
    for (j in which(idx != c(1,2,3))) {
        s = s +  log_omit1(M_logits, j)
    }
    s = s / 2
    s <- s - max(s)
    p <- exp(s)
    p / sum(p)
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
e3[["half1"]] = function(x) half1(x, 1)
e3[["half2"]] = function(x) half1(x, 2)
e3[["half3"]] = function(x) half1(x, 3)


log3 <- new.env(parent = emptyenv())
log3[["omit23"]] = function(x) log_omit1(x, 1)
log3[["omit13"]] = function(x) log_omit1(x, 2)
log3[["omit12"]] = function(x) log_omit1(x, 3)

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


lda_binary <- function(dfs, subclasses = 1:max(dfs$train$class_id), tol = 0.05) {
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
        tol2 <- 0
        model <- NULL

        while (is.null(model)) {
            non_constant_cols <- remove_constant_within_groups(df_ij, group_col, tol2)
            # Subset the data to keep only non-constant columns + group column
            cleaned_data <- df_ij[, non_constant_cols]
            pca <- prcomp(cleaned_data,          # we do PCA to try to remove collinearity
                                        tol = tol, 
                                        scale. = F,
                                        center = F)
            newd <- as.data.frame(predict(pca, cleaned_data)) # project data to lower dim space
            # print(head(newd))
            # print(class(newd))
            newd$class_id = df_ij$class_id
            # print(dim(newd))
            # print(head(newd))
            # stopifnot(F)

            m1 <- NULL
            m1 <- tryCatch(lda(class_id == i ~ ., 
                                data = newd), 
            error = function(e) {
                tol2 <<- if (tol2 == 0) 1e-8 else 2 * tol2
                return(NULL)
            }, 
            finally = if (tol2 > 0) print(sprintf("Tol2 increased to %e for %d %d", tol2, i,j)))
            # print(is.null(m1))
            model <- m1
        }
        # project using the same projection as in training
        dft <- predict(pca, dfs$test[, non_constant_cols]) |> as.data.frame()
        dft$class_id = dfs$test$class_id
        
        pred <- predict(model, dft)
        dfb_ij <- filter(dft, class_id %in% c(i,j))
        pred_ij <- predict(model, dfb_ij)

        bacc <- mean((dfb_ij$class_id == i) ==  pred_ij$class)
        ece <- getECE(dfb_ij$class_id == i, pred_ij$posterior[,2])
        pr = predict(model, dft )
        prp <- pr$posterior
        n1 <- sum(newd$class_id == i)
        n2 <- sum(newd$class_id == j)
        m1 <- mean(predict(model, newd)$x[newd$class_id == i])
        m2 <- mean(predict(model, newd)$x[newd$class_id == j])
        v1 <- var(predict(model, newd)$x[newd$class_id == i])
        v2 <- var(predict(model, newd)$x[newd$class_id == j])
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
             tol = tol,
             tol2 = tol2, 
             # r = r,
             col_nonconstant = length(non_constant_cols),
             col_used = ncol(dft) - 1,
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
                i = li$i, j = li$j, bacc = li$bacc, ece = li$ece, col_used = li$col_used, tol = li$tol, tol2 = li$tol2)
    }) |> list_rbind()


    return(list(r = R, truth = dfs$test$class_id, summary = binary, subclasses = subclasses))
}

lda_pred3 <- function(dfs, i, j, k) {
    triple <- c(i,j,k)
    v <- lda_binary(dfs, subclasses = triple)
    K <- max(dfs$train$class_id)
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
}

log_pred3 <- function(dfs, i, j, k) {
    triple <- c(i,j,k)
    v <- lda_binary(dfs, subclasses = triple)
    K <- max(dfs$train$class_id)
    truth <- dfs$test$class_id
    r <- v$r
    idx <- truth %in% triple
    dft <- dfs$test
    N <- nrow(dft)
    t2 <- truth[idx]
    q <- sapply(1:length(t2), function(k) {
            which(t2[k] == triple) 
        })
    map (c(ls(log3)), function(m) {
        p <- sapply(1:N, function(k) {
            fn <- log3[[m]]
            # print("done 3")
            fn(r[,,k])
        }) |> t()
        colnames(p) <- c("p1", "p2", "p3")
        pred <- apply(p[idx,],1, which.max)
        data.frame(p[idx,], 
                    id = 1:sum(idx), orig_id = (1:N)[idx], 
                    truth = q, pred = pred, method = m, 
                    orig_truth = truth[idx])
    }) |>  list_rbind()
}

lda_triples <- function(dfs) {
    v <- lda_binary(dfs)
    K <- max(dfs$train$class_id)
    # print(K)

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
                    p <- sapply(1:N, function(s) {
                        fn <- e[[m]]
                        vecp <- fn(r[,,s])
                        if (sum(is.na(vecp)) > 0) {
                            print(sprintf("Triple %d %d %d is bad for %s", 
                                            i, j, k, m))
                            print(r[,,s])
                            print(vecp)
                            stopifnot(F)
                        }
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

multi_pred <- function(dfs, method, tol) {
    dft <- dfs$test
    v <- lda_binary(dfs, tol = tol)
    r <- v$r
    N <- nrow(dft)
    K <- max(dft$class_id)
    run = dfs$run
    n = dfs$n
    fn <- e[[method]]

    probs <- sapply(1:N, function(k) {
            fn(r[,,k])
            }) |> t()
    p <- apply(probs, 1, which.max)
    print(sprintf("Accuracy %f", mean(p == dft$class_id)))
    #data.frame(n = n, K = K, method = m, correct = sum(p == dft$class_id))
    probs
}

score_matrix <- function(p, truth, score) {
    if (nrow(p) != length(truth)) {
        print(dim(p))
        print(length(truth))
        stopifnot(F)

    }

    if (score == "brier") {
        data.factor <- as.factor(truth)
        onehot <- model.matrix(~ data.factor - 1)
        apply(p - onehot, 1, function(x) -sum(x * x))
    } else if (score == "log") {
        idx = cbind(1:length(truth), truth)
        log(p[idx]) # TODO: replac with logsumexp
    } else if (score == "acc") {
        maxes <- apply(p, 1, max)
        mv <- p[cbind(1:nrow(p), truth)]
        cts <- apply(p == maxes,1, sum)
        ifelse(mv == maxes, 1 / cts, 0)
        # as.numeric(apply(p, 1, which.max) == truth)
    } else {
        stopifnot(F)
    }
}
 
lda_multi <- function(dfs, tols, score = "acc") {

   map(tols, function(tol) {
        dft <- dfs$test
        v <- lda_binary(dfs, tol = tol)
        r <- v$r
        N <- nrow(dft)
        K <- max(dft$class_id)
        run = dfs$run
        n = dfs$n

        multi = map (ls(e), function(m) {
            p <- sapply(1:N, function(k) {
                    fn <- e[[m]]
                    fn(r[,,k])
                    }) |> t()
            scores <- score_matrix(p, dft$class_id, score)
            data.frame(n = n, K = K, method = m, correct = mean(scores))
        } )|> list_rbind()
        # Now, let's do Hinton's oracle
        mr <- sapply(1:N, function(i) {
            true_cl <- dft$class_id[i]
            v <- r[,true_cl,i]
            maxv <- max(v)
            v1 <- exp(v - maxv)
            v1 / sum(v1)
        }) |> t()
        print(dim(mr))
        scores <- score_matrix(mr, dft$class_id, score)
        multi <- rbind(multi,
            data.frame(n = n, K = K, method = "Hinton's oracle", correct =
                mean(scores)))
                
#        multi <- rbind(multi,
            #data.frame(n = n, K = K, method = "oracle", correct = sum(mr == K-1)))

        multi$dataset = dfs$dataset
        multi$run = run
        multi$tol = tol
        list(multi = multi, binary = v$summary)
   }) |> my_pmap()
}

write_multi <- function(runs = 20, workers = 11, tol = 1/2^(2:12), score = "acc") {
    plan(multicore, workers = workers)
    df <- map(dataset_names(), function(f) {
        future_map(0:(runs - 1), function (r) {
            dfs <- read_wlws(800, f, r)
            lda_multi(dfs, tols = tol, score = score)$multi
        }) |> list_rbind()
    }) |> list_rbind() |> 
        group_by(dataset, tol, method) |> summarize(acc = mean(correct)) 
# |>
 #       pivot_wider(names_from = "method", values_from = "acc")
    df$score = score
    write.csv(df, file = sprintf("data/multi-%s.csv", score), row.names = F, quote = F)
    df
}

write_triples <- function(runs = 20, workers = 11) {
    plan(multicore, workers = workers)
    df <- map(dataset_names(), files, function(f) {
        future_map(0:(runs - 1), function (r) {
            dfs <- read_wlws(800, f, r)
            df1 <- lda_triples(dfs) |> 
                pivot_wider(names_from = "method", values_from = "acc")
            df1$dataset = f
            df1$run = r
            df1
        }) |> list_rbind()
    }) |> list_rbind() 
    write.csv(df, file = "data/triples.csv", row.names = F, quote = F)
}

# Generate probability distributions on qL classses with denominators = div
#
gen_W <- function(n = 200, div = 100, qL = 4) {
  sapply( 1:n,  function(i) { 
  s <- sample(0:div , qL - 1, replace = T) |> 
       sort()
  if (qL == 2)
      c(s[1], div - s[1]) / div
  else 
      c(s[1], 
      sapply(1:(qL - 2), function(i) (s[i + 1] - s[i])),
      (div - s[qL - 1])) / div
  })  |> t()
}


gen_W2 <- function(div = 50) {
    df <- expand.grid(w1 = 0:div, w2 = 0:div) |>
                filter(w1 + w2 <= div) |> 
                filter((w1 + w2 > 0)) |>
                filter((div - w2) > 0) |>
                filter((div - w1) > 0) 
    W <- matrix(0, nrow = nrow(df), ncol = 3)
    W[,c(1:2)] = as.matrix(df)
    W[,3] = div - W[,1] - W[,2]
    W / div
}

gen_W3 <- function(n = 1000, lbd = 0, ubd = 1) {
    d = ubd - lbd
    matrix((runif(3 * n) + lbd /d) * (ubd - lbd),
        nrow = n, ncol = 3)


}


# use triangular grid (same as gen_W2)
# probably should remove in order not to confuse 
#
gen_W4 <- function(div = 40) {
    M <- matrix(0, nrow = (div + 2) * (div + 1) /2, ncol = 3)
    A <- c(1, 0, 0)
    u <- c(-1, 0, 1) / div 
    v <- c(-1, 1, 0) / div 
    row <- 1
    for (i in 0:div ) {
        for (j in 0: div ) {
            if (i + j > div )
                next
            stopifnot(row <= nrow(M))
            M[row, ] = A + i * u + j * v
            row <- row + 1
        }
    }
    s <- apply(M, 1, sum)
    M / s
}

# Get a list of Bayes-covariant methods's predictions for a triple
#
get_bcp <- function(dfs, i, j, k) {
    pred <- log_pred3(dfs, i, j, k)
    # q1 <- filter(pred, method =="normal") 
    q2 <- filter(pred, method =="omit12") 
    q3 <- filter(pred, method =="omit13") 
    q4 <- filter(pred, method =="omit23") 
    ql <- lapply(list(q2, q3, q4), function(df) 
            as.matrix(dplyr::select(df, 1:3)))

    # ql <- map(ls(log3), function(fnn) {
        # filter(pred, method == fnn) |> as.matrix(dplyr

    # })
    names <- list("omit12", "omit13", "omit23")
    list(names = names, pl = ql, truth  = q2$truth, 
        dataset = dfs$dataset, run = dfs$run, triple = c(i,j,k))
}

# Select a subset of Bayes-covariant methods
# 
select_bcp <- function(bcp, predictors = 1:2) {
    list(pl = bcp$pl[predictors], 
         names = bcp$pl[predictors],
         truth = bcp$truth)
}

# Use cross-validation to estimate accuracy of stacking ensemble
# 
cv_stack <- function(bcp,  W = gen_W(qL = length(bcp$pl)), 
                    f = 5, reps = 20, score = "brier", acc = NULL) {
    if (is.null(acc))
        acc <- eval_scores(bcp, W = W, score = "acc")$scores
    scores <- if(score == "acc") acc else
                eval_scores(bcp, W = W, score = score)$scores
    truth <- bcp$truth
    # print(dim(scores))

    df <- map (1:reps, function(rep) {
        folds <- make_folds(truth, f)
        ws <- map (1:f, function(j) {
            ix = folds != j
            score1 <- apply(scores[ix,], 2, mean)
            stopifnot(length(score1) == nrow(W))
            mscore1 <- max(score1)
            s <- sample(which(score1 == mscore1), 1)
            iy = folds == j
            correct2 <- sum(acc[iy, s])
            acc2 <- mean(acc[iy, s])
            data.frame(rep = rep, best_w = s, score1 = mscore1, acc2 = acc2, correct2 = correct2)
        }) |> list_rbind()
    }) |> list_rbind()
}

cv_acc <- function(...) {
    cv_stack(..., score = "acc")
}

# computes crossvalidation accuracy of stacking ensemble
#
distW <- function(bcp, W = gen_W(qL = length(bcp$pl)), ratio = 0.2, reps = 10, score = "score") {
    # stopifnot(ncol(W) == length(predictors))
    truth <- bcp$truth
    scores <- eval_scores(bcp, W = W, score = score)$scores

    # triple <- c(i, j, k)
    Nsamp <- round(nrow(scores) * ratio)
    res <- vector("integer", nrow(W))
    print(Nsamp)

    for (r in 1:reps) {
        ix = sample(1:nrow(scores), Nsamp)
        score1 <- apply(scores, 2, function(i) sum(scores[ix,i]))
        max_score = max(score1)
        r <- sample(which(max_score == score1), 1)
        res[r] <- res[r] + 1
    }
    print(res)
    cbind(dist = res, W = W)
}

eval_scores <- function(bcp, W = gen_W(qL = length(bcp$pl)), score = "brier") {

    # We use convention from
    # https://www.tandfonline.com/doi/abs/10.1198/016214506000001437
    # in particular, it's better to have higher score

    truth <- bcp$truth # note that this is 1,2, or 3, 
    stopifnot(all(truth %in% c(1,2,3)))
    stopifnot(length(bcp$names) == ncol(W))

    colnames(W) <- bcp$names
    ql <- bcp$pl
    qL <- length(ql)

    # Computation of scores for each sample
    scores <- sapply(1:nrow(W), function(i) {
            logp <- W[i,1] * ql[[1]]
            for (j in 2:qL) {
                logp <- logp + W[i,j] * ql[[j]]
            }
            # cp <- exp(logp)
            m = apply(logp, 1, max)
            p1 <- exp(logp - m)
            sums <- apply(p1, 1, sum)
            p <- p1 / sums
            # p <- apply(logp, 1, function(r) {
                #m = max(r)
                #p1 <- exp(r - m)
                #p1 / sum(p1)
            #}) |> t()
            # if (i == 1) print(dim(p))
            if (sum(is.na(p)) > 0) {
                print(logp)
                print(cp)
                print(p)
                stopifnot(F)
            }
            if (score == "brier") {
                data.factor <- as.factor(truth)
                onehot <- model.matrix(~ data.factor - 1)
                apply(p - onehot, 1, function(x) -sum(x * x))
            } else if (score == "log") {
                idx = cbind(1:length(truth), truth)
                log(p[idx]) # TODO: replac with logsumexp
            } else if (score == "acc") {
                as.numeric(apply(p, 1, which.max) == truth)
            } else {
                stopifnot(F)
            }
        
            } ) # |> t()
    list(scores = scores, W = W, score = score)
}

par_triples <- function(div = 10, score = "acc", workers = 12, limit = 100, seed = 123) {
    df <- read.csv("data/triples.csv")
    # df1 <- df |> mutate(m = pmax(normal, omit12, omit13, omit23)) |> 
    df1 <- df |>  filter(wlw2 > normal)
    print(sprintf("Remains %d out of %d", nrow(df1), nrow(df)))
    # print(head(df1))
    
    W <- gen_W2(div)

    #df2 <- map(1:nrow(df1), function(i) {
    set.seed(seed)
    rows <- if (nrow(df1) > limit) {
        sample(1:nrow(df1), limit)
    } else 1:nrow(df1)
    cv_triples(df1[rows,], W = W, score = score, workers = workers, seed = seed)
}

cv_triples <- function(df1, W, score = "acc", workers = 12, seed = 123) {
    plan(multicore, workers = workers)
    df2 <- future_map(1:nrow(df1), function(i) {
        set.seed(seed + i)
        dfs <- read_wlws(800, df1$dataset[i], df1$run[i]) 
        bcp <- get_bcp(dfs, df1$i[i], df1$j[i], df1$k[i])
        acc <- eval_scores(bcp, W = W, score = "acc")$scores
        df3 <- cv_stack(bcp, W = W, score = score, acc = acc)
        df4 <- df3 |> group_by(rep) |> summarize(m = sum(correct2)) |> 
                    summarize(cv_acc = mean(m) / nrow(bcp$pl[[1]]))
        df5 <- df3 |> group_by(best_w) |> summarize(n = n())|> arrange(desc(n))

        df4$dataset = df1$dataset[i]
        df4$i = df1$i[i]
        df4$j = df1$j[i]
        df4$k = df1$k[i]
        df4$run = df1$run[i]
        df4$wlw2 = df1$wlw2[i]
        df4$best_w = df5$best_w[1]
        df4$best_n = df5$n[1]
        df4$max_bc = max(apply(acc, 2, mean))
        df4
    }) |> list_rbind()
    list(total = nrow(df), remained = nrow(df1) , 
        df = df2) 
}

bc_all <- function(dfs, add = 2) {
    K <- max(dfs$train$class_id)
    triples <- combn(K, 3)

    map(1:ncol(triples), function(m) {
        bc_triple(dfs, triples[1,m], triples[2,m], triples[3,m], add)
    }) |> list_rbind()
}

bc_triple <- function(dfs, i, j, k, add = 2) {
    print(sprintf("%s %d %d %d", dfs$dataset, i, j, k))
    df1 <- dfs$train
    df2 <- dfs$test
    # pred <- lda_pred3(dfs, i, j, k)
    # p <- as.matrix(pred[,1:3])
    K <- max(dfs$train$class_id)

    triple <- c(i,j,k)
    # dfbc <- dfs
    # dfbc$train = dfs$train |> filter(class_id %in% triple)
    # dfbc$test = dfs$test |> filter(class_id %in% triple)

    v <- lda_binary(dfs, subclasses = triple)

    res1 <- map(triple, function(n) {
        {
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
                fn <- e[[m]]
                if (is.null(fn)) fn<- e3[[m]]
                p <- sapply(1:N, function(k) {
                    fn(r[,,k])
                }) |> t()

                # colnames(p) <- c("p1", "p2", "p3")
                w <- rep(1,3)
                w[which(n == triple)] = (1 + add)
                # p <- w * p
                wp <- p
                it <- which(n == triple)     
                wp[, it] = (1 + add) * wp[,it]
                # also could use t(w * t(p))
                predictions <- apply(wp, 1, which.max)
                weight <- if_else (t2 == n, 1 + add, 1)
                # TODO: fix pred$truth
                # correct <- weight * (predictions == pred$truth)
                correct <- weight * (predictions[idx] == q)
                acc <- sum(correct) / sum(weight[idx])

                p2 <- sapply(1:N, function(k) {
                    r2 <- r[,,k]
                    for (s in 1:3) {
                        if (s == it) 
                            next
                        r2[it, s] = r2[it,s] + log((1 + add))
                        r2[s, it] = -r2[it,s]
                    }
                    fn(r2)
                }) |> t()
                # colnames(p2) <- c("p1", "p2", "p3")

                pred <- apply(p2[idx,],1, which.max)
                correct2 <- weight * (pred == q)
                acc2 <- sum(correct2) / sum(weight[idx])
                # acc2 <- mean(pred == q)
                data.frame(changed = n, factor = (1 + add) , binary = acc2, multi = acc, method = m)
                # data.frame(p[idx,], 
                        # id = 1:sum(idx), orig_id = (1:N)[idx], 
                        # truth = q, pred = pred, method = m, 
                    # orig_truth = truth[idx])
            }) |>  list_rbind()
        }
        # ix1 <- dfs$train$class_id == n
        # ix2 <- dfs$test$class_id == n
        # for( w in 1:add) {
            # dfbc$train <- rbind(dfbc$train, dfs$train[ix1,])
            # dfbc$test <- rbind(dfbc$test, dfs$test[ix2,])
        # }
        # stopifnot(sum(dfbc$train$class_id %in% triple) == nrow(dfbc$train))
        # dfbc$train$class_id = sapply(dfbc$train$class_id, function(s) which(s == triple))
        # dfbc$test$class_id = sapply(dfbc$test$class_id, function(s) which(s == triple))
        # dfm <- lda_triples(dfbc)
        # stopifnot(0 == sum(is.na(dfm$acc)))
        # dfm$enlarged = n
        # dfm$bc = "binary"
        # dfm$i = triple[1]
        # # stopifnot(F)
        # dfm$j = triple[2]
        # dfm$k = triple[3]
        # w <- rep(1,3)
        # w[which(n == triple)] = (1 + add)
        # p <- w * p
        # predictions <- apply(p, 1, which.max)
        # weight <- if_else (pred$truth == n, 1 + add, 1)
        # correct <- weight * (predictions == pred$truth)
        # for (method in unique(pred$method))  {
            # ix = pred$method == method
            # acc <- sum(correct[ix]) / sum(weight[ix])
            # dfm[nrow(dfm) + 1, ] = list(i = triple[1], j = triple[2], k = triple[3], 
                    # method = method, acc = acc, 
                    # enlarged = n, bc = "multi")
        # }
        # dfm
    }) |> list_rbind() 
}


# Function to evaluate deviation from Hinton's oracle
#
hinton_triple <- function(dfs, i = 1 , j = 2 ,k = 3) {
    triple <- c(i,j,k)
    print(sprintf("%s %d %d %d", dfs$dataset, i, j, k))
    v <- lda_binary(dfs, subclasses = triple) 
    gt <- matrix(0, nrow = length(v$truth), ncol = 3)
    N <- length(v$truth)
    truth <- sapply(1:N, function(i) which.max(triple == v$truth[i]))
    gt <- sapply(1:N, function(i) {
        omit1(v$r[,,i], truth[i])
    }) |> t()
    # print(dim(gt))
    brier <- matrix(nrow = N, ncol = length(ls(e)))
    acc <- matrix(nrow = N, ncol = length(ls(e)))
    colnames(brier) <- ls(e)
    colnames(acc) <- ls(e)
    for (i in 1:N) {
        gti <- which.max(gt[i,])
        for (m in ls(e)) {
            fn <- e[[m]]
            p <- fn(v$r[,,i])
            brier[i, m] = sum( (p - gt[i,])^2)
            acc[i, m] = which.max(p) == gti
        }
    }
    list(truth = truth, brier = brier, acc = acc)
}

hinton_multi <- function(dfs) {
    v <- lda_binary(dfs)
    K <- length(v$subclasses)
    gt <- matrix(0, nrow = length(v$truth), ncol = K)
    N <- length(v$truth)
    truth <- sapply(1:N, function(i) which.max(v$subclasses == v$truth[i]))

    gt <- sapply(1:N, function(i) {
        # pred <- vector("numeric", K)
        pred <- v$r[,truth[i],i]
        pred <- pred - max(pred)
        p1 <- exp(pred) 
        p1 / sum(p1)
    }) |> t()
    brier <- matrix(nrow = N, ncol = length(ls(e)))
    acc <- matrix(nrow = N, ncol = length(ls(e)))
    colnames(brier) <- ls(e)
    colnames(acc) <- ls(e)
    for (i in 1:N) {
        gti <- which.max(gt[i,])
        for (m in ls(e)) {
            fn <- e[[m]]
            p <- fn(v$r[,,i])
            brier[i, m] = sum( (p - gt[i,])^2)
            acc[i, m] = which.max(p) == gti
        }
    }
    list(truth = truth, brier = brier, acc = acc)
}

hinton <- function(dfs) {
    com <- combn(max(dfs$train$class_id), 3)
    triple <- map(1:ncol(com), function(m) {
        i = com[1,m]
        j = com[2,m]
        k = com[3,m]
        h <- hinton_triple(dfs, i,j,k)
        l1 <- as.list(apply(h$brier, 2, mean))
         #print(l1)
        data.frame(i = i, j = j, k = k,  l1) # l1$normal, l1$radial, l1$wlw2)
    }) |> list_rbind() 
    triple1 <- triple |> dplyr::select(-(1:3)) |> apply(2, mean) 
    multi  <- hinton_multi(dfs)$brier |> apply(2, mean) |> as.list() |> data.frame()
    list(triple = triple1, multi = multi)
}

bc_tuple <- function(dfs, tuple, ns = tuple, add = 2) {
    stopifnot(all(as.logical(sapply(ns, function(n) n %in% tuple))))
    df1 <- dfs$train
    df2 <- dfs$test
    K <- max(dfs$train$class_id)

    v <- lda_binary(dfs, subclasses = tuple)

    truth <- dfs$test$class_id
    r <- v$r
    idx <- truth %in% tuple
    dft <- dfs$test
    N <- nrow(dft)
    t2 <- truth[idx]
    q <- sapply(1:length(t2), function(k) {
            which(t2[k] == tuple) 
        })

    map(ns, function(n) {
        print(n)
        it <- which(n == tuple)     
        
        map (ls(e), function(m) {
            fn <- e[[m]]
            p <- sapply(1:N, function(k) {
                fn(r[,,k])
            }) |> t()

            #w <- rep(1,3)
            #w[which(n == triple)] = (1 + add)
            # p <- w * p
            wp <- p
            old_predictions <- apply(p, 1, which.max)
            wp[, it] = (1 + add) * wp[,it]
            # also could use t(w * t(p))
            predictions <- apply(wp, 1, which.max)
            weight <- if_else (t2 == n, 1 + add, 1)
            correct <- weight * (predictions[idx] == q)
            acc <- sum(correct) / sum(weight[idx])
            flipped = mean(predictions[idx] != old_predictions[idx])

            p2 <- sapply(1:N, function(k) {
                r2 <- r[,,k]
                for (s in 1:3) {
                    if (s == it) 
                        next
                    r2[it, s] = r2[it,s] + log((1 + add))
                    r2[s, it] = -r2[it,s]
                }
                fn(r2)
            }) |> t()

            pred <- apply(p2[idx,],1, which.max)
            correct2 <- weight * (pred == q)
            acc2 <- sum(correct2) / sum(weight[idx])
            data.frame(changed = n, factor = (1 + add) , binary = acc2, multi = acc, method = m, flipped = flipped)
        }) |>  list_rbind()
    }) |>  list_rbind()
}
