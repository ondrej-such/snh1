library(e1071)
suppressPackageStartupMessages(library(dplyr))
library(Rcpp)
library(ParBayesianOptimization)
suppressPackageStartupMessages(library(furrr))
library(future)
library(magrittr)
library(readr)
library(purrr)

sourceCpp("decode.cpp")

# FLAGS:
#   
#   cv2x   - do inner loop CV to get unbiased decision values ?
#   minC   - select (C, sigma) with minimal C and random sigma ?
#   lapl   - do laplace version of logistic fitting ?
#   cent   - use crosscentropy
#   wlw2   - use Wu-Lin-Weng's coupling  ?
#   grid   - use grid search to find optimal C, sigma ? 

read_wlw <- function(f = "unzips/n300/letter.scale-0", min_col = NA, max_col = NA) {
	stopifnot(file.exists(f))
	lines <- readLines(f)
	chunks <- strsplit(lines, " ") 

  # This does not work
	# elts1 <- strsplit(chunks[[1]], ":") 
	# ncol <- length(elts1) - 1

  if (is.na(min_col)) {
    # This assumes that the keys are in increasing order
    lk <- sapply(chunks, function(chunk) {
            elts <- strsplit(chunk, ":")
            # print(elts)
            as.integer((tail(elts, 1)[[1]]))}
            )
    min_col <- min(lk)
    max_col <- max(lk)
  }

  ncol <- max_col - min_col + 1

	M <- matrix(0, nrow = length(lines), ncol = ncol)
	classes <- vector(mode = "integer", length(lines))

	for (i in 1:length(chunks)) {
    # print(chunks[[i]])
    # print(i)
		chunk <- chunks[[i]]
	  elts <- strsplit(chunk, ":")
	    # print(elts)
    keys <- sapply(2:length(elts), function(i) as.numeric(elts[[i]][[1]])) + 1 - min_col
	  vals <- sapply(2:length(elts), function(i) as.numeric(elts[[i]][[2]]))
    # print(keys)
    # print(vals)
	  M[i,keys] <- vals
	  classes[i] <- as.integer(elts[[1]])
	}
	colnames(M) <- sprintf("f%d", 1:ncol)

  # some datasets have classess starting from 0! We need to fix this.
  us <- sort(unique(classes))
  uclasses = vector("integer", length(classes))
  for (i in 1:length(us)) {
       uclasses[classes == us[i]] = i 
  }
	df <- data.frame(class_id = uclasses, M)
  list(df = df, min_col = min_col, ncol = ncol, max_col = max_col)
}

read_wlws <- function(nsamp, dataset, i) {
        name1 <- sprintf("unzips/n%d/%s.scale-%d", nsamp, dataset, i)
        print(name1)
        r1 <- read_wlw(name1)
        name2 <- sub("-", ".t-", name1)
        r2 <- read_wlw(name2)

        if ((r1$min_col != r2$min_col) | (r1$max_col != r2$max_col)) {
            mc = min(r1$min_col, r2$min_col)
            xc = max(1 + r1$ncol - r1$min_col, 1 + r2$ncol - r2$min_col)
            r1 <- read_wlw(name1, min_col = mc, max_col = xc)
            r2 <- read_wlw(name2, min_col = mc, max_col = xc)
        }
        list(train = r1$df, test = r2$df, n = nsamp, dataset = dataset, run = i)
}

sah1 <- function(nsamp = 300 , dataset = "usps", lapl = NA, cent = c(T)) {
    map(0:19, function (i) {
        print(sprintf("%s-%d iter %d", dataset, nsamp, i))
        dfs <- read_wlws(nsamp, dataset, i)
        d1 <- dfs$d1
        d2 <- dfs$d2

        flags <- list(grid = T, minC = T, cv2x = F, cent = cent)
        if (!all(is.na(lapl))) {
            flags$lapl <- lapl
        }
        rs1 <- opt(d1, flags)
        rs1$stage = "search"

        df.factors <- if (all(is.na(lapl)))
                            expand.grid(wlw2 = c(T,F), lapl = c(T,F), cent = cent)
                      else 
                            expand.grid(wlw2 = c(T), lapl = lapl, cent = cent)

        df.eval <- map(1:nrow(df.factors), function (j) {
            res <- subset(rs1, wlw2 == df.factors$wlw2[j] & 
                               lapl == df.factors$lapl[j] &
                               cent == df.factors$cent[j]) |> na.omit()
            stopifnot(nrow(res) > 0)
            max_acc <- max(res$correct)
            rs2 <- subset(res, correct == max_acc)
            minC <- min(rs2$log2C)
            if (nrow(rs2) == 0) {
                print(sum(which(rs2$log2C == minC)))
                print(minC)
                print(dim(rs2))
                print(max_acc)
                print(dim(res))
                print("------")
                print(res)
            }
            ix <- sample(which(rs2$log2C == minC), size = 1)
            flags$log2C = rs2$log2C[ix]
            flags$log2gamma = rs2$log2gamma[ix]
            flags$lapl = rs2$lapl[ix]
            flags$wlw2 = rs2$wlw2[ix]
            flags$cent = rs2$cent[ix]
            dl <- fit_dvs(d1, d2, flags)
            logits <- dvs2logits(dl$old, d1$class_id, dl$new, flags)
            correct <- count_correct(d2$class_id, logits, flags)
            df <- data.frame( lapl = flags$lapl,
                               wlw2 = flags$wlw2,
                            correct = correct,
                            cent = flags$cent,
                            log2C = flags$log2C,
                            log2gamma = flags$log2gamma,
                            acc = correct/nrow(d2),
                            stage = "final")
            df
        }) |> list_rbind()
        rbind(rs1, df.eval) |> mutate(iter = i, samples = nsamp, dataset = dataset)
    }) |> list_rbind()
}

# Selection of the best values of C and gamma
# Heurestic is noted in a foonote of Wu-Lin-Weng paper
# 
opt <- function(data, flags, ... ) {
    stopifnot(!is.null(flags$grid))
    stopifnot(!is.null(flags$minC))
    res <- if (flags$grid) grid_opt(data, flags, ...) else 
                           bayes_opt(data, flags)
                        
    res
    # list(log2C = rs2$log2C[ix], log2gamma = rs2$log2gamma[ix], cv_acc = max_acc)
}

# Grid search for optimal hyperparameters. The grid is the one used in
# https://jmlr.csail.mit.edu/papers/volume5/wu04a/wu04a.pdf
# 
grid_opt <- function(data, flags, workers = 12) {
  stopifnot(flags$grid == TRUE)
	lat <- expand.grid(log2C=seq(-5, 15, by = 2), 
			   log2gamma = seq(-15,15, by = 2))
	lat$method = "grid"
  plan(multicore, workers = workers)

	res <- 1:nrow(lat) %>% future_map( ~ {
    opts <- flags
    opts$log2C <- lat$log2C[.]
    opts$log2gamma <- lat$log2gamma[.]
		data.frame(cv_eval(data = data, opts = opts))
    }
	) %>% future_map_dfr(~ .)
  # cbind(lat, res)
  res
}

# Bayesian optimization of hyperpamaters
# The boundaries of the grid are from 
# https://jmlr.csail.mit.edu/papers/volume5/wu04a/wu04a.pdf
# 
bayes_opt <- function(data, flags, ...) {
    opt_fun <- function(logC, logGamma) {
		   opts <- c(flags, log2C = logC, log2gamma = logGamma)
       list(Score = cv_eval(data, opts))
    }

    br <- ParBayesianOptimization::bayesOpt(opt_fun
           , bounds = list(logC = c(-5,15), logGamma = c(-5,15)) 
           , initPoints = 10
           , iters.n = 137 # so we have 1/3 of evaluations of the grid search
           , gsPoints = 10 #TODO: why 10????
           , verbose = 0)
}

cv_eval <- function(data, opts, ...) {
  stopifnot(!is.null(opts$cv2x))
  if (is.null(data$fold)) {
    # data$fold = sample(1:5, replace = T, size = nrow(data))
    data$fold = make5folds(data$class_id)
  }
	if (opts$cv2x)
		cv2(data = data, opts = opts, ...)
	else
		cv1(data = data, opts = opts, ...)
}



# Creates stratified 5 folds by a given group. 
#
make5folds <- function(grp) {
    gs <- unique(grp)
    K <- length(gs)
    res <- vector("integer", length(grp))
    for (g in gs) {
        N <- sum(grp == g) 
        N5 <- N %/% 5
        ubds <- N5 * seq(1, 5)
        R5 <- N %% 5
        v <- as.integer(rep(0, 5))
        if (R5 > 0) {
            pos <- sample(1:5, size = R5, replace = F)
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


# Two level crossvalidation. Aims to replicate scheme from Figure 4
# in https://jmlr.csail.mit.edu/papers/volume5/wu04a/wu04a.pdf
# 
cv2 <- function(data, opts, ...) {
  stopifnot(F) # not used right now
	correct = 0 
  K = max(data$class_id)
  stopifnot(1 == min(data$class_id))
  stopifnot(!is.null(data$fold))
  res <- expand.grid(lapl = c(F,T), cent = c(T))
  res$correct <- 0
	for (i in 1:5) {
		idx = data$fold == i
		df <- dplyr::filter(data, fold != i) |> dplyr::select(-fold)
    fold2 <- make5folds(df$class_id)

    dvs <- array(dim = c(K, K, nrow(df)))
    for (j in 1:5) {
        df2 <- subset(df, fold2 != j) 
        df2.val <- subset(df, fold2 == j)
        
		    dvs.list <- fit_dvs(data = df2, newdata = df2.val, opts = opts)
        dvs[,,fold2 == j] <- dvs.list$new
    }

		df.val <- dplyr::filter(data, fold == i) |> dplyr::select(-fold)
    dvs.big <- fit_dvs(data = df, newdata = df.val, opts = opts)
    for (l in c(FALSE, TRUE)) {
        opts$lapl = l
        for  (ce in c(TRUE)) {
            opts$cent <- ce
            ab <- get_ab(dvs, df, opts)
            logits <- dvs2logits_ab(ab, dvs.big$new)
            opts$wlw2 <- T
		        correct <- count_correct(df.val$class_id, logits, opts)
            ix <- which(res$lapl == l & res$cent == ce)
            stopifnot(length(ix) == 1)
            res$correct[ix] = res$correct[ix] + correct
        }
    }
	}	
  res |> mutate(acc = correct / nrow(data)) |>
        mutate(error_rate = 1 - acc)
}

dvs2logits_ab <- function(ab, dvs) {
   K = dim(ab$probA)[1]  
   stopifnot(K == dim(ab$probA)[2])
   stopifnot(K == dim(ab$probB)[1])
   stopifnot(K == dim(ab$probB)[2])
   stopifnot(K == dim(dvs)[1])
   stopifnot(K == dim(dvs)[2])

   res <- array(0, dim = dim(dvs))
   for (i in 1:(K-1)) {
       for (j in (i+1):K) {
           res[i,j,] = -(ab$probA[i,j] * dvs[i,j,] + ab$probB[i,j])
           res[j,i,] =  - res[i,j,]
       }
   }
   res
}

# One level crossvalidation. Wu-Lin-Weng seem to suggest this yields biased
# decision values
# 
cv1 <- function(data, opts, ...) {
  stopifnot(is.null(opts$wlw2))
  stopifnot(!is.null(data$fold))

  res <- if (is.null(opts$lapl)) {
      expand.grid(lapl = c(T,F), wlw2 = c(T,F), cent = opts$cent)
  } else  {
      expand.grid(lapl = opts$lapl, wlw2 = T, cent = opts$cent)
  }
  res$correct <- 0
  res$log2C <- opts$log2C
  res$log2gamma <- opts$log2gamma

  for (i in 1:5) {
    df <- subset(data, fold != i) |> dplyr::select(-fold)
    df.val <- subset(data, fold == i) |> dplyr::select(-fold)
    
    dvs.list <- fit_dvs(data = df, newdata = df.val, opts = opts)
    for (l in unique(res$lapl)) {
        opts$lapl = l
        for (cent in unique(res$cent)) {
            opts$cent = cent
            logits <- dvs2logits(dvs.list$old, df$class_id, dvs.list$new, opts)
            for (w in unique(res$wlw2[res$lapl == l])) {
                opts$wlw2 <- w
                new_correct <- count_correct(df.val$class_id, logits, opts)
                idx <- which(res$lapl == l & res$wlw2 == w & res$cent == cent)
                res$correct[idx] <- new_correct + res$correct[idx]
           }
        }
    }	
  }
  res |> mutate(acc = correct / nrow(data))
}

# Evaluates predictions on unseen data using a chosen coupling funciton
# 
count_correct <- function(class_ids, logits, opts) {
    # print("doing count_correct")
    stopifnot(!is.null(opts$wlw2)) # precondition
    cfun <- if (opts$wlw2) wu2_d else stratified_d
    N <- dim(logits)[3]
    stopifnot(N == length(class_ids))
    K <- dim(logits)[1]
    stopifnot(K == dim(logits)[2])
    preds <- matrix(NA, nrow = N, ncol = K)
    for (i in 1:N) {
        preds[i,] <- cfun(logits[,,i])
    }
    cpreds <- apply(preds, 1, which.max)
    sum(cpreds == class_ids)
}

# Fits logistic function to dvs, class_ids data 
# based on the fit, computes logits for array new_dvs decision values 
# 
# Currently not used
dvs2logits <- function(dvs, class_ids, new_dvs, opts) {
  stopifnot(!is.null(opts$lapl)) # precondition
	K = dim(dvs)[1]
	N = dim(dvs)[3]
	stopifnot(N == length(class_ids))
	stopifnot(K == dim(dvs)[2])
  stopifnot(K == dim(new_dvs)[1])
  stopifnot(K == dim(new_dvs)[2])
	res <- array(0, dim = c(K,K, dim(new_dvs)[3]))
	for (i in 1:(K - 1)) {
		for (j in (i + 1):K) {
        targets <- get_targets(i, j, class_ids, opts)
        idx1 <- class_ids == i 
        idx2 <- class_ids == j 
				idx <- idx1 | idx2
        # print(c(sum(idx1), sum(idx2)))
        # print(targets)
				df <- data.frame(target = targets[idx], dv = dvs[i,j,idx])

	      model <- glm(target ~ dv, data = df, family = binomial())
        a <- coef(model)[2] # this is the slope

        # print(opts$cent)
        if ((opts$cent == T) | is.na(a)) {
            pred <- predict(model, data.frame(dv = new_dvs[i,j,]))
				    res[i,j,] = pred
        } else {
            # print('doing Brier')
            fn <- loss_fn(dvs[i,j,idx], targets[idx], opts)
            b <- coef(model)[1] # this is the y-intercept
            # print(summary(model))
            # print(c(a,b))
            ores <- optim(c(a,b), fn, method = "BFGS",
                        control = list(maxit = 250))
            A <- ores$par[1]
            B <- ores$par[2]
            # print(sprintf("slope %.3f %.3f", a, A))
            # print(sprintf("interc %.3f %.3f", b, B))
            # print(ores$counts)
            if (ores$convergence != 0) {
             # stopifnot(F)
               if (fn(c(A,B)) > fn(c(a,b))) {
                    print("reverting to logistic values")
                    A <- a
                    B <- b
               }
            }
            res[i,j,] = A * new_dvs[i,j,] + B
        }
	      res[j,i,] = -res[i,j,]
		}
	}
	res
}

loss_fn <- function(dv, targets, opts) {
	N = length(dv)
  stopifnot(!is.null(opts$cent))

  df <- data.frame(t1 = targets, 
                  t2 = 1 - targets, 
                  dv = dv)
  fn <- if (opts$cent) {
     # This should be crossentropy computation
      function(u) {
          a = u[1]
          b = u[2]

          lin <- a * df$dv + b
          v5 <-  (targets - 1) * lin + log1p(exp(lin))
          v6 <-  targets * lin + log1p(exp(-lin))
          sum(ifelse(lin > 0, v6, v5))
      }
    } else {
        # This is quadratic score (Brier loss)
        # stopifnot(F) # not used currently, because we do not know how to
                     # to optimize reliably
        function(u) {
          a = u[1]
          b = u[2]

          lin <- a * df$dv + b
          p1 <- plogis(lin)
          p2 <- plogis(-lin)
          s1 <- sum(df$t1 * (p1 - p1 * p1 - p2 * p2))
          s2 <- sum(df$t2 * (p2 - p1 * p1 - p2 * p2))
          s1 + s2
        }
    }
    fn
}

get_targets <- function(i, j, class_ids, opts) {
  stopifnot(!is.null(opts$lapl))
  idx1 <- class_ids == j 
  idx2 <- class_ids == i 
  N1 <- sum(idx1)
  N2 <- sum(idx2)
  targets <- vector("numeric", length(class_ids)) 
  idx <- idx1 | idx2 # this  is a logical vector (TRUE if class is either
                     # i or j
  if (is.logical(opts$lapl)) {
      if (opts$lapl) {
        # here we use Bayes prior (Beta(1,1))
        # print("We use Bayes prior")
        targets[idx1] = 1 / (N1 + 2)
        targets[idx2] = (N2 + 1) / (N2 + 2)
      } else {
        # here we use Jeffrey's prior per wikipedia 
        # article on beta distribution
        targets[idx1] = 1 / (2 * (N1 + 1))
        targets[idx2] = (2 * N2 + 1) / (2 * (N2 + 1))
      }
  } else {
      s <- 2^opts$lapl
      targets[idx1] = 1 / (s * N1 + 2)
      targets[idx2] = (1 + s * N2) / (s * N2 + 2)
  }
  targets
}

optim_wlw <- function(dv, ti, fFirst, opts) {
    max_iter <- 100
    sigma <- 1e-12
    eps <- 1e-5
    min_step <- 1e-10
    N1 = sum(fFirst)
    N2 = length(dv) - N1
    N = N1 + N2

    A <- 0
    B <- log((N1 + 1) / (N2 + 1))
    fval <- 0

    fn <- loss_fn(dv, ti, opts)
    fval <- fn(c(A,  B))

    for (iter in 1:max_iter) {
        h11 = sigma
        h22 = sigma
        h21 = 0
        g1 = 0
        g2 = 0
        for (i in 1:N) {
            fApB <- dv[i] *  A + B
            if (fApB > 0) {
                p <- exp(-fApB) / (1 + exp(-fApB))
                q <- 1 / (1 + exp(-fApB))
            } else {
                p <- 1 / (1 + exp(fApB))
                q <- exp(fApB) / (1 + exp(fApB))
            }
            d2 <- p * q
            h11 <- h11 + dv[i] * dv[i] * d2
            h22 <- h22 + d2
            h21 <- h21 + dv[i] * d2
            d1 <- ti[i] - p
            g1 <- g1 + dv[i] * d1    
            g2 <- g2 + d1
        }
        if (abs(g1) < eps & abs(g2) < eps)
            break

        det <- h11 * h22 - h21 * h21
        dA <- -(h22 * g1 - h21 * g2) / det
        dB <- -(-h21 * g1 + h11  * g2) / det
        gd <- g1 * dA + g2 * dB

        stepsize <- 1
        while (stepsize >= min_step) {
            newA <- A + stepsize * dA
            newB <- B + stepsize * dB
            newf <- fn(c(newA, newB))
            if (newf < fval + 0.0001 * stepsize * gd) {
                A <- newA
                B <- newB
                fval <- newf
                # print(c(A,B))
                # print(fval)
                break
            }
            stepsize <- stepsize / 2
        }
        stopifnot (stepsize >= min_step)
    }
    stopifnot(iter < max_iter)
    # print(sprintf("%d iterations", iter))
    return (c(A,B))
}

get_ab <- function(dvs, data, opts) {
	ids <- as.integer(unique(data$class_id))
  stopifnot(1 == min(ids))
	K <- length(ids)
	N = dim(dvs)[3]
  stopifnot(K == max(ids))
  stopifnot(!is.null(opts$lapl))
  stopifnot(!is.null(opts$cent))
  stopifnot(N == nrow(data))
  stopifnot (opts$cent) 

  targets <- vector("numeric", N)
  probA <- matrix(0, nrow = K, ncol = K)
  probB <- matrix(0, nrow = K, ncol = K)
  class_ids <- data$class_id
  mit <- 0
  n_converged <- 0

	for (i in 1:(K - 1)) {
		for (j in (i + 1):K) {
        idx1 <- class_ids == j 
        idx2 <- class_ids == i 
        idx <- idx1 | idx2 # this  is a logical vector 

        targets <- get_targets(i, j, data$class_id, opts)
        wo <- optim_wlw(dvs[i,j,idx], targets[idx], data$class_id[idx] == i, opts)

        probA[i,j] = wo[1]
        probB[i,j] = wo[2]

        probA[j,i] = -probA[i,j]
        probB[j,i] = -probB[i,j]
		}
	}
  return(list(probA = probA, probB = probB, maxit = mit)) 
}

# Creates SVM model for given parameters C, sigma and evaluates the resulting 
# decision values on newdata dataset
# 
fit_dvs <- function(data, newdata, opts) {
	ids <- unique(data$class_id)	
  stopifnot(1 == min(ids))
	K <- length(ids)
  stopifnot(K == max(ids))
	dvs.new <- array(dim = c(K, K, nrow(newdata)))
	dvs.old <- array(dim = c(K, K, nrow(data)))

  df.1 <- data
  df.1$class_id = factor(df.1$class_id)
  df.2 <- df.1 |> droplevels()
  probA <- matrix(0, nrow = K, ncol = K)
  probB <- matrix(0, nrow = K, ncol = K)

  model <- svm(class_id ~ . , 
      data = df.2,
      kernel = "radial", 
      scale = FALSE,
      probability = TRUE,
      gamma = 2^opts$log2gamma,
      cost = 2^opts$log2C)
  
  stopifnot(model$type == 0) # type should be C-classification
  pred <- predict(model, newdata, decision.values = T)
  dv <- attr(pred, "decision.values")
  for (k in 1:ncol(dv)) {
      keys <- strsplit(colnames(dv)[k], "/") [[1]]
      i <- as.integer(keys[[1]])
      j <- as.integer(keys[[2]])
      dvs.new[i,j,] = dv[,k]
      dvs.new[j,i,] = dv[,k]
  }

  pred <- predict(model, data, decision.values = T)
  dv <- attr(pred, "decision.values")
  for (k in 1:ncol(dv)) {
      keys <- strsplit(colnames(dv)[k], "/")[[1]]
      i <- as.integer(keys[[1]])
      j <- as.integer(keys[[2]])
      dvs.old[i,j,] = dv[,k]
      dvs.old[j,i,] = dv[,k]
      probA[i,j] <- model$probA[k]
      probA[j,i] <- -probA[i,j]
      probB[i,j] <- model$probB[k]
      probB[j,i] <- -probB[i,j]
  }
     
  list(new = dvs.new, old = dvs.old, probA = probA, probB = probB)
}

