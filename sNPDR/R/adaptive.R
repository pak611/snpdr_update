# =========================================================================#
#' vwok
#'
#' Variable-Wise Optimized k
#' method for optimizing NPDR scores for each attribute as a function of k
#' Computes p x k beta and P value matrices for a data set with p attributes
#'
#' @param dats m x (p+1) data set of m instances and p attributes with 1 binary outcome or m x [p(p - 1) + 1] with p(p-1) correlations and 1 outcome. Outcome is last column for standard m x (p + 1) and first column for m x [p(p - 1) + 1] (no good reason for the difference).
#' @param k.grid increasing sequence of k values used as looping index. Default is seq(1,(nrow(dats)-1),by=1).
#' @param verbose logical indicating whether to print progress with loop. Default is FALSE, but TRUE also does not give anything useful.
#' @param attr.diff.type character indicating the type of attribute diff to use. Default is 'numeric-abs' for standard continuous data. Use 'correlation-data' for rs-fMRI data.
#' @param corr.attr.names character indicating names of ROIs for attr.diff.type='correlation-data'. Default is NULL.
#' @param signal.names variable names with simulated signals
#' @param separate.hitmiss.nbds logical indicating whether to compute hit/miss neighborhoods separately. Default is FALSE.
#' @param label character indicating type of response. Default is "class" and should not change as of yet.
#'
#' @return A list with:
#' \describe{
#'   \item{vwok.out}{p x 4 data.frame of sorted beta coefficients, atts, vwok ks, and p-values from NPDR}
#'   \item{best.auPRC.k}{1 x 2 data.frame of auPRC-optimal fixed k (if signal.names provided) and corresponding auPRC}
#' }
#'
#' @examples
#' dats <- do.call(rbind, case.control.3sets[c("train", "holdout")])
#'
#' # run Variable-Wise Optimized k function
#' \dontrun{
#' out <- vwok(
#'   dats = dats,
#'   k.grid = NULL,
#'   verbose = T,
#'   attr.diff.type = "numeric-abs",
#'   label = "class"
#' )
#' }
#'
#' @export
#'
vwok <- function(dats = NULL,
                 k.grid = NULL,
                 verbose = F,
                 attr.diff.type = "numeric-abs",
                 corr.attr.names = NULL,
                 signal.names = NULL,
                 separate.hitmiss.nbds = FALSE,
                 label = "class") {
  check_installed("foreach", reason = "for fast parallel computing with `foreach()` and `%dopar%`")
  check_installed("doParallel", reason = "for `registerDoParallel()`")
  check_installed("parallel", reason = "for `makeCluster()`, `detectCores()`, and `stopCluster()`")
  `%dopar%` <- foreach::`%dopar%`

  start_time <- Sys.time()

  m <- nrow(dats)
  p <- ncol(dats) - 1

  num.attr <- if (attr.diff.type == "correlation-data") { # corrdata
    ceiling(sqrt(p))
  } else {
    ncol(dats) - 1
  }

  if (attr.diff.type == "correlation-data") { # corrdata
    attr.idx.list <- list()
    for (i in seq.int(num.attr)) {
      lo.idx <- (i - 1) * (num.attr - 1) + 1
      hi.idx <- i * (num.attr - 1)
      attr.idx.list[[i]] <- seq(lo.idx, hi.idx)
    }
  }

  ks <- k.grid %||% seq.int((m - 1)) # if not separating hit/miss

  if (attr.diff.type == "correlation-data") {
    pheno.vec <- as.factor(dats[, 1])
    dats <- as.matrix(dats[, -1])
  }

  avai.cors <- parallel::detectCores() - 2
  cl <- parallel::makeCluster(avai.cors)
  doParallel::registerDoParallel(cl)

  beta.mat <- matrix(0, nrow = num.attr, ncol = length(ks))
  pval.mat <- matrix(0, nrow = num.attr, ncol = length(ks))
  out.list <- foreach::foreach(
    k = ks,
    .export = c("npdr"),
    .packages = c("npdr")
  ) %dopar% {

    # cat("k = ",k,"\n")

    if (attr.diff.type != "correlation-data") {
      npdr.cc.results <- npdr::npdr(label, dats,
        regression.type = "binomial", attr.diff.type = "numeric-abs",
        nbd.method = "relieff", nbd.metric = "manhattan", msurf.sd.frac = .5, knn = k,
        neighbor.sampling = "none", separate.hitmiss.nbds = separate.hitmiss.nbds,
        padj.method = "bonferroni", verbose = verbose
      )
    } else {
      # pheno.vec <- dats[,1]
      # dats <- dats[,-1]
      npdr.cc.results <- npdr::npdr(pheno.vec, dats,
        regression.type = "binomial", attr.diff.type = "correlation-data",
        nbd.method = "relieff", nbd.metric = "manhattan", msurf.sd.frac = .5, knn = k,
        neighbor.sampling = "none", separate.hitmiss.nbds = separate.hitmiss.nbds,
        padj.method = "bonferroni", verbose = verbose,
        corr.attr.names = corr.attr.names,
        dopar.nn = T, dopar.reg = T
      )
    }

    out.mat <- npdr.cc.results %>%
      select(att, beta.Z.att, pval.adj) %>%
      mutate(att = as.character(att)) %>%
      arrange(att)

    beta.mat <- as.numeric(as.character(out.mat[, "beta"]))
    pval.mat <- as.numeric(as.character(out.mat[, "pval"]))
    att.mat <- as.character(out.mat[, 1])

    if (!is.null(signal.names)) {
      functional.vars <- signal.names

      idx.func <- which(c(att.mat %in% functional.vars))
      func.betas <- beta.mat[idx.func] # functional variable betas
      neg.betas <- beta.mat[-idx.func] # noise variable betas

      # precision-recall curve and area
      pr.npdr <- PRROC::pr.curve(
        scores.class0 = func.betas,
        scores.class1 = neg.betas,
        curve = T
      )

      auPRC <- pr.npdr$auc.integral # area under the precision-recall curve

      list(
        betas = beta.mat,
        pvals = pval.mat,
        atts = att.mat,
        auPRC = auPRC
      )
    } else {
      list(
        betas = beta.mat,
        pvals = pval.mat,
        atts = att.mat
      )
    }
  }
  parallel::stopCluster(cl)

  avai.cors <- parallel::detectCores() - 2
  cl <- parallel::makeCluster(avai.cors)
  doParallel::registerDoParallel(cl)

  beta.mat <- NULL
  out.betas <- foreach::foreach(k = 1:length(ks), .combine = "cbind") %dopar% {
    out.list[[k]]$betas
  }
  parallel::stopCluster(cl)

  avai.cors <- parallel::detectCores() - 2
  cl <- parallel::makeCluster(avai.cors)
  doParallel::registerDoParallel(cl)

  pval.mat <- NULL
  out.pvals <- foreach::foreach(k = 1:length(ks), .combine = "cbind") %dopar% {
    out.list[[k]]$pvals
  }
  parallel::stopCluster(cl)

  if (!is.null(signal.names)) {
    avai.cors <- parallel::detectCores() - 2
    cl <- parallel::makeCluster(avai.cors)
    doParallel::registerDoParallel(cl)

    auPRC.mat <- NULL
    out.auPRC <- foreach::foreach(k = 1:length(ks), .combine = "rbind") %dopar% {
      matrix(c(k, out.list[[k]]$auPRC), nrow = 1, ncol = 2)
    }
    parallel::stopCluster(cl)
  }

  beta.mat <- as.data.frame(out.betas) %>%
    `colnames<-`(paste0("k.", ks)) %>%
    `row.names<-`(as.character(out.list[[1]]$atts))

  pval.mat <- as.data.frame(out.pvals) %>%
    `colnames<-`(paste0("k.", ks)) %>%
    `row.names<-`(as.character(out.list[[1]]$atts))

  if (!is.null(signal.names)) {
    auPRC.mat <- as.data.frame(out.auPRC) %>%
      `colnames<-`(c("k", "auPRC"))
  }

  betas <- beta.mat # num.variables x (num.samples - 1) matrix of betas
  pvals <- pval.mat # num.variables x (num.samples - 1) matrix of p-values
  best.ks <- apply(betas, 1, which.max) # best k's computed from max beta for each attribute
  best.betas <- numeric() # betas corresponding to best k's
  best.pvals <- numeric() # p-values corresponding to best k's

  for (j in seq.int(nrow(betas))) {
    best.betas[j] <- betas[j, as.numeric(best.ks[j])]
    best.pvals[j] <- pvals[j, as.numeric(best.ks[j])]
  }

  # data frame of attributes, betas, and p-values from vwak
  df.betas <- data.frame(
    att = row.names(betas),
    best.ks = best.ks,
    betas = best.betas,
    pval.att = best.pvals
  ) %>%
    # sort data frame by decreasing beta
    arrange(desc(betas))

  if (!is.null(signal.names)) {
    best.auPRC.k <- as.data.frame(matrix(auPRC.mat[which.max(auPRC.mat[, "auPRC"]), ], nrow = 1, ncol = 2)) %>% 
      `colnames<-`(c("k", "auPRC"))
  }

  end_time <- Sys.time()
  elapsed <- end_time - start_time
  cat("Elapsed: ", elapsed, "\n")

  if (!is.null(signal.names)) {
    return(list(vwok.out = df.betas, best.auPRC.k = best.auPRC.k))
  }

  list(vwok.out = df.betas)
}
