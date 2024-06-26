# =========================================================================#
#' npdrDiff
#'
#' A diff is a function that computes the difference of values for an attribute between two instances.
#' It is used for attribute selection for attribute diffs and phenotype diffs.
#' This function is vectorized: input a and b can be two vectors of values for one attribute.
#'
#' @param a value of attribute for first instance. Vector for correlation-data.
#' @param b value of attribute for second instance. Vector for correlation-data.
#' @param diff.type Metric for the difference computation.
#' @param norm.fac Normalization factor.
#'
#' @return Value or vector of differences between two vectors element-wise.
#' @export
npdrDiff <- function(a, b, diff.type = c("manhattan", "numeric-abs", "numeric-sqr", "allele-sharing", "match-mismatch", "correlation-data", "binomial-surv"), norm.fac = 1) {
  if (is.null(diff.type)) diff.type <- "numeric-abs"
  diff.type <- match.arg(diff.type)

  #browser()
  
  val <- switch(diff.type,
                `numeric-sqr` = abs(a - b)^2 / norm.fac, # numeric squared difference
                `allele-sharing` = abs(a - b) / 2, # snps
                `match-mismatch` = ifelse(a == b, 0, 1), # hit pairs = 0, miss = 1
                `correlation-data` = rowSums(abs(a - b) / norm.fac), # a and b are matrices
                `numeric-abs` = abs(a - b) / norm.fac, # numeric abs difference
                `manhattan` = abs(a - b) / norm.fac, # same as numeric-abs
                `binomial-surv` = ifelse(a[,2] == b[ ,2] & a[,1] == b[,1], 1, 0) # hit pairs = 0, miss = 1
            
  )
  # For correlation data, a and b are matrices
  # with m*k rows and numvars-1 cols.
  # m*k rows because looking at all neighbor pairs
  # (fixed k not required).
  # nvars-1 because for a given var,
  # we are looking at all other correlation partners.
  # a represents the first of neighbor pairs
  # b represents the second of neighbor pairs
  # See Eq. 157 and Fig. 9 from
  # https://doi.org/10.1371/journal.pone.0246761
  
  val
}

# =========================================================================#
#' npdrDistances
#'
#' Create m x m distance matrix from m instances and p attributes using different metrics. Used by nearestNeighbors().
#' Note: Probably best to standardize data before manhattan and euclidean.
#'
#' @param attr.mat m x p matrix of m instances and p attributes
#' @param metric for distance matrix between instances (default: \code{"manhattan"},
#' others include \code{"euclidean"}, versions scaled by max-min,
#' \code{"relief-scaled-manhattan"} and \code{"relief-scaled-euclidean"},
#' and for GWAS \code{"allele-sharing-manhattan"}).
#' @param fast.dist whether or not distance is computed by faster algorithm in wordspace, default as F
#'
#' @return  matrix of m x m (instances x intances) pairwise distances.
#' @export
#' @examples
#' train_dat <- case.control.3sets$train
#' dist.mat <- npdrDistances(
#'   train_dat[, names(train_dat) != "class"],
#'   metric = "manhattan"
#' )
npdrDistances <- function(attr.mat, metric = "manhattan", fast.dist = FALSE) {
  npdr.dist.fn <- dist
  if (fast.dist) {
    check_installed("wordspace", reason = "for fast distance computation with `dist.matrix()`")
    npdr.dist.fn <- wordspace::dist.matrix
  }
  
  # Compute distance matrix between all samples (rows)
  # default is numeric manhattan ("manhattan"), max-min scaling is only needed for relief
  if (metric == "hamming") {
    return(as.matrix(hamming.binary(attr.mat)))
  }
  
  if (metric == "allele-sharing-manhattan") {
    # allele-sharing-manhattan, AM for SNPs
    
    attr.mat.scale <- attr.mat / 2
    distance.mat <- npdr.dist.fn(attr.mat.scale, method = "manhattan")
    return(as.matrix(distance.mat))
  }
  
  if (metric == "relief-scaled-euclidean" ||
      metric == "relief-scaled-manhattan") {
    # value of metric, euclidean, manhattan or maximum
    
    method <- strsplit(metric, "-")[[1]][3]
    maxminVec <- attr.range(attr.mat)
    minVec <- apply(attr.mat, 2, min)
    attr.mat.centered <- t(attr.mat) - minVec
    attr.mat.scale <- t(attr.mat.centered / maxminVec)
    distance.mat <- npdr.dist.fn(attr.mat.scale, method = method)
    return(as.matrix(distance.mat))
  }
  
  if (metric == "euclidean" || metric == "manhattan") {
    distance.mat <- npdr.dist.fn(attr.mat, method = metric)
    return(as.matrix(distance.mat))
  }
}

# =========================================================================#
#' nearestNeighbors
#'
#' Find nearest neighbors of each instance using relief.method
#' Used for npdr (no hits or misses specified in neighbor function).
#'
#' @param attr.mat m x p matrix of m instances and p attributes, or m x m distance matrix if metric is `precomputed`
#' @param nbd.metric used in npdrDistances for distance matrix between instances. parameter can be `precomputed` if user-supplied distance matrix.
#' default to `manhattan` (numeric input).
#' @param nbd.method neighborhood method `multisurf` or `surf` (no k) or `relieff` if you want to specify k (required k).
#' @param sd.vec vector of standard deviations for surf
#' @param sd.frac multiplier of the standard deviation from the mean distances, subtracted from mean distance to create for SURF or multiSURF radius. The multiSURF default "dead-band radius" is sd.frac=0.5: mean - sd/2
#' @param k number of constant nearest hits/misses for \code{"relieff"} (fixed k).
#' The default k=0 means use the expected SURF theoretical k with sd.frac (.5 by default) for relieff nbd.
#' @param neighbor.sampling "none" or \code{"unique"} if you want to return only unique neighbor pairs
#' @param att_to_remove attributes for removal (possible confounders) from the distance matrix calculation.
#' @param fast.dist whether or not distance is computed by faster algorithm in wordspace, default as F
#' @param dopar.nn whether or not neighborhood is computed in parallel, default as F
#' @return  Ri_NN.idxmat, matrix of Ri's (first column) and their NN's (second column)
#'
#' @export
#' @examples
#'
#' # multisurf neighborhood with sigma/2 (sd.frac=0.5) "dead-band" boundary
#' neighbor.pairs.idx <- nearestNeighbors(
#'   predictors.mat,
#'   nbd.method = "multisurf",
#'   nbd.metric = "manhattan",
#'   sd.frac = 0.5
#' )
#' head(neighbor.pairs.idx)
#'
#' # reliefF (fixed-k) neighborhood using default `k` equal to
#' # theoretical surf expected value.
#' # One can change the theoretical value by changing sd.frac (default 0.5).
#' neighbor.pairs.idx <- nearestNeighbors(
#'   predictors.mat,
#'   nbd.method = "relieff",
#'   nbd.metric = "manhattan"
#' )
#' head(neighbor.pairs.idx)
#'
#' # reliefF (fixed-k) neighborhood with a user-specified k
#' neighbor.pairs.idx <- nearestNeighbors(
#'   predictors.mat,
#'   nbd.method = "relieff",
#'   nbd.metric = "manhattan",
#'   k = 10
#' )
#' head(neighbor.pairs.idx)
nearestNeighbors <- function(attr.mat,
                             nbd.method = "multisurf",
                             nbd.metric = "manhattan",
                             sd.vec = NULL, sd.frac = 0.5, k = 0,
                             neighbor.sampling = "none",
                             att_to_remove = c(), 
                             fast.dist = FALSE, dopar.nn = FALSE) {
  
  #browser()
  if (dopar.nn) {
    check_installed("foreach", reason = "for fast parallel computing with `foreach()` and `%dopar%`")
    check_installed("doParallel", reason = "for `registerDoParallel()`")
    check_installed("parallel", reason = "for `makeCluster()`, `detectCores()`, and `stopCluster()`")
    `%dopar%` <- foreach::`%dopar%`
  }
  
  # Goal is to create a matrix with num.samp rows and two columns.
  # First column is sample Ri, second is Ri's nearest neighbors (NN)
  # this is basically an edge list.
  
  if (nbd.metric == "precomputed"){
    # allow user to input their own distance matrix
    # in which case attr.mat will be a distance matrix
    num.samp <- nrow(attr.mat)
    dist.mat <- attr.mat %>%  as.data.frame() %>%
      `colnames<-`(seq.int(num.samp))
  } else{
    # if distance is not precomputed, use a metric
    num.samp <- nrow(attr.mat)
    if (!is.null(att_to_remove)) {
      # remove attributes (possible confounders) 
      #from distance matrix calculation
      tryCatch(
        attr.mat <- attr.mat %>% data.frame() %>%
          select(-att_to_remove),
        error = function(c) "The attribute to remove does not exist."
      )
    }
    dist.mat <- attr.mat %>%
      as.matrix() %>%
      unname() %>%
      npdrDistances(metric = nbd.metric, fast.dist = fast.dist) %>%
      as.data.frame() %>%
      `colnames<-`(seq.int(num.samp))
  } # end if precomputed
  # now find neighbors
  if (nbd.method == "relieff") {
    if (k == 0) { # if no k specified or value 0
      # replace k with the theoretical expected value for SURF (close to multiSURF)
      erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
      # theoretical surf k (sd.frac=.5) for regression problems (does not depend on a hit/miss group)
      k <- floor((num.samp - 1) * (1 - erf(sd.frac / sqrt(2))) / 2) # uses sd.frac
    }
    
    if (dopar.nn) {
      avai.cors <- parallel::detectCores() # - 2
      #cl <- parallel::makeCluster(avai.cors)
      #doParallel::registerDoParallel(cl)
      doParallel::registerDoParallel(cores=avai.cors)
      Ri_NN.idxmat <- foreach::foreach(
        Ri.int = seq.int(num.samp), .combine = "rbind", .packages = c("dplyr")
      ) %dopar% {
        Ri <- as.character(Ri.int)
        Ri.nearest.idx <- get_Ri_nearest(dist.mat, Ri, nbd.method, k = k)
        
        return(data.frame(Ri_idx = Ri.int, NN_idx = Ri.nearest.idx))
      }
      ## [1] 1.000000 1.414214 1.732051
      #parallel::stopCluster(cl)
    } else {
      Ri.nearestPairs.list <- vector("list", num.samp)
      for (Ri in colnames(dist.mat)) { # for each sample Ri
        Ri.int <- as.integer(Ri)
        Ri.nearest.idx <- get_Ri_nearest(dist.mat, Ri, nbd.method, k = k)
        
        if (!is.null(Ri.nearest.idx)) { # if neighborhood not empty
          # bind automatically repeated Ri, make sure to skip Ri self
          Ri.nearestPairs.list[[Ri.int]] <- data.frame(Ri_idx = Ri.int, NN_idx = Ri.nearest.idx)
        }
      }
      Ri_NN.idxmat <- bind_rows(Ri.nearestPairs.list)
    }
  } else {
    if (nbd.method == "surf") {
      num.pair <- num.samp * (num.samp - 1) / 2 # number of paired distances
      radius.surf <- sum(dist.mat) / (2 * num.pair) # const r = mean(all distances)
      sd.const <- sd(dist.mat[upper.tri(dist.mat)])
      # bam: orignal surf does not subtract sd-frac but should for fair multisurf comparison
      Ri.radius <- rep(radius.surf - sd.frac * sd.const, num.samp)
      names(Ri.radius) <- as.character(1:num.samp)
    }
    if (nbd.method == "multisurf") {
      
      if (is.null(sd.vec)) {
        sd.vec <- sapply(1:num.samp, function(x) sd(dist.mat[-x, x]))
        
      }
      Ri.radius <- colSums(dist.mat) / (num.samp - 1) - sd.frac * sd.vec # use adaptive radius
    }
    if (dopar.nn) {
      print('here1')
      avai.cors <- parallel::detectCores() - 2
      cl <- parallel::makeCluster(avai.cors)
      doParallel::registerDoParallel(cl)
      Ri_NN.idxmat <- foreach::foreach(
        Ri.int = seq.int(num.samp), .combine = "rbind", .packages = c("dplyr")
        
      ) %dopar% {
        Ri <- as.character(Ri.int)
        Ri.nearest.idx <- get_Ri_nearest2(dist.mat, Ri, nbd.method, Ri.radius = Ri.radius)
        #browser()
        
        return(data.frame(Ri_idx = Ri.int, NN_idx = Ri.nearest.idx))
      }
      #parallel::stopCluster(cl)
    } else {
      print('here2')
      # put each Ri's nbd in a list then rbind them at the end with bind_rows()
      Ri.nearestPairs.list <- vector("list", num.samp) # initialize list
      
      for (Ri in colnames(dist.mat)) { # for each sample Ri
        Ri.int <- as.integer(Ri)
        Ri.nearest.idx <- get_Ri_nearest2(dist.mat, Ri, nbd.method, k, Ri.radius = Ri.radius)
        
        if (!is.null(Ri.nearest.idx)) { # similar to relieff
          Ri.nearestPairs.list[[Ri.int]] <- data.frame(Ri_idx = Ri.int, NN_idx = Ri.nearest.idx)
        }
      }
      Ri_NN.idxmat <- bind_rows(Ri.nearestPairs.list)
    }
  }
  
  
  if (neighbor.sampling == "unique") {
    # if you only want to return unique neighbors
    Ri_NN.idxmat <- uniqueNeighbors(Ri_NN.idxmat)
  }
  
  # matrix of Ri's (first column) and their NN's (second column)
  Ri_NN.idxmat
}

# =========================================================================#
#' nearestNeighborsSeparateHitMiss
#'
#' Find nearest neighbors of each instance using relief.method.
#' Treat the hit and miss distributions separately to circument potential hit bias.
#' ReliefF version makes hit/miss neighborhoods balanced. Surf and MultiSurf are still imbalanced.
#' Used for npdr (no hits or misses specified in neighbor function).
#'
#' @param attr.mat m x p matrix of m instances and p attributes
#' @param pheno.vec vector of class values for m instances
#' @param nbd.metric used in npdrDistances for distance matrix between instances, default: \code{"manhattan"} (numeric)
#' @param nbd.method neighborhood method \code{"multisurf"} or \code{"surf"} (no k) or \code{"relieff"} (specify k)
#' @param k number of constant nearest hits/misses for \code{"relieff"} (fixed k).
#' The default k=0 means use the expected SURF theoretical k with sd.frac (.5 by default) for relieff nbd.
#' @param sd.frac multiplier of the standard deviation from the mean distances, subtracted from mean distance to create for SURF or multiSURF radius. The multiSURF default "dead-band radius" is sd.frac=0.5: mean - sd/2
#' @param neighbor.sampling "none" or \code{"unique"} if you want to return only unique neighbor pairs
#' @param att_to_remove attributes for removal (possible confounders) from the distance matrix calculation.
#' @param fast.dist whether or not distance is computed by faster algorithm in wordspace, default as F
#' @param dopar.nn whether or not neighborhood is computed in parallel, default as F
#' @return  Ri_NN.idxmat, matrix of Ri's (first column) and their NN's (second column)
#'
#' @export
#' @examples
#' # reliefF (fixed-k) neighborhood using default k equal to theoretical surf expected value
#' # One can change the theoretical value by changing sd.frac (default 0.5)
#' neighbor.pairs.idx <- nearestNeighborsSeparateHitMiss(
#'   predictors.mat, case.control.3sets$train$class, # need attributes and pheno
#'   nbd.method = "relieff", nbd.metric = "manhattan",
#'   sd.frac = .5, k = 0
#' )
#' head(neighbor.pairs.idx)
nearestNeighborsSeparateHitMiss <- function(attr.mat, pheno.vec,
                                            nbd.method = "relieff",
                                            nbd.metric = "manhattan",
                                            sd.frac = 0.5, k = 0,
                                            neighbor.sampling = "none",
                                            att_to_remove = c(), fast.dist = FALSE, dopar.nn = FALSE) {
  if (dopar.nn) {
    check_installed("foreach", reason = "for fast parallel computing with `foreach()` and `%dopar%`")
    check_installed("doParallel", reason = "for `registerDoParallel()`")
    check_installed("parallel", reason = "for `makeCluster()`, `detectCores()`, and `stopCluster()`")
    `%dopar%` <- foreach::`%dopar%`
  }
  
  # create a matrix with num.samp rows and two columns
  # first column is sample Ri, second is Ri's nearest neighbors
  num.samp <- nrow(attr.mat)
  pheno.vec <- as.numeric(as.character(pheno.vec))
  majority.pheno <- which.max(table(pheno.vec)) %>%
    names() %>%
    as.integer()
  majority.frac <- max(table(pheno.vec)) / length(pheno.vec)
  
  if (!is.null(att_to_remove)) {
    # remove attributes (possible confounders) from distance matrix calculation
    tryCatch(
      attr.mat <- attr.mat %>% data.frame() %>%
        select(-att_to_remove),
      error = function(c) "The attribute to remove does not exist."
    )
  }
  
  dist.mat <- attr.mat %>%
    as.matrix() %>%
    unname() %>%
    npdrDistances(metric = nbd.metric, fast.dist = fast.dist) %>%
    as.data.frame() %>%
    `colnames<-`(as.integer(seq.int(num.samp)))
  
  if (nbd.method == "relieff") {
    if (k == 0) { # if no k specified or value 0
      # replace k with the theoretical expected value for SURF (close to multiSURF)
      erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
      # theoretical surf k (sd.frac=.5) for regression problems (does not depend on a hit/miss group)
      # k.alpha <- function(m,alpha){
      #  floor((m - 1) * (1 - erf(alpha / sqrt(2))) / 2)
      # }
      k <- knnSURF(num.samp, sd.frac) # uses sd.frac
      # we will use different k for imbalanced data
    }
    
    if (dopar.nn) {
      avai.cors <- parallel::detectCores() - 2
      cl <- parallel::makeCluster(avai.cors)
      doParallel::registerDoParallel(cl)
      Ri.nearestPairs.list <- foreach::foreach(
        Ri.int = seq.int(num.samp), .packages = c("dplyr")
      ) %dopar% {
        Ri.distances <- dist.mat[Ri.int, ] # all distances to sample Ri
        Ri.nearest <- order(Ri.distances, decreasing = F) # closest to farthest
        # consider distance distributions of hits and misses separately
        Ri.hits <- Ri.nearest[pheno.vec[Ri.int] == pheno.vec[Ri.nearest]]
        Ri.misses <- Ri.nearest[pheno.vec[Ri.int] != pheno.vec[Ri.nearest]]
        
        # fix imbalance 7-28-21
        m.hits <- length(Ri.hits) - 1
        m.miss <- length(Ri.misses)
        
        pheno.tab <- as.numeric(table(as.character(pheno.vec)))
        if (pheno.tab[1] == pheno.tab[2]) {
          k.hits <- floor(0.5 * knnSURF(num.samp - 1, sd.frac))
          k.miss <- k.hits
        } else {
          k.hits <- knnSURF(m.hits, sd.frac)
          k.miss <- knnSURF(m.miss, sd.frac)
        }
        
        Ri.nearest.idx <- c(Ri.hits[2:(k.hits + 1)], Ri.misses[1:k.miss])
        
        # 7-28-21
        # make hit and miss neighborhoods the same size
        # depending on whether Ri is majority or minority class, the number of hits/misses changes
        
        if (!is.null(Ri.nearest.idx)) { # if neighborhood not empty
          # bind automatically repeated Ri, make sure to skip Ri self
          return(data.frame(Ri_idx = Ri.int, NN_idx = Ri.nearest.idx))
        }
      }
      ## [1] 1.000000 1.414214 1.732051
      #parallel::stopCluster(cl)
    } else { # begin non-parallel version
      Ri.nearestPairs.list <- vector("list", num.samp)
      for (Ri in colnames(dist.mat)) { # for each sample Ri
        Ri.int <- as.integer(Ri)
        Ri.distances <- dist.mat[Ri, ] # all distances to sample Ri
        Ri.nearest <- order(Ri.distances, decreasing = F) # closest to farthest
        # consider distance distributions of hits and misses separately
        Ri.hits <- Ri.nearest[pheno.vec[Ri.int] == pheno.vec[Ri.nearest]]
        Ri.misses <- Ri.nearest[pheno.vec[Ri.int] != pheno.vec[Ri.nearest]]
        # 7-28-21 better way for imbalance
        m.hits <- length(Ri.hits) - 1
        m.miss <- length(Ri.misses)
        pheno.tab <- as.numeric(table(as.character(pheno.vec)))
        if (pheno.tab[1] == pheno.tab[2]) {
          k.hits <- floor(0.5 * knnSURF(num.samp - 1, sd.frac))
          k.miss <- k.hits
        } else {
          k.hits <- knnSURF(m.hits, sd.frac)
          k.miss <- knnSURF(m.miss, sd.frac)
        }
        Ri.nearest.idx <- c(Ri.hits[2:(k.hits + 1)], Ri.misses[1:k.miss])
        
        # for misses, option to use farthest is not a good idea because it makes all variables appear
        # different between groups, even null variables
        # if (miss.ordering=="farthest"){ # choose misses that are farthest from Ri
        #  Ri.misses <- rev(Ri.misses)
        #    }
        #
        # make hit and miss neighborhoods the same size (balanced)
        # depending on whether Ri is majority or minority class, the number of hits/misses changes
        
        # 7-28-21 replacing with above method for imbalance
        # hits.frac <- if (pheno.vec[Ri.int] == majority.pheno) majority.frac else (1 - majority.frac)
        # Ri.nearest.idx <- c(
        #  Ri.hits[2:floor(hits.frac * k + 1)], # (2) skip Ri self
        #  Ri.misses[1:floor((1 - hits.frac) * k + 1)]
        # )
        
        if (!is.null(Ri.nearest.idx)) { # if neighborhood not empty
          # bind automatically repeated Ri, make sure to skip Ri self
          Ri.nearestPairs.list[[Ri.int]] <- data.frame(Ri_idx = Ri.int, NN_idx = Ri.nearest.idx)
        }
      } # end for
      # Ri_NN.idxmat <- bind_rows(Ri.nearestPairs.list)
    } # end else dopar.nn
    
    Ri_NN.idxmat <- bind_rows(Ri.nearestPairs.list)
  } else { # surf or multisurf...
    
    # For treating hit/miss distance distributions separately, compute separate hit and miss radii
    # User might want to shrink alpha standard deviation fraction. Unlike relieff, the hit and miss
    # neighborhoods are not balanced.
    if (nbd.method == "surf") { # compute surf radii
      
      hit.dist.rows <- vector("list", num.samp)
      for (i in seq(1, num.samp)) {
        hit.mask <- pheno.vec[i] == pheno.vec
        hit.mask <- hit.mask[-i] # remove self
        hit.dist.rows[[i]] <- dist.mat[i, hit.mask]
      }
      hit.dist.vec <- unlist(hit.dist.rows)
      # average of all hit neighbors
      Ri.hit.radii <- rep(mean(hit.dist.vec) - sd.frac * sd(hit.dist.vec), num.samp)
      names(Ri.hit.radii) <- as.character(1:num.samp)
      
      miss.dist.rows <- vector("list", num.samp)
      for (i in seq(1, num.samp)) {
        miss.mask <- pheno.vec[i] != pheno.vec
        miss.dist.rows[[i]] <- dist.mat[i, miss.mask]
      }
      miss.dist.vec <- unlist(miss.dist.rows)
      # average of all miss neighbors
      Ri.miss.radii <- rep(mean(miss.dist.vec) - sd.frac * sd(miss.dist.vec), num.samp)
      names(Ri.miss.radii) <- as.character(1:num.samp)
    } # end surf radius calc
    
    if (nbd.method == "multisurf") { # compute multisurf radii
      
      Ri.hit.radii <- vector("numeric", num.samp)
      Ri.miss.radii <- vector("numeric", num.samp)
      for (i in seq.int(num.samp)) {
        # grab neighbors that are hits of Ri
        hit.mask <- pheno.vec[i] == pheno.vec
        hit.mask <- hit.mask[-i] # remove self
        hit.dist.row <- as.numeric(dist.mat[i, hit.mask])
        Ri.hit.radii[i] <- mean(hit.dist.row) - sd.frac * sd(hit.dist.row)
        
        # grab neighbors that are misses of Ri
        miss.mask <- pheno.vec[i] != pheno.vec
        miss.dist.row <- as.numeric(dist.mat[i, miss.mask])
        Ri.miss.radii[i] <- mean(miss.dist.row) - sd.frac * sd(miss.dist.row)
      }
      
      names(Ri.hit.radii) <- names(Ri.miss.radii) <- as.character(1:num.samp)
    } # end multisurf radii calc
    
    if (dopar.nn) {
      avai.cors <- parallel::detectCores() - 2
      cl <- parallel::makeCluster(avai.cors)
      doParallel::registerDoParallel(cl)
      Ri.nearestPairs.list <- foreach::foreach(
        Ri.int = seq.int(num.samp), .packages = c("dplyr")
      ) %dopar% {
        Ri.distances <- dist.mat[Ri.int, ]
        Ri.nearest.hits <- which((pheno.vec[Ri.int] == pheno.vec) & (Ri.distances < Ri.hit.radii[Ri.int]) &
                                   (Ri.distances > 0)) # skip Ri self (dist=0)
        Ri.nearest.misses <- which((pheno.vec[Ri.int] != pheno.vec) & (Ri.distances < Ri.miss.radii[Ri.int]))
        # join hit and miss into one nbd
        Ri.nearest.idx <- c(Ri.nearest.hits, Ri.nearest.misses)
        
        if (!is.null(Ri.nearest.idx)) { # similar to relieff
          return(data.frame(Ri_idx = Ri.int, NN_idx = Ri.nearest.idx))
        }
      }
      #parallel::stopCluster(cl)
    } else {
      # put each Ri's nbd in a list then rbind them at the end with bind_rows()
      Ri.nearestPairs.list <- vector("list", num.samp) # initialize list
      
      for (Ri in colnames(dist.mat)) { # for each sample Ri
        Ri.int <- as.integer(Ri)
        Ri.distances <- dist.mat[Ri.int, ]
        Ri.nearest.hits <- which((pheno.vec[Ri.int] == pheno.vec) & (Ri.distances < Ri.hit.radii[Ri.int]) &
                                   (Ri.distances > 0)) # skip Ri self (dist=0)
        Ri.nearest.misses <- which((pheno.vec[Ri.int] != pheno.vec) & (Ri.distances < Ri.miss.radii[Ri.int]))
        # join hit and miss into one nbd
        Ri.nearest.idx <- c(Ri.nearest.hits, Ri.nearest.misses)
        
        if (!is.null(Ri.nearest.idx)) { # similar to relieff
          Ri.nearestPairs.list[[Ri.int]] <- data.frame(Ri_idx = Ri.int, NN_idx = Ri.nearest.idx)
        }
      } # end for
    } # end else dopar.nn
    Ri_NN.idxmat <- bind_rows(Ri.nearestPairs.list)
  }
  
  
  if (neighbor.sampling == "unique") {
    # if you only want to return unique neighbors
    Ri_NN.idxmat <- uniqueNeighbors(Ri_NN.idxmat)
  }
  
  # matrix of Ri's (first column) and their NN's (second column)
  Ri_NN.idxmat
}

# =========================================================================#
#' uniqueNeighbors
#'
#' Find pairs of unique nearest neighbors pairs from possible redundant pairs.
#' Used as options (neighbor.sampling="unique") in nearestNeighbors and npdr functions.
#'
#' @param neighbor.pairs two columns of (possibly redundant) "i,j" pairs from nearestNeighbors function
#' @return two columns of unique pairs from the possibly redundant input, i.e.,
#' new neighborhood pair matrix of only unique neighbor pairs
#' @export
#'
#' @examples
#' neighbor.pairs.idx <- nearestNeighbors(
#'   predictors.mat,
#'   nbd.method = "multisurf",
#'   nbd.metric = "manhattan",
#'   sd.frac = 0.5
#' )
#' head(uniqueNeighbors(neighbor.pairs.idx))
uniqueNeighbors <- function(neighbor.pairs) {
  # return:
  # sort and make create redundant vector of "i,j" pairs
  # e.g., pairs 1  36 and 36  1 both become 1  36
  pair_str <- data.frame(
    xmin = pmin(neighbor.pairs[, 1], neighbor.pairs[, 2]),
    xmax = pmax(neighbor.pairs[, 1], neighbor.pairs[, 2])
  ) %>%
    transmute(combined = paste0(xmin, ",", xmax))
  
  neighbor.pairs[!duplicated(pair_str), ]
}

# =========================================================================#
#' knnVec
#'
#' Number of neighbors for each sample (vector) from a neighbor-pair matrix.
#'
#' @param neighbor.pairs.mat two columns of redundant "i,j" pairs from nearestNeighbors function
#' @return  knn.vec vector number of nearest neighbors for each instance
#'
#' @examples
#' neighbor.pairs.idx <- nearestNeighbors(
#'   predictors.mat,
#'   nbd.method = "multisurf",
#'   nbd.metric = "manhattan",
#'   sd.frac = 0.5
#' )
#' mean(knnVec(neighbor.pairs.idx)) # average number of neighbors
#' @export
#'
knnVec <- function(neighbor.pairs.mat) {
  data.frame(neighbor.pairs.mat) %>%
    count(Ri_idx) %>%
    pull(n)
}


#' Get the indices of samples nearest to Ri
#'
#' @param dist.mat Distance matrix
#' @param Ri Sample Ri index.
#' @param nbd.method neighborhood method \code{"multisurf"} or \code{"surf"} (no k) 
#' or \code{"relieff"} (specify k). Used by nearestNeighbors().
#' @param k Integer for the number of neighbors (\code{"relieff"} method).
#' @param Ri.radius Radius of the neighborhood (other methods).
#'
#' @return Numeric vector of nearest indices.
#' @export
#' 
get_Ri_nearest <- function(dist.mat, Ri, nbd.method, k = 0, Ri.radius = NULL) {
  dist.mat %>%
    select(!!Ri) %>% # select the column Ri, hopefully reduce processing power
    rownames2columns() %>% # push the neighbors from rownames to a column named rowname
    {if (nbd.method == "relieff"){
      slice_min(., order_by = !!sym(Ri), n = k + 1) %>% # select the k closest neighbors, include self
        filter(., (!!sym(Ri)) > 0)
      #arrange(!!sym(Ri)) # sort by increasing distance
    } else {
      filter(., ((!!sym(Ri)) < Ri.radius[Ri]) & ((!!sym(Ri)) > 0))
      #arrange(!!sym(Ri)) # sort by increasing distance
    }} %>% # top_n does not sort output, so make sure remove self
    pull(rowname) %>% # get the neighbors
    as.integer() # convert from string (rownames - not factors) to integers
}


#' Get the indices of samples nearest to Ri
#'
#' @param dist.mat Distance matrix
#' @param Ri Sample Ri index.
#' @param nbd.method neighborhood method \code{"multisurf"} or \code{"surf"} (no k) 
#' or \code{"relieff"} (specify k). Used by nearestNeighbors().
#' @param k Integer for the number of neighbors (\code{"relieff"} method).
#' @param Ri.radius Radius of the neighborhood (other methods).
#'
#' @return Numeric vector of nearest indices.
#' @export
#' 
get_Ri_nearest2 <- function(dist.mat, Ri, nbd.method, k = 0, Ri.radius = NULL) {
  
  df <- dist.mat %>%
    select(!!Ri) %>% 
    tibble::rownames_to_column(var = "rowname")
  
  if (nbd.method == "relieff") {
    df <- df %>%
      slice_min(order_by = .data[[Ri]], n = k + 1) %>%
      filter(.data[[Ri]] > 0) %>%
      arrange(.data[[Ri]]) # sort by increasing distance
  } else {
    df <- df %>%
      filter((.data[[Ri]] < Ri.radius[Ri]) & (.data[[Ri]] > 0)) %>%
      arrange(.data[[Ri]]) # sort by increasing distance
  }
  
  neighbors <- df %>% pull(rowname) %>% as.integer()
  
  return(neighbors)
}


