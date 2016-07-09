# Code for estimating Markov random field parameters for k = 5 classes, i.e. 3 species, 
# using the simulated field algorithm described in Celeux 2003.
# The only input that is required is a list of priors, where each element of 
# the list is a matrix of priors, with classes in columns and genes in rows. 
# A toy example can be found in the file ex.priors.rda. 

library(LaplacesDemon)
##########
# to install Lapalces Demon:
# install.packages("devtools")
# library(devtools)
# install_github("ecbrown/LaplacesDemon")

run.example <- function() {
  load("ex.priors.rda")
  # solve for MRF parameters
  res <- solveMRF(ex.priors)
  
  # get posterior probabilities of class membership
  ex.post <- posterior(res$a, res$b, ex.priors)
  
  # get the posterior best classifications
  postb <- best.mod.list(ex.post)
  
  return(list(res = res, pp = ex.priors, postp = ex.post, postb = postb))
}

####################################################
# solves for the MRF parameters for case k = 5
# pp: prior probabilities
# maxSteps: the maximum number of parameter estimation steps
# fix.a: whether or not the alphas should be estimated, or fixed at 0
solveMRF <- function(pp, maxSteps = 20, fix.a = TRUE) {
  a = rep(0, 4)
  b = 0
  
  if (fix.a == TRUE) print(paste("Fixing a = ", paste(a, collapse = " "), sep = ""))
  
  # gets initial best states
  Z.init <- best.mod.list(pp)
  
  # gets initial neighbors
  U.init <- countNeighbors(Z.init)
  
  Z <- Z.init
  U <- U.init
  
  # newton raphson steps; average results over last 5
  last.a <- matrix(NA, ncol = 4, nrow = 5)
  last.b <- rep(NA, 5)
  
  for (i in 1:(maxSteps + 5)) {
    coeffs <- list()
    for (r in names(U)) {
      # get MRF probabilities
      mrfps <- t(sapply(X = 1:nrow(U[[r]]), function(i) mrfProb(a, b, U[[r]][i, ])))
      
      # get full model probabilities
      coeffs[[r]] <- fullProb(mrfps, pp[[r]])
    }
    
    # solve for updated parameters a, b
    sol <- evalFullGradHess(a, b, U, coeffs)
    vars <- c(a, b) - solve(sol$hess) %*% sol$grad
    
    if (fix.a != TRUE) a <- vars[-5]
    b <- vars[5]
    
    last.a[(i-1) %% 5 + 1, ] <- a
    last.b[(i-1) %% 5 + 1] <- b
    if (i %% 5 == 1) print(c(a, b))
    
    # update Z with new a, b
    Z <- updateZ(a, b, Z, pp)
    U <- countNeighbors(Z)
  }
  
  list(a = colMeans(last.a), b = mean(last.b))
}

# get MAP estimates given estimated parameters and prior probabilities
posterior <- function(a, b, pp, burnin = 50, sample = 100) {
  Z <- best.mod.list(pp)
  if (!is.null(colnames(Z))) {
    Regions <- colnames(Z)
  } else {
    Regions <- paste("Region", 1:ncol(Z), sep = "_")
    colnames(Z) <- Regions
  }
  
  print(paste("Burn-in: ", burnin, sep = ""))
  for (i in 1:burnin) {
    Z <- updateZ(a, b, Z, pp)
  }
  
  print(paste("Sample: ", sample, sep = ""))
  
  # posterior samples is taken as a list
  # list has Nregions items
  # each item is a matrix with 5 columns (for k = 5)
  # and Ngenes rows
  postSamp <- list()
  for (r in Regions) {
    postSamp[[r]] <- matrix(0, nrow = nrow(Z), ncol = 5)
  }
  
  for (i in 1:sample) {
    Z <- updateZ(a, b, Z, pp)
    pc <- countPost(Z)
    
    for (r in Regions) {
      postSamp[[r]] <- postSamp[[r]] + pc[[r]]
    }
  }
  
  postp <- lapply(X = postSamp, function(x) x / sample)
  for (r in names(pp)) {
    rownames(postp[[r]]) <- rownames(pp[[r]])
    colnames(postp[[r]]) <- colnames(pp[[r]])
  }
  postp
}

# counts posterior samples and organizes into different format
countPost <- function(Z) {
  L <- list()
  for (r in 1:ncol(Z)) {
    L[[r]] <- t(sapply(X = Z[, r], function(i) {x <- rep(0, 5); x[i] <- 1; return(x)}))
  }
  names(L) <- colnames(Z)
  L
}

# the neighbors matrix is a list of Nregions matrices,
# neighbors of Z of class k = (1, 2, 4, 5) - 3
countNeighbors <- function(Z) {
  U <- list()
  
  for (r in 1:ncol(Z)) {
    nbs <- t(sapply(X = 1:nrow(Z), function(i) countNeighborsbyRow(Z[i, ], r)))
    colnames(nbs) <- paste("u_k", c(1, 2, 4, 5), sep = "")
    U[[r]] <- nbs
  }
  
  names(U) <- colnames(Z)
  U
}

# get neighbor counts for a single row/gene
countNeighborsbyRow <- function(z, r) {
  sapply(X = c(1, 2, 4, 5), function(k) length(which(z[-r] == k)) - length(which(z[-r] == 3)))
}

# given parameters a, b, initial probabilities pp, and 
# state matrix Z, updates Z sequentially via the simulated field algorithm
updateZ <- function(a, b, Z, pp) {
  newZ <- sapply(X = 1:nrow(Z), 
                 function(i) {
                   ppr <- lapply(pp, function(x) x[i, ])  
                   updateZbyRow(a, b, Z[i, ], ppr)
                 })
  
  t(newZ)
}

# given parameters a, b, initial probabilities for a single gene/row, and 
# state vector z for the row, updates z sequentially via the simulated field algorithm
updateZbyRow <- function(a, b, z, ppr) {
  newz <- z
  for (j in 1:length(newz)) {
    u <- sapply(X = c(1, 2, 4, 5), 
                FUN = function(k) length(which(newz[-j] == k)) - length(which(newz[-j] == 3)))
    
    mrfp <- mrfProb(a, b, u)        # update MRF prob
    newz[j] <- classSamp(mrfp, ppr[[j]])   # sample new Z
  }
  newz
}

# sample new class
classSamp <- function(mrfp, fp) {
  p <- fp * mrfp / (sum(fp * mrfp))
  rcat(n = 1, p = p)
}

# get MRF probability for each class k, given a, b, u, and ws
# here u is a vector of length 4, describing the 
# neighbors, for k = {1, 2, 4, 5} - 3
mrfProb <- function(a, b, u) {
  F <- sapply(X = 1:4, function(i) {exp(a[i] + b * u[i])})
  F <- c(F[1:2], 1, F[3:4])
  F / sum(F)
}

# gets full probability given MRF probs mrfp0 and 
# density probs fp0. Both have 5 columns for k = 5
fullProb <- function(mrfp0, fp0) {
  fp <- mrfp0 * fp0
  sweep(fp, MARGIN = 1, FUN = "/", rowSums(fp))
}

# evalutes gradient and hessian over all genes and regions
# i.e., sums
# U: list of neighbors matrices
# coeffs: list of full probability matrices
evalFullGradHess <- function(a, b, U, coeffs) {
  if (length(which(!names(U) %in% names(coeffs))) > 0) 
    return("Region names not the same!")
  
  grad <- rep(0, 5)
  hess <- matrix(0, 5, 5)
  
  for (r in names(U)) {
    g <- t(sapply(1:nrow(U[[r]]), function(i) gradient(a, b, U[[r]][i, ], as.numeric(coeffs[[r]][i, -3]))))
    grad <- grad + colSums(g)
    
    h <- lapply(X = 1:nrow(U[[r]]), function(i) hessian(a, b, U[[r]][i, ]))
    hess <- hess + Reduce("+", h)
  }
  
  list(grad = grad, hess = hess)
}

gradient <- function(a, b, u, coeff) {
  x = exp(a + b * u)
  ws = 1 + sum(x)
  
  # alpha terms
  g <- coeff - x / ws
  
  # beta term
  g <- c(g, sum(coeff * u) - sum(x * u)/ws)
  g  
}

hessian <- function(a, b, u) {
  x = exp(a + b * u)
  ws = 1 + sum(x)
  
  # the hessian matrix is symmetric, so we only have to 
  # calculate the top/bottom triangular matrix
  h <- matrix(NA, 5, 5)
  
  # diagonal entries for alpha
  diag(h)[1:4] <- sapply(X = 1:4, FUN = function(i) -x[i] * (1 + sum(x[-i])) / (ws ^ 2))
  
  # diagonal entry for beta
  h[5, 5] <- (sum(x * u) ^ 2 - sum(x * u ^ 2) * ws) / ws ^ 2
  
  # off diagonal entries for alpha
  p <- combn(1:4, 2)
  for (i in 1:ncol(p)) {
    h[p[1, i], p[2, i]] <- x[p[1, i]] * x[p[2, i]] / (ws ^ 2)
  }
  
  # off diagonal entries for beta
  for (i in 1:4) {
    h[i, 5] <- -x[i] * (u[i] + sum(x[-i] * (u[i] - u[-i]))) / (ws ^ 2)
  }
  
  # make symmetric
  h[lower.tri(h)] <- t(h)[lower.tri(h)]
  h
}

best.mod.list <- function(plist) {
  b <- lapply(plist, best.mod)
  do.call(cbind, b)
}

# gets the best model based on relative p-value
# in the case of ties, randomly select
best.mod <- function(pmat) {
  b <- apply(X = pmat, MARGIN = 1, FUN = function(x) bestr(x)) 
  unlist(b)
}

bestr <- function(r) {
  i <- which(r == max(r))
  if (length(i) > 1) i <- sample(i, 1)
  i
}