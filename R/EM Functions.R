find.outermost.vertices = function(x) {
  # Accepts a matrix with vertices as columns or a list of vectors.
  # Returns indices ordered by decreasing L2 distance from the centre.
  if (is.list(x)) x = do.call(cbind, x)
  centre = rowMeans(x)
  distances = colSums((x - centre)^2)
  return(order(distances, decreasing = TRUE))
}

find.convex.combination = function(ext, relint, min.wt = 0) {
  # Both parameters can be a matrix with vertices as columns or a list of 
  # vectors. Returns a matrix of coefficients that minimize the L2 norm, where 
  # each column is the coefficients corresponding to columns of relint.
  if (is.list(ext)) ext = do.call(cbind, ext)
  if (is.list(relint)) ext = do.call(cbind, relint)
  n.ext = ncol(ext)
  d = nrow(ext)
  quadratic = crossprod(ext)
  if (d < n.ext) {
    # If rank-deficient, we will also minimize x^T x.
    diag(quadratic) = diag(quadratic) + 1
  }
  
  coefs = apply(relint, 2, function(x) {
    solve.QP(quadratic,
             colSums(x * ext),
             cbind(1, diag(n.ext)),
             c(1, rep(min.wt, n.ext)),
             1)$solution
  })
  return(coefs)
}

weighted.sum = function(list, weights) {
  return(Reduce(`+`, mapply("*", list, weights, SIMPLIFY = FALSE)))
}

compute.z.weights = function(data, etas, lambdas, Pi, S_mat = apply(data, 1, tcrossprod)) {
  if (ncol(etas) != dim(lambdas)[3]) stop("Incorrect inputs.")
  K = ncol(etas)
  
  row_addition = log(as.vector(Pi)) + (-0.9189385332046727) * ncol(data)
  for (i in 1:K) {
    logdetlambda = log(det(lambdas[,,i]))
    h1 = etas[,i] %*% solve(lambdas[,,i], etas[,i])
    row_addition[i] = row_addition[i] + 0.5 * (logdetlambda - h1)
  }
  
  Lmat = sapply(1:K, function(k) lambdas[,,k])
  new_z_nk = data %*% etas - (crossprod(S_mat, Lmat) / 2)
  new_z_nk = t(t(new_z_nk) + row_addition)
  new_z_nk = new_z_nk - rowMaxs(new_z_nk, value = TRUE)
  new_z_nk = prop.table(exp(new_z_nk), 1)
  return(new_z_nk)
}

number.of.parameters = function(state) {
  Kp = ncol(state$alphas)
  Kc = nrow(state$alphas)
  d = ncol(state$data)
  
  return(Kp * (d + d*(d+1)/2) + Kc * min(Kp - 1, d + d*(d+1)/2) + (Kp + Kc - 1))
}

get.bic = function(state) {
  return(2 * tail(na.omit(as.numeric(state$ll_trace)), 1) - 
           number.of.parameters(state) * log(nrow(state$data)))
}

initialisation.kmeans = function(data, Kp, Kc) {
  K = Kp + Kc
  d = ncol(data)
  S_mat = apply(data, 1, tcrossprod)
  
  kmeans.fit = kmeans(data, K)
  mus = t(kmeans.fit$centers)
  distances = find.outermost.vertices(mus)
  
  alphas = t(find.convex.combination(mus[,distances[1:Kp], drop = FALSE],
                                     mus[,-distances[1:Kp], drop = FALSE]))
  
  cov.estimate = lapply(1:K, function(k) {
    cov.wt(data[kmeans.fit$cluster == k, ], method = "unbiased")$cov
  })
  Pi = prop.table(kmeans.fit$size)
  
  # Assumption of equal covariance throughout
  # Could be bumped to spherical
  wt.cov.equal = weighted.sum(cov.estimate, Pi)
  lambda.equal = spdinv(wt.cov.equal)
  
  etas = lambda.equal %*% mus[,distances[1:Kp]]
  etas = cbind(etas, tcrossprod(etas, alphas))
  lambdas = array(lambda.equal, dim = c(d, d, K))
  z_nk = compute.z.weights(data, etas, lambdas, Pi, S_mat)
  sigmas = array(wt.cov.equal, dim = c(d, d, K))
  Pi = colMeans(z_nk)
  
  state = list(gamma = z_nk, Pi = Pi, etas = etas, lambdas = lambdas,
               sigmas = sigmas, alphas = alphas, data = data, S_mat = S_mat,
               ll_trace = numeric(0), gamma_trace = numeric(0),
               delta_V_trace = numeric(0))
  return(state)
}

initialisation.mclust = function(data, Kp, Kc, modelname = "VVV") {
  K = Kp + Kc
  d = ncol(data)
  S_mat = apply(data, 1, tcrossprod)
  
  mclust.fit = Mclust(data, K, modelname, verbose = FALSE)
  if (is.null(mclust.fit)) stop("Failed to initialise with mclust!")
  mclust.mus = mclust.fit$parameters$mean
  mclust.sigmas = asplit(mclust.fit$parameters$variance$sigma, 3)
  
  etas = sapply(1:K, function(k) {solve(mclust.sigmas[[k]], mclust.mus[,k])})
  
  distances = find.outermost.vertices(etas)
  etas.P = etas[,distances[1:Kp], drop = FALSE]
  etas.C = etas[,-distances[1:Kp], drop = FALSE]
  alphas = t(find.convex.combination(etas.P, etas.C))
  
  mclust.sigmas = mclust.sigmas[distances]
  
  etas = cbind(etas.P, tcrossprod(etas.P, alphas))
  lambdas = sigmas = array(NA, dim = c(d, d, K))
  for (p in 1:Kp) {
    lambdas[,,p] = spdinv(mclust.sigmas[[p]])
    sigmas[,,p] = mclust.sigmas[[p]]
  }
  for (c in 1:Kc) {
    lambdas[,,Kp+c] = weighted.sum(asplit(lambdas[,,1:Kp], 3), alphas[c,])
    sigmas[,,Kp+c] = spdinv(lambdas[,,Kp+c])
  }
  Pi = mclust.fit$parameters$pro[distances]
  z_nk = compute.z.weights(data, etas, lambdas, Pi, S_mat)
  Pi = colMeans(z_nk)
  
  state = list(gamma = z_nk, Pi = Pi, etas = etas, lambdas = lambdas,
               sigmas = sigmas, alphas = alphas, data = data, S_mat = S_mat,
               ll_trace = numeric(0), gamma_trace = numeric(0),
               delta_V_trace = numeric(0))
  return(state)
}

visualise.state = function(state, pairs.cols = c(1,2)) {
  Kp = ncol(state$alphas)
  Kc = nrow(state$alphas)
  cbbPalette = palette.colors(palette = "classic tableau")
  ll_trace = state$ll_trace
  i = length(ll_trace)
  assign = rowMaxs(state$gamma)
  print(Table(assign))
  old.par = par(mfrow = c(1, 3))
  lowerbound = quantile(ll_trace[1:i][!is.infinite(ll_trace[1:i])], 0.05, na.rm = TRUE)
  plot(c(1,length(ll_trace)), c(lowerbound, max(ll_trace[1:i])), col = "transparent", xlab = "Iteration", ylab = "Log-Likelihood", 
       main = paste0("LL: ", round(ll_trace[i], 2)))
  lines(ll_trace[1:i], lwd = 2)
  plot(state$data[,pairs.cols], col = cbbPalette[assign %% length(cbbPalette) + 1], 
       pch = c(20,21)[(assign > Kp) + 1], cex = 2,
       main = paste0("BIC: ", round(get.bic(state), 2)))
  plot(c(1,length(ll_trace) - 1), c(1e-13, 1), col = "transparent", log = 'y', xlab = "Iteration", ylab = "Delta Log-Likelihood", main = i)
  lines(Pmax(diff(ll_trace[1:i]), rep(1e-13, i)), lwd = 2)
  par(old.par)
}

diagnose.convergence = function(state) {
  cbbPalette = palette.colors(palette = "classic tableau")
  ll_trace = state$ll_trace
  i = length(ll_trace)
  old.par = par(mfrow = c(1, 3))
  plot(c(1,length(ll_trace) - 1), c(1e-13, 1), col = "transparent", log = 'y', xlab = "Iteration", ylab = "Delta Log-Likelihood", 
       main = paste0("Last Loglik Change: ", signif(diff(tail(ll_trace, 2)),2)))
  lines(Pmax(diff(ll_trace), rep(1e-13, i)), lwd = 2)
  plot(c(1,length(ll_trace) - 1), c(1e-13, 1), col = "transparent", log = 'y', xlab = "Iteration", ylab = "|Delta Z|_1", 
       main = paste0("Last sum(abs(Delta Z)): ", signif(tail(state$gamma_trace, 1), 2)))
  lines(state$gamma_trace, lwd = 2)
  plot(c(1,length(ll_trace) - 1), c(1e-13, 1), col = "transparent", log = 'y', xlab = "Iteration", ylab = "|Delta theta|_1 / # Parameters", 
       main = paste0("Last mean(abs(Delta theta)): ", signif(tail(state$delta_V_trace, 1), 2)))
  lines(state$delta_V_trace, lwd = 2)
  par(old.par)
}

multiple.initialisations = function(data, Kp, Kc, niters = 1000,
                                    modelnames = c("EII", "VII", "EEI", "VEI", 
                                                   "EVI", "VVI", "EEE", "EEV", 
                                                   "VEV", "VVV", "EVE", "VVE", 
                                                   "VEE", "EVV", "kmeans"),
                                    nr.maxit = 100L,
                                    nr.eps = 1e-12,
                                    hold.z.iters = 0L,
                                    verbosity = 0L,
                                    return.all.models = FALSE,
                                    min.eigenvalue = 1e-4,
                                    min.zg = 1e-8) {
  fitted.models = lapply(modelnames, function(modelname) {
    if (verbosity > 0) cat("Starting fit with ", modelname, "...\n")
    state = tryCatch({
      if (modelname == "kmeans") 
        state = initialisation.kmeans(data, Kp, Kc)
      else 
        state = initialisation.mclust(data, Kp, Kc, modelname)
      state = runEM(niters, state, show_progress = verbosity > 1, 
                    holdziters = hold.z.iters, nr_maxit = nr.maxit, nr_eps = nr.eps)
      
      sigma.eigenvalues = sapply(apply(state$sigmas, 3, eigen, only.values = TRUE), `[[`, "values")
      if (# Handle corner cases with failures.
        is.null(state) ||
        !is.finite(tail(state$ll_trace, 1)) ||     # NaN parameters
        any(sigma.eigenvalues < min.eigenvalue) || # Degenerate covariances
        any(colsums(state$gamma) < min.zg)         # Empty assignments
      ) state = NULL 
      state
    }, error = function(e) {
      if (verbosity > 0) print(e)
      NULL
    })
    if (verbosity > 0) {
      if (!is.null(state)) cat("=> BIC: ", round(get.bic(state), 2), "\n") 
      else cat("=> Bad solution.\n")
    }
    return(state)
  })
  names(fitted.models) = modelnames
  
  bics = sapply(fitted.models, function(x) if (!is.null(x)) get.bic(x) else NA)
  if (all(is.na(bics))) {
    warning("No initialisations produced feasible solutions.")
    return("NULL")
  }
  if (verbosity > 0) cat("==> Best Model: ", modelnames[which.max(bics)], "\n")
  if (return.all.models) return(list(
    names = modelnames,
    models = fitted.models,
    bics = bics,
    best.index = which.min(bics)
  )) else return(fitted.models[[which.max(bics)]])
}

continue.em = function(state,
                       niters = NULL,
                       nr.maxit = 100L,
                       nr.eps = 1e-12,
                       hold.z.iters = 0L,
                       verbosity = 0L) {
  state = runEM(niters, state, show_progress = verbosity > 0, 
                holdziters = hold.z.iters, nr_maxit = nr.maxit, nr_eps = nr.eps)
  return(state)
}

get.ari = function(state, true.classes) {
  adjustedRandIndex(rowMaxs(state$gamma), true.classes)
}

chimeralclustering = function(data, Kp, Kc, 
                              niters = 5000L,
                              hold.z.iters = 0L,
                              mini.em.iters = 1000L,
                              modelnames = c("EII", "VII", "EEI", "VEI", 
                                             "EVI", "VVI", "EEE", "EEV", 
                                             "VEV", "VVV", "EVE", "VVE", 
                                             "VEE", "EVV", "kmeans"),
                              nr.maxit = 100L,
                              nr.eps = 1e-12,
                              verbosity = 0L,
                              min.eigenvalue = 1e-4,
                              min.zg = 1e-8) {
  
}