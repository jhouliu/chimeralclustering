generate.k.wreath = function(k,
                             n,
                             r = 64 / (2 * tan(pi / k)),
                             twist = pi / 3,
                             wts = c(1 / k, (1 - 1 / k) / 2)
) {
  rotation.matrix = function(theta) 
    matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)), 2, 2, byrow = TRUE)
    
  Kp = k
  angles = seq(0, 2*pi, length.out = Kp + 1)[-1]
  mu.list = lapply(angles, function(angle) as.numeric(rotation.matrix(-angle) %*% c(0, r)))

  cov.list = lapply(1:Kp, function(p) {
    rot = rotation.matrix(pi/2 - angles[p] + twist)
    sigma = rot %*% diag(c(32, 2)) %*% t(rot)
    return(sigma)
  })
  
  eta.P = lapply(1:Kp, function(p) {
    as.numeric(rotation.matrix(pi/2 - angles[p] + twist) %*% diag(1/c(32, 2)) %*% 
      t(rotation.matrix(pi/2 + twist)) %*% c(0, r))
  })
  lambda.P = lapply(cov.list, solve)
  
  alphas = diag(wts[2], Kp)
  alphas[alphas == 0] = wts[1]
  alphas[row(alphas) == col(alphas) + 1] = alphas[1, Kp] = wts[2]

  # alphas = diag(1, Kp)
  # alphas[alphas == 0] = 1
  # alphas[col(alphas) - 1 == (row(alphas) - 1) %% Kp] = 1/Kp
  # alphas[col(alphas) - 1 == (row(alphas) + 0) %% Kp] = 1/Kp
  # alphas[col(alphas) - 1 == (row(alphas) + 1) %% Kp] = 1/Kp
  
  # alphas = matrix(0, Kp, Kp)
  # if (Kp %% 2 == 0) {
  #   # Even, pick two opposites
  #   alphas = diag(0.5, Kp)
  #   # alphas[col(alphas) - 1 == (row(alphas) + Kp/2 - 0) %% Kp] = (1 - 1/Kp^0.5) / 2
  #   alphas[col(alphas) - 1 == (row(alphas) + Kp/2 - 1) %% Kp] = 0.5
  #   # alphas[col(alphas) - 1 == (row(alphas) + Kp/2 - 2) %% Kp] = (1 - 1/Kp^0.5) / 2
  # } else {
  #   # Odd, pick one and two opposites
  #   # alphas = diag(1/Kp, Kp)
  #   # alphas[col(alphas) - 1 == (row(alphas) + floor(Kp/2) - 1) %% Kp] = (1 - 1/Kp) / 2
  #   # alphas[col(alphas) - 1 == (row(alphas) + ceiling(Kp/2) - 1) %% Kp] = (1 - 1/Kp) / 2
  #   alphas = diag(0.2, Kp)
  #   alphas[col(alphas) - 1 == (row(alphas) + floor(Kp/2) - 1) %% Kp] = 0.4
  #   alphas[col(alphas) - 1 == (row(alphas) + ceiling(Kp/2) - 1) %% Kp] = 0.4
  # }

  alphas = prop.table(alphas, 1)
  
  eta.C = lapply(asplit(alphas, 1), weighted.sum, list = eta.P)
  lambda.C = lapply(asplit(alphas, 1), weighted.sum, list = lambda.P)
  
  etas = c(eta.P, eta.C)
  lambdas = c(lambda.P, lambda.C)
  
  data = lapply(1:(Kp * 2), function(i) {
    mvtnorm::rmvnorm(n, solve(lambdas[[i]], etas[[i]]), solve(lambdas[[i]]))
  })
  data = do.call(rbind, data)
  true.classes = rep(1:(Kp * 2), rep(n, Kp * 2))
  
  gamma = matrix(0, 2 * n * Kp, 2 * Kp)
  gamma[cbind(1:(2*n*Kp), true.classes)] = 1
  
  return(list(
    data = data,
    S_mat = apply(data, 1, tcrossprod),
    true.classes = true.classes,
    gamma = gamma,
    k = k,
    etas = do.call(cbind, etas),
    lambdas = `dim<-`(do.call(cbind, lambdas), c(2, 2, 2 * Kp)),
    sigmas = `dim<-`(do.call(cbind, lapply(lambdas, spdinv)), c(2, 2, 2 * Kp)),
    alphas = alphas,
    Pi = rep(1 / (2 * Kp), 2 * Kp),
    ll_trace = numeric(0),
    gamma_trace = numeric(0),
    delta_V_trace = numeric(0)
  ))
}

generate.d.radioactive = function(d, n) {
  nullspace = function(x) {
    qr = qr(x)
    indices = if (qr$rank == 0) 1:ncol(x) else -(1:qr$rank)
    return(qr.Q(qr, complete = TRUE)[, indices, drop = FALSE])
  }
  
  # Generate mu's by creating a regular simplex with d + 1 vertices.
  
  # Function translated from MATLAB
  # Source: http://wiki.gis.com/wiki/index.php/Simplex
  v = matrix(0, d + 1, d)
  v[1,1] = 1
  v[2:(d+1), 1] = -1 / d
  if (d > 1) {
    for (i in 2:d) {
      c = sqrt(1 - sum(v[i, 1:(i-1)]^2))
      v[i, i] = c
      w = v[i, 1:(i-1)]
      a = -1/c * (1/d + sum(w^2))
      v[(i+1):(d+1), i] = a
    }
  }
  
  mu.list = asplit(sqrt(200 * d) * v, 1)
  cov.list = lapply(mu.list, function(mu) {
    if (all(mu == 0)) return(diag(d))
    mu = mu / sqrt(sum(mu^2))
    # diag(1000, d) - 999 * tcrossprod(mu)
    rot = cbind(mu, nullspace(mu))
    rot %*% diag(c(1, rep(100, d - 1))) %*% t(rot)
  })
  
  prototype.data = lapply(1:(d+1), function(i) 
    rmvnorm(n, mu.list[[i]], cov.list[[i]]))
  
  eta.P = mapply(solve, cov.list, mu.list, SIMPLIFY = FALSE)
  lambda.P = lapply(cov.list, solve)
  
  alphas = diag(3, d + 1)
  alphas[alphas == 0] = 1
  alphas = rbind(alphas, 1)
  alphas = prop.table(alphas, 1)
  
  eta.C = lapply(asplit(alphas, 1), weighted.sum, list = eta.P)
  lambda.C = lapply(asplit(alphas, 1), weighted.sum, list = lambda.P)
  
  etas = c(eta.P, eta.C)
  lambdas = c(lambda.P, lambda.C)
  
  data = lapply(1:(2 * d + 3), function(i) {
    mvtnorm::rmvnorm(n, solve(lambdas[[i]], etas[[i]]), solve(lambdas[[i]]))
  })
  data = do.call(rbind, data)
  true.classes = rep(1:(2 * d + 3), rep(n, 2 * d + 3))
  
  gamma = matrix(0, n * (2 * d + 3), 2 * d + 3)
  gamma[cbind(1:(n * (2 * d + 3)), true.classes)] = 1
  
  return(list(
    data = data,
    S_mat = apply(data, 1, tcrossprod),
    true.classes = true.classes,
    gamma = gamma,
    d = d,
    etas = do.call(cbind, etas),
    lambdas = `dim<-`(do.call(cbind, lambdas), c(d, d, 2 * d + 3)),
    sigmas = `dim<-`(do.call(cbind, lapply(lambdas, spdinv)), c(d, d, 2 * d + 3)),
    alphas = alphas,
    Pi = rep(1 / (2 * d + 3), (2 * d + 3)),
    ll_trace = numeric(0),
    gamma_trace = numeric(0),
    delta_V_trace = numeric(0)
  ))
}
