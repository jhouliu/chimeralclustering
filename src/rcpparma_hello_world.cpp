// -*- mode: C++; c-indent-level: 2; c-basic-offset: 2; indent-tabs-mode: nil; -*-

#define ARMA_DONT_USE_OPENMP
#define ARMA_DONT_PRINT_ERRORS
#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
#include "riccatiCareSolution.h"

// [[Rcpp::depends(RcppProgress)]]
#include "progress.hpp"

// // [[Rcpp::export]]
// arma::mat care_inv(const arma::mat& aaa,
//                    const arma::mat& qqq,
//                    const arma::mat& ggg) {
//   const unsigned int n = aaa.n_rows;
//   arma::mat x = arma::join_rows(arma::join_cols(aaa, -qqq), arma::join_cols(-ggg, -aaa.t()));
//   
//   arma::mat AA, BB, Q, Z;
//   arma::qz(AA, BB, Q, Z, arma::eye(size(x)), x, "lhp");
//   arma::mat result = arma::solve((Z.submat(n, 0, 2*n - 1, n - 1)).t(),
//                                  (Z.submat(0, 0, n - 1, n - 1)).t(),
//                                  solve_opts::refine + solve_opts::equilibrate);
//   return result;
// }

//[[Rcpp::export]]
void updateE (const arma::mat& etas, const arma::cube& lambdas, const arma::cube& sigmas,
              const arma::mat& data, const arma::mat& S_mat, 
              arma::rowvec& pi, arma::mat &z_nk, 
              const bool& hold_z, double &loglik) {
  arma::mat Lmat = (arma::mat(lambdas.memptr(), lambdas.n_rows*lambdas.n_cols, lambdas.n_slices));
  arma::rowvec row_addition = log(pi) + -0.9189385332046727 * data.n_cols;
  for (unsigned int i = 0; i < pi.size(); i++) {
    double logdetlambda = log(det(lambdas.slice(i)));
    if (std::isnan(logdetlambda)) {row_addition(i) = -INFINITY;
    } else {row_addition(i) += 0.5 * (logdetlambda - arma::as_scalar(etas.col(i).t() * sigmas.slice(i) * etas.col(i)));}
  }
  arma::mat new_z_nk = (data * etas - (S_mat.t() * Lmat) / 2);
  new_z_nk.each_row() += row_addition;
  loglik = arma::sum(arma::log(arma::sum(exp(new_z_nk), 1)));
  if (!hold_z) {// [[Rcpp::export]]
    new_z_nk.each_col() -= arma::max(new_z_nk, 1);
    z_nk = arma::normalise(exp(new_z_nk), 1, 1);
    pi = arma::sum(z_nk) / z_nk.n_rows;
  }
  return;
}

//[[Rcpp::export]]
void updateAllEtaP (const arma::mat& data, const arma::cube& lambdas, const arma::cube& sigmas,
                    arma::mat& etas, const arma::mat& alphas, const arma::mat& z_nk) {
  unsigned int Kp = alphas.n_cols;
  unsigned int Kc = alphas.n_rows;
  unsigned int d = data.n_cols;
  
  arma::mat ygs = data.t() * z_nk;
  arma::rowvec p_g = arma::sum(z_nk, 0);
  arma::mat lhs(d * Kp, d * Kp, fill::zeros);
  arma::vec rhs(d * Kp);
  
  for (unsigned int g = 0; g < Kp; g++) {
    rhs.subvec(g * d, (g + 1) * d - 1) = ygs.cols(Kp, Kp + Kc - 1) * alphas.col(g) + ygs.col(g);
    for (unsigned int p = g; p < Kp; p++) {
      if (g == p) {lhs.submat(p * d, g * d, (p + 1) * d - 1, (g + 1) * d - 1) += sigmas.slice(p) * p_g(p);}
      for (unsigned int i = 0; i < Kc; i++) 
        lhs.submat(p * d, g * d, (p + 1) * d - 1, (g + 1) * d - 1) += alphas(i, g) * alphas(i, p) * p_g(Kp + i) * sigmas.slice(Kp + i);
    }
  }
  
  lhs = arma::symmatl(lhs);
  arma::vec neweta = arma::solve(lhs, rhs, solve_opts::likely_sympd);
  etas.head_cols(Kp) = arma::reshape(neweta, d, Kp);
  return;
}

//[[Rcpp::export]]
void updateOneLambdaP (const unsigned int p,
                       const arma::mat& data, const arma::mat& S_mat,
                       arma::cube& lambdas, arma::cube& sigmas,
                       const arma::mat& etas, const arma::mat& alphas, 
                       const arma::mat& gamma, const arma::rowvec& pi) {
  const unsigned int Kp = alphas.n_cols;
  const unsigned int Kc = alphas.n_rows;
  const unsigned int d = data.n_cols;
  const unsigned int N = gamma.n_rows;
  
  arma::uvec pdummy = {p};
  arma::uvec gamma_col_indices = arma::join_cols(pdummy, arma::regspace<uvec>(Kp, Kp + Kc - 1));
  
  arma::rowvec p_g = arma::sum(gamma);
  arma::mat S_g = S_mat * gamma.cols(gamma_col_indices);
  arma::cube S_g_cube(S_g.memptr(), d, d, 1 + Kc);
  
  arma::mat A(d, d, fill::zeros);
  A.diag().fill(-0.5 * p_g(p));
  arma::mat Q = S_g_cube.slice(0);
  arma::mat BB(d, d, fill::zeros);
  
  for (unsigned int c = 0; c < Kc; c++) {
    arma::vec temp1 = sigmas.slice(Kp + c) * etas.col(Kp + c);
    Q += alphas(c, p) * S_g_cube.slice(1 + c);
    BB += (alphas(c, p) * p_g(Kp + c)) * (temp1 * temp1.t() + sigmas.slice(Kp + c));
  }
  
  BB = lambdas.slice(p) * BB * lambdas.slice(p);
  BB += p_g(p) * (etas.col(p) * etas.col(p).t());
  
  arma::mat tobesolved = arma::join_rows(arma::join_cols(A, -Q), arma::join_cols(-BB, -A.t()));
  arma::mat solved_C(d, d, fill::zeros);
  care_inv_C(tobesolved.memptr(), solved_C.memptr(), d); // re-use A matrix
  // arma::mat solved_C = care_inv(A, Q, BB);

  lambdas.slice(p) = (solved_C + solved_C.t())/2;
  sigmas.slice(p) = arma::inv(lambdas.slice(p));
  return;
}

//[[Rcpp::export]]
void updateOneAlphaC (const unsigned int c,
                      const arma::mat& data, const arma::mat& S_mat,
                      const arma::cube& lambdas, const arma::cube& sigmas,
                      const arma::mat& etas, arma::mat& alphas, 
                      const arma::mat& gamma, const arma::rowvec& pi,
                      const unsigned int& nr_maxit, 
                      const double& nr_eps = 1e-16) {
  const unsigned int Kp = alphas.n_cols;
  // const unsigned int Kc = alphas.n_rows;
  const unsigned int d = data.n_cols;
  // const unsigned int n = data.n_rows;
  
  arma::mat Hmat = etas.head_cols(Kp).t();
  arma::mat Lmat = (arma::mat(lambdas.memptr(), d*d, Kp)).t();
  arma::rowvec p_g = arma::sum(gamma);
  
  double Nc = p_g(Kp + c);
  if (Nc == 0) return;
  
  arma::vec yc = data.t() * gamma.col(Kp + c);
  arma::mat Wc = arma::reshape(S_mat * gamma.col(Kp + c), d, d);
  
  arma::vec par = log(arma::clamp(alphas.row(c).t(), 2.2e-16, 1));
  
  double logdetlambda = log(det(lambdas.slice(Kp + c)));
  double before_surrogate_objective = 
    Nc * (0.5 * (logdetlambda - arma::as_scalar(etas.col(Kp + c).t() * sigmas.slice(Kp + c) * etas.col(Kp + c)))) +
    arma::dot(etas.col(Kp + c), yc) - 0.5 * arma::dot(lambdas.slice(Kp + c), Wc);
  if (std::isnan(logdetlambda)) {before_surrogate_objective = -INFINITY;}
  
  for (unsigned int nr_iter = 0; nr_iter < nr_maxit; nr_iter++) {
    arma::vec eval_alpha = arma::normalise(exp(par), 1);
    
    arma::mat Lc(d, d, arma::fill::zeros);
    arma::vec ec(d, arma::fill::zeros);
    
    for (unsigned int i = 0; i < Kp; i++) {
      Lc += eval_alpha[i] * lambdas.slice(i);
      ec += eval_alpha[i] * etas.col(i);
    }
    
    // Precompute some values
    arma::mat Sc = arma::inv(Lc);   // Sc = Sigma_c = Lambda_c inverse 
    arma::vec Sec = Sc * ec;        // Sec = Sigma_c %*% eta_c
    arma::mat tSec = Sec * Sec.t(); // tSec = tcrossprod(Sec)
    double sey = arma::dot(ec, yc); // sey = <eta_c, y_c>
    double sLW = arma::dot(Lc, Wc); // sLW = Tr(Lambda_c %*% W_c)
    // Some complicated scalar value that gets its' own temporary variable.
    double tempscalar = 0.5 * sLW - sey + 0.5 * Nc * arma::dot(ec, Sec) - 0.5 * Nc * d;
    
    arma::vec grad = -eval_alpha % (
      Hmat * (yc - Nc * Sec) +
        0.5 * Lmat * arma::vectorise(Nc * (Sc + tSec) - Wc) + 
        tempscalar
    );
    
    arma::mat hess(Kp, Kp);
    hess.fill(-2 * tempscalar - 0.5 * Nc * d);
    hess.each_col() += 2 * Hmat * (Nc * Sec - yc) + Lmat * arma::vectorise(Wc - Nc * tSec);
    hess += Nc * (2 * Lmat * arma::kron(Sec, Sc) * Hmat.t() -
      Lmat * arma::kron(Sc, tSec + 0.5 * Sc) * Lmat.t() -
      Hmat * Sc * Hmat.t());
    hess = 0.5 * (hess + hess.t());
    hess = -hess % (eval_alpha * eval_alpha.t());
    hess.diag() += grad;
    
    arma::mat nullsp = arma::null(arma::ones(1, Kp) / sqrt(Kp));
    arma:;vec delta_par = arma::solve(nullsp.t() * hess * nullsp, nullsp.t() * grad);
    if (arma::accu(arma::abs(delta_par)) < nr_eps) 
    {
      break;
    }
    par = nullsp * ((nullsp.t() * par) - delta_par);
    par = par - max(par);
  }
  
  // Check if objective improved
  arma::vec eval_alpha = arma::normalise(exp(par), 1);
  
  arma::mat Lc(d, d, arma::fill::zeros); arma::vec ec(d, arma::fill::zeros);
  for (unsigned int i = 0; i < Kp; i++) {Lc += eval_alpha[i] * lambdas.slice(i); ec += eval_alpha[i] * etas.col(i);}
  
  logdetlambda = log(det(Lc));
  double after_surrogate_objective = 
    Nc * (0.5 * (logdetlambda - arma::dot(ec, arma::solve(Lc, ec)))) +
    arma::dot(ec, yc) - 0.5 * arma::dot(Lc, Wc);
  if (std::isnan(logdetlambda)) {after_surrogate_objective = -INFINITY;}
  
  if (after_surrogate_objective < before_surrogate_objective) {
    // Bad! Throw error and don't update alpha_c.
  } else {
    // Good! Update alpha_c.
    alphas.row(c) = eval_alpha.t();
  }
}

//[[Rcpp::export]]
void equilibriateAllAlphaC (arma::cube& lambdas, arma::cube& sigmas,
                            arma::mat& etas, arma::mat& alphas, 
                            const double& eqb_tol = 1e-8) {
  const unsigned int Kp = alphas.n_cols;
  const unsigned int Kc = alphas.n_rows;
  const unsigned int d = etas.n_rows;
  
  arma::mat Vmat(d + d * (d + 1) / 2, Kp);
  for (unsigned int p = 0; p < Kp; p++) {
    
  }
  
  return;
}

/*
//   Vmat = rbind(state$etas[,1:Kp], apply(state$lambdas[,,1:Kp], 3, vech))
//   Vc = rbind(state$etas[,Kp + 1:Kc], apply(state$lambdas[,,Kp + 1:Kc], 3, vech))
//   Amat = cbind(1, diag(Kp), t(Vmat), -t(Vmat))
//   Dmat = diag(Kp)
//   dvec = numeric(Kp)
//   new = sapply(1:Kc, function(c) {
//     bvec = c(1, rep(0, Kp), Vc[,c] - tol, -Vc[,c] - tol)
//     identified = solve.QP(Dmat, dvec, Amat, bvec, 1)
//     return(identified$solution)
//   })
//   return(t(new))
// }
*/

//[[Rcpp::export]]
arma::vec vech (arma::mat x) {
  arma::uvec indices = arma::trimatl_ind(size(x));
  arma::vec vech = x(indices);
  return vech;
}

//[[Rcpp::export]]
void updateLambdaC (arma::cube& lambdas, arma::cube& sigmas,
                    const arma::mat& alphas) {
  unsigned int Kp = alphas.n_cols;
  unsigned int Kc = alphas.n_rows;
  
  for (unsigned int c = 0; c < Kc; c++) {
    lambdas.slice(Kp+c).zeros();
    
    for (unsigned int p = 0; p < Kp; p++)
      lambdas.slice(Kp+c) += alphas(c, p) * lambdas.slice(p);
    
    sigmas.slice(Kp+c) = arma::inv(lambdas.slice(Kp+c));
  }
}

//[[Rcpp::export]]
void updateEtaC (arma::mat& etas, const arma::mat& alphas) {
  unsigned int Kp = alphas.n_cols;
  unsigned int Kc = alphas.n_rows;
  etas.tail_cols(Kc) = etas.head_cols(Kp) * alphas.t();
}

//[[Rcpp::export]]
arma::vec parameterVector (const arma::mat& etas, const arma::cube& lambdas, 
                           const arma::mat& alphas, const arma::rowvec& pi,
                           const arma::mat& data) {
  const unsigned int Kp = alphas.n_cols;
  const unsigned int Kc = alphas.n_rows;
  const unsigned int d = data.n_cols;
  const unsigned int n = data.n_rows;
  
  unsigned int redundant_n_pars = (Kp + Kc) + Kp * (d + d*(d+1)/2) + Kc * (Kp);
  arma::vec V(redundant_n_pars);
  // pi, alphas by c, etas by p, Lambdas by p
  
  // pi
  unsigned int index = 0;
  V.subvec(index, index + Kp + Kc - 1) = pi.t();
  index += Kp + Kc;
  
  // alphas
  for (unsigned int c = 0; c < Kc; c++) {
    V.subvec(index, index + Kp - 1) = alphas.row(c).t();
    index += Kp;
  }
  
  // etas
  for (unsigned int p = 0; p < Kp; p++) {
    V.subvec(index, index + d - 1) = etas.col(p);
    index += d;
  }
  
  arma::uvec upper_indices = arma::trimatu_ind(arma::size(lambdas.slice(0)));

  // Lambdas
  for (unsigned int p = 0; p < Kp; p++) {
    arma::mat lambda = lambdas.slice(p);
    V.subvec(index, index + d*(d+1)/2 - 1) = lambda(upper_indices);
    index += d*(d+1)/2;
  }
  
  return V;
}

// [[Rcpp::export]]
Rcpp::List runEM (unsigned int niter, Rcpp::List state, 
                  const unsigned int nr_maxit = 1,
                  const double nr_eps = 1e-16,
                  const unsigned int& holdziters = 0,
                  const bool show_progress = true) {
  
  arma::mat data     = state["data"];
  arma::mat S_mat    = state["S_mat"];
  arma::mat z_nk     = state["gamma"];
  arma::rowvec pi    = state["Pi"];
  arma::mat etas     = state["etas"];
  arma::cube lambdas = state["lambdas"];
  arma::cube sigmas  = state["sigmas"];
  arma::mat alphas   = state["alphas"];
  const unsigned int Kp = alphas.n_cols;
  const unsigned int Kc = alphas.n_rows;
  
  arma::vec ll_trace(niter, fill::zeros);
  ll_trace.fill(datum::nan);
  arma::vec delta_z_nk_trace(niter, fill::zeros);
  arma::vec delta_V_trace(niter, fill::zeros);
  Progress progress_bar (niter, show_progress);
  
  bool holdzflag = true;
  double loglik = -INFINITY;
  
  arma::vec old_V = parameterVector(etas, lambdas, alphas, pi, data);
  // arma::uvec subtiming(4, arma::fill::zeros);
  
  for (unsigned int iter = 0; iter < niter; iter++) {
    // auto stage1 = std::chrono::high_resolution_clock::now();
    
    arma::mat old_z_nk = z_nk;
    if (iter >= holdziters) {holdzflag = false;}
    updateE(etas, lambdas, sigmas, data, S_mat, pi, z_nk, holdzflag, loglik);
    ll_trace(iter) = loglik;
    if (std::isnan(loglik)) {stop("Parameters became NaN");}
    // auto stage2 = std::chrono::high_resolution_clock::now();
    
    updateAllEtaP(data, lambdas, sigmas, etas, alphas, z_nk);
    updateEtaC(etas, alphas);
    if (!holdzflag) updateE(etas, lambdas, sigmas, data, S_mat, pi, z_nk, holdzflag, loglik);
    // auto stage3 = std::chrono::high_resolution_clock::now();
    
    for (unsigned int p = 0; p < Kp; p++) {
      updateOneLambdaP(p, data, S_mat, lambdas, sigmas, etas, alphas, z_nk, pi);
      updateLambdaC(lambdas, sigmas, alphas);
      if (!holdzflag) updateE(etas, lambdas, sigmas, data, S_mat, pi, z_nk, holdzflag, loglik);
    }
    // auto stage4 = std::chrono::high_resolution_clock::now();
    
    for (unsigned int c = 0; c < Kc; c++) {
      updateOneAlphaC(c, data, S_mat, lambdas, sigmas, etas, alphas, z_nk, pi, nr_maxit, nr_eps);
      updateEtaC(etas, alphas);
      updateLambdaC(lambdas, sigmas, alphas);
      if (!holdzflag) updateE(etas, lambdas, sigmas, data, S_mat, pi, z_nk, holdzflag, loglik);
    }
    // auto stage5 = std::chrono::high_resolution_clock::now();
    
    double delta_z_nk = arma::accu(arma::abs(z_nk - old_z_nk));
    delta_z_nk_trace(iter) = delta_z_nk;
    
    arma::vec new_V = parameterVector(etas, lambdas, alphas, pi, data);
    double delta_V = arma::mean(abs(old_V - new_V));
    old_V = new_V;
    delta_V_trace(iter) = delta_V;
    
    progress_bar.increment();
    
    // subtiming(0) += std::chrono::duration_cast<std::chrono::microseconds>(stage2 - stage1).count();
    // subtiming(1) += std::chrono::duration_cast<std::chrono::microseconds>(stage3 - stage2).count();
    // subtiming(2) += std::chrono::duration_cast<std::chrono::microseconds>(stage4 - stage3).count();
    // subtiming(3) += std::chrono::duration_cast<std::chrono::microseconds>(stage5 - stage4).count();
  }
  
  // subtiming.print("Subtimings");
  
  state["gamma"] = z_nk;
  state["Pi"] = pi;
  state["etas"] = etas;
  state["lambdas"] = lambdas;
  state["sigmas"] = sigmas;
  state["alphas"] = alphas;
  state["ll_trace"] = arma::join_vert(as<arma::vec>(state["ll_trace"]), ll_trace);
  state["gamma_trace"] = arma::join_vert(as<arma::vec>(state["gamma_trace"]), delta_z_nk_trace);
  state["delta_V_trace"] = arma::join_vert(as<arma::vec>(state["delta_V_trace"]), delta_V_trace);
  return state;
}