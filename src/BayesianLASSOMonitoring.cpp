// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#define ARMA_WARN_LEVEL 1
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// simple example of creating two matrices and
// returning the result of an operatioon on them
//
// via the exports attribute we tell Rcpp to make this function
// available from R
//

using namespace Rcpp;
using namespace arma;

arma::colvec rinvgaussiancpp(int n, double mu, double lambda){
  Rcpp::Environment pkg = Rcpp::Environment::namespace_env("VGAM");
  
 // Picking up Matrix() function from Matrix package
  Rcpp::Function rinvgaussian = pkg["rinv.gaussian"];
  
  Rcpp::NumericVector  tmp = rinvgaussian(n, mu, lambda);
 arma::colvec out = Rcpp::as<arma::colvec>(tmp); 
  return out;
}

// [[Rcpp::export]]
arma::mat getV(arma::colvec Y, int q) {
  int T = Y.n_elem;
 arma::mat out(T, q);
  out.zeros();
  
  int tmp;
  
  int i;
  int j;
  for (i = 1; i < T; i++) {
    
    if (i > q) {
      tmp = q;
    } else {
      tmp = i;
    }
    
    for (j = 0; j < tmp; j++) {
      out(i, j) = Y(i - j - 1);
    }
  }
  return(out);
}

arma::colvec rmvnorm(arma::colvec Mean,arma::mat Sigma) {
  
  int q = Mean.n_elem;
 arma::colvec Z =arma::randn(q);
 arma::colvec out = Mean +arma::chol(Sigma) * Z;
  return(out);
  
}

//' Absolute-value-constrained normal distribution
//' 
//' gets a sample from a normal distribution whose absolute observations are constrained.
//'
//' @param n is sample size.
//' @param a is the lower bound in absolute value.
//' @param b is the upper bound in absolute value.
//' @param mean is a mean.
//' @param sd is a standard deviation.
//' @export
//' @examples
//' rtwosegnorm(10, 1, 2, 0, 1)
// [[Rcpp::export]]
arma::colvec rtwosegnorm(int n, double a, double b, double mean, double sd) {
  arma::colvec out(n);
  arma::colvec U =arma::randu(n);
   
   double a1 = (-a - mean) / sd;
   double b1 = (-b - mean) / sd;
   double a2 = (a - mean) / sd;
   double b2 = (b - mean) / sd;
   
   double p1 =arma::normcdf(a1) -arma::normcdf(b1);
   double p2 =arma::normcdf(b2) -arma::normcdf(a2);
   double p = p1 + p2;
   
   int i;
   if (p > 0) {
     for (i = 0; i < n; i++) {
       if (U(i) <= p1 / p) {
         out(i) = R::qnorm5(U(i) * p +arma::normcdf(b1), 0.0, 1.0, 1, 0) * sd + mean;
       } else {
         out(i) = R::qnorm5(U(i) * p +arma::normcdf(a2) - p1, 0.0, 1.0, 1, 0) * sd + mean;
       }
     }
     
   }
   return(out);
   
 }



// [[Rcpp::export]]
arma::mat getGMat(int T, int q) {
 arma::mat tmp(T, T);
  tmp.eye();
  
 arma::mat out(T - q, T);
  out = tmp.submat(q, 0, T - 1, T - 1);
  return(out);
}

// [[Rcpp::export]]
arma::mat getPhiMat(arma::colvec Phi, int T) {
  int q = Phi.n_elem;
 arma::mat tmp(T, T);
  tmp.zeros();
  
  int i;
  for (i = 0; i < q; i++) {
    tmp.diag(-i - 1).fill(Phi(i));
  }
  
 arma::mat out(T - q, T);
  out = tmp.submat(q, 0, T - 1, T - 1);
  return(out);
}

//' Design Matrix (MT)
//' 
//' gets a design matrix as that in MT
//'
//' @param T is length of a process.
//' @param q is the number of lags.
//' @export
//' @examples
//' getHMatMT(100, 5)
// [[Rcpp::export]]
arma::mat getHMatMT(int T, int q) {
  arma::mat tmp(T, T);
   tmp.ones();
  arma::mat L =arma::trimatl(tmp);
   
  arma::mat out(T, T - q);
   out = L.submat(0, q, T - 1, T - 1);
   return(out);
   
 }

//' Design Matrix for Sustained Shift (CM)
//' 
//' gets a design matrix for sustained shift as that in CM
//'
//' @param T is length of a process.
//' @param q is the number of lags.
//' @export
//' @examples
//' getHMatSustained(100, 5)
// [[Rcpp::export]]
arma::mat getHMatSustained(int T, int q) {
   int w = 1;
  arma::mat tmp(T, T);
   tmp.ones();
  arma::mat L =arma::trimatl(tmp);
  arma::mat out = L.cols(q + w, T - 2 + 1);
   int nn = T - 2 - w + 1 - q;
  arma::mat out1(T, nn);
   int rr = 0;
   int kk = 0;
   int i;
   for (i = 0; i < nn; i++) {
     
     if (rr == 0) {
       out1.col(kk) = out.col(i);
       kk = kk + 1;
       rr = rr - w;
     }
     rr = rr + 1;
   }
   int tmpn = floor(nn / w);
   return(out1.cols(0, tmpn - 1));
 }

//' Design Matrix for Isolated Shift (CM)
//' 
//' gets a design matrix for isolated shift as that in CM
//'
//' @param T is length of a process.
//' @param q is the number of lags.
//' @export
//' @examples
//' getHMatIsolated(100, 5)
// [[Rcpp::export]]
arma::mat getHMatIsolated(int T, int q) {
   int w = 1;
  arma::mat tmp(T, T);
   
   int i;
   for (i = 0; i < w; i++) {
     tmp.diag(-i).ones();
   }
   
  arma::mat out = tmp.cols(q, T - 1);
   int nn = T - 1 - q + 1;
  arma::mat out1(T, nn);
   int rr = 0;
   int kk = 0;
   
   for (i = 0; i < nn; i++) {
     
     if (rr == 0) {
       out1.col(kk) = out.col(i);
       kk = kk + 1;
       rr = rr - w;
     }
     rr = rr + 1;
   }
   int tmpn = floor(nn / w);
   return(out1.cols(0, tmpn - 1));
 }

//' gets a design matrix for gradual shift
//'
//' @param T is length of a process.
//' @param q is the number of lags.
//' @export
//' @examples
//' getHMatGradual(100, 5)
// [[Rcpp::export]]
arma::mat getHMatGradual(int T, int q) {
  arma::mat tmp(T, T);
   int w = 1;
   int i;
   for (i = 0; i < T; i++) {
     tmp.diag(-i).fill(i);
   }
   
  arma::mat out = tmp.cols(q, T - 2);
   int nn = T - 2 - q + 1;
  arma::mat out1(T, nn);
   int rr = 0;
   int kk = 0;
   
   for (i = 0; i < nn; i++) {
     
     if (rr == 0) {
       out1.col(kk) = out.col(i);
       kk = kk + 1;
       rr = rr - w;
     }
     rr = rr + 1;
   }
   int tmpn = floor(nn / w);
   return(out1.cols(0, tmpn - 1));
   
 }

//' gets a design matrix for Fourier series
//'
//' @param T is length of a process.
//' @param s is the number of period
//' @param n is the number of Fourier series
//' @export
//' @examples
//' getXSeasonalityFS(100, 15, 10)
// [[Rcpp::export]]
arma::mat getXSeasonalityFS(int T, double s, int n) {
   double pi = 3.141592653589793238462643383280;
   
  arma::mat X(T, 1);
  arma::mat out(T, 2 * n);
   
   int i;
   for (i = 0; i < T; i++) {
     X(i, 0) = i + 1.0;
   }
   
   int k = 0;
   int j;
   for (j = 0; j < n; j++) {
     out.col(k) =arma::cos(2.0 * pi * (j + 1.0) * X / s);
     k = k + 1;
     
     out.col(k) =arma::sin(2.0 * pi * (j + 1.0) * X / s);
     k = k + 1;
   }
   return(out);
 }

arma::mat getInv(arma::mat A) {
  int nn = A.n_cols;
 arma::mat out(nn, nn);
  bool success = true;
  
  success =arma::inv_sympd(out, A);
  if (success == false) {
    success =arma::inv(out, A);
    if (success == false) {
      success =arma::pinv(out, A);
    }
  } 
  
  return(out);
}

/////////////////////////////

arma::mat removeRow(arma::mat A, int a) {
  
 arma::mat out(A.n_rows - 1, A.n_cols);
  int n = A.n_rows;
  int j = 0;
  for (int i = 0; i < n; i++) {
    if (i != a) {
      out.row(j) = A.row(i);
      j = j + 1;
    }
  }
  return(out);
} 


arma::mat removeCol(arma::mat A, int a) {
  
 arma::mat out(A.n_rows, A.n_cols - 1);
  int n = A.n_cols;
  int j = 0;
  for (int i = 0; i < n; i++) {
    if (i != a) {
      out.col(j) = A.col(i);
      j = j + 1;
    }
  }
  return(out);
} 


arma::mat removeCols(arma::mat A,arma::uvec a) {
  
  int mm = a.n_elem;
  
 arma::mat out(A.n_rows, A.n_cols - mm);
  int n = A.n_cols;
  int j = 0;
  int k = 0;
  int i = 0;
  double tmp;
  for (i = 0; i < n; i++) {
    //Rcpp::Rcout << k << std::endl;
    //Rcpp::Rcout << a(k) << std::endl;
    tmp = i + 0.0;
    if (tmp != a(k)) {
      out.col(j) = A.col(i);
      j = j + 1;
    } else {
      if (k + 1 < mm - 1) {
        k = k + 1;
      } else {
        k = mm - 1;
      }
    }
  }
  return(out);
} 


arma::mat removeRows(arma::mat A,arma::uvec a) {
  
  int mm = a.n_elem;
  
 arma::mat out(A.n_rows - mm, A.n_cols);
  int n = A.n_rows;
  int j = 0;
  int k = 0;
  int i = 0;
  double tmp;
  for (i = 0; i < n; i++) {
    tmp = i + 0.0;
    if (tmp != a(k)) {
      out.row(j) = A.row(i);
      j = j + 1;
    } else {
      if (k + 1 < mm - 1) {
        k = k + 1;
      } else {
        k = mm - 1;
      }
    }
  }
  return(out);
} 


Rcpp::List arimacpp(arma::colvec Y, int q){
  Rcpp::Environment pkg = Rcpp::Environment::namespace_env("stats");
  
 // Picking up Matrix() function from Matrix package
  Rcpp::Function arima = pkg["arima"];
  
  Rcpp::NumericVector ord = Rcpp::NumericVector::create(q, 0, 0);
  
  return arima(Y, Rcpp::Named("order") = ord, Rcpp::Named("method") = "CSS");
}

Rcpp::List arimaxcpp(arma::colvec Y, int q, arma::mat X) {
  Rcpp::Environment pkg = Rcpp::Environment::namespace_env("stats");
  
  // Picking up Matrix() function from Matrix package
  Rcpp::Function arima = pkg["arima"];
  
  Rcpp::NumericVector ord = Rcpp::NumericVector::create(q, 0, 0);
  
  return arima(Y, Rcpp::Named("order") = ord, Rcpp::Named("xreg") = X, Rcpp::Named("method") = "CSS");
}

// [[Rcpp::export]]
Rcpp::List initGibbs(arma::colvec Y, arma::mat X, int q, int m, 
                     double tol, Rcpp::String method){
  int T = Y.n_rows;
  int betap = X.n_cols;
  Rcpp::List model = arimaxcpp(Y, q, X);
  Rcpp::NumericVector tmpcoef = model["coef"];
 
  //Rcpp::Rcout << 1 << std::endl;
 
  arma::colvec coef(q + betap);
  
  int i;
  for (i = 0; i < q; i++) {
    coef(i) = tmpcoef(i);
  }
  
  double mu0 = tmpcoef(q);
  
  for (i = (q + 1); i < (q + betap + 1); i++) {
    coef(i - 1) = tmpcoef(i);
  }
  
  arma::colvec Phi = coef.rows(0, q - 1);
  arma::colvec Beta = coef.rows(q, q + betap - 1);
  
  //Rcpp::Rcout << 2 << std::endl;
  
  double sigma2 = model["sigma2"];
  arma::colvec Tau(m);
  Tau.zeros();
  arma::colvec Gamma(m);
  Gamma.fill(tol);
  arma::colvec Mu(T);
  Mu.fill(mu0);
  Mu = Mu + X * Beta;
  double pho = tol;
  arma::colvec eta2 = arma::pow(coef, 2);
  arma::colvec lambda2(q + betap);
  
  //Rcpp::Rcout << 3 << std::endl;
  
  int gg;
  if ((method == "LASSO")) {
    lambda2.fill(pow(q * sqrt(sigma2) /arma::accu(arma::abs(coef)), 2));
  } else if ((method == "ALASSO")) {
    for (gg = 0; gg < (q + betap); gg++) {
      lambda2(gg) = pow((sqrt(sigma2) / abs(coef(gg))), 2);
    }
  }
  
  //Rcpp::Rcout << 4 << std::endl;
  
  Rcpp::NumericVector resi = model["residuals"];
  
  Rcpp::Rcout << 5 << std::endl;
  
  arma::colvec residuals = Rcpp::as<arma::colvec>(resi);
  
  Rcpp::Rcout << 6 << std::endl;
  
  residuals.shed_rows(0, q - 1);
  
  Rcpp::Rcout << 7 << std::endl;
  
  Rcpp::List out1;
  
  out1 = Rcpp::List::create(
    Rcpp::_["Phi"] = Phi,
    Rcpp::_["Beta"] = Beta,
    Rcpp::_["sigma2"] = sigma2,
    Rcpp::_["Tau"] = Tau,
    Rcpp::_["Gamma"] = Gamma,
    Rcpp::_["mu0"] = mu0,
    Rcpp::_["Mu"] = Mu,
    Rcpp::_["pho"] = pho,
    Rcpp::_["eta2"] = eta2,
    Rcpp::_["lambda2"] = lambda2,
    Rcpp::_["residuals"] = residuals
  );
  
  return(out1);
  
}

arma::mat checkSym(arma::mat S) {
 arma::mat out;
  if (S.is_symmetric() == false) {
    S = (S.t() + S) / 2;
  }
  out = S;
  return(out);
}

arma::colvec updatePhi(arma::mat V,arma::mat Vas, 
                      arma::mat A,arma::colvec oldPhi, 
                       double sigma2,arma::mat inveta2mat, 
                       double bound0, double boundqplus1,
                       int phimono, Rcpp::String method) {
  
  //int n = V.n_rows;
  int q = Vas.n_cols;
  
 // Initialize 
 arma::mat tVasVas(q, q);
 arma::mat invtVasVas(q, q);
 arma::mat Phihat(q, 1);
 arma::mat Phi = oldPhi;
 arma::mat S(q, q);
 arma::mat tmpS(q, q);
 arma::colvec M(q);
  
 arma::mat Sgg(1, q);
 arma::mat Snotgg(q - 1, q);
 arma::mat Snotnot(q - 1, q - 1);
 arma::mat invSnotgg(q - 1, q - 1);
  
 arma::mat tmpMat(1, 1); 
 arma::mat tmpSS(q, q);
  int gg;
  
  double mj;
  double sj;
  
  double bound1;
  double bound2;
  
 arma::colvec tmp; 
  
  //Rcpp::Rcout << Vas << std::endl;
  
 // Get tV V
  tVasVas = Vas.t() * Vas;
  //Rcpp::Rcout << tVasVas << std::endl;
  
 // Get Phi hat
  invtVasVas = getInv(tVasVas);
  Phihat = invtVasVas* Vas.t() * V;
  
 // update Phi
  
    if (method == "MT") {
      tmpSS = tVasVas / sigma2 + A;
      tmpSS = checkSym(tmpSS);
      S = getInv(tmpSS);
      S = checkSym(S);
      M = S * ((tVasVas / sigma2) * Phihat);
    } else if (method == "regression") {
      tmpSS = tVasVas + A;
      tmpSS = checkSym(tmpSS);
      tmpS = getInv(tmpSS);
      S = tmpS * sigma2;
      S = checkSym(S);
      M = tmpS * (tVasVas * Phihat);
      
    } else if ((method == "LASSO") || (method == "ALASSO")) {
      tmpSS = tVasVas + inveta2mat;
      tmpSS = checkSym(tmpSS);
      tmpS = getInv(tmpSS);
      S = tmpS * sigma2;
      S = checkSym(S);
      M = tmpS * (tVasVas * Phihat);
    }
  
    
  if (phimono == 0) {
    Phi = rmvnorm(M, S);
  } else {
    for (gg = 0; gg < q; gg++) {
        Sgg = S.row(gg);
        Snotgg = removeRow(S, gg);
        Snotnot = removeCol(Snotgg, gg);
        invSnotgg = getInv(Snotnot);
        
        tmpMat = M(gg) + removeCol(Sgg, gg) * invSnotgg * (removeRow(Phi, gg) - removeRow(M, gg));
        mj = tmpMat(0);
        
        tmpMat = Sgg(gg) - removeCol(Sgg, gg) * invSnotgg * Snotgg.col(gg);
        sj = tmpMat(0);
        
        if (gg == 0) {
          bound1 = abs(bound0);
          bound2 = abs(Phi(1));
        } else if (gg == q - 1) {
          bound1 = abs(Phi(q - 2));
          bound2 = abs(boundqplus1);
        } else {
          bound1 = abs(Phi(gg - 1));
          bound2 = abs(Phi(gg + 1));
        }
        
        if (bound2 < bound1) {
          tmp = rtwosegnorm(1, bound2, bound1, mj, sqrt(sj));
        } else {
          tmp = rtwosegnorm(1, boundqplus1, bound1, mj, sqrt(sj));
        }
        
        Phi(gg) = tmp(0);
      }
  }
    
  return(Phi);
}

arma::colvec updateBeta(arma::mat Y, arma::mat X, arma::mat A, 
                        arma::mat H, arma::mat Tau, arma::mat Gamma,
                        double mu0, double sigma2, arma::mat inveta2mat, 
                        Rcpp::String method, arma::mat D, int Hflg) {
  
  //int n = V.n_rows;
  int p = X.n_cols;
  
 // Initialize 
 arma::mat resiY = Y - mu0;
 if (Hflg == 1) {
   resiY = resiY - H * (Tau % Gamma);
 }
 
 arma::mat tXD = X.t() * D;
 arma::mat tXDX = tXD * X;
 arma::mat invtXDX = getInv(tXDX);
 
 arma::colvec Betahat = invtXDX * tXD * (resiY);
 
 arma::mat S(p, p);
 arma::mat tmpS(p, p);
 arma::colvec M(p);
  
 arma::mat tmpSS(p, p);

 arma::colvec tmp; 
  
  // update Phi
  
    if (method == "MT") {
      tmpSS = tXDX / sigma2 + A;
      tmpSS = checkSym(tmpSS);
      S = getInv(tmpSS);
      S = checkSym(S);
      M = S * ((tXDX / sigma2) * Betahat);
    } else if (method == "regression") {
      tmpSS = tXDX + A;
      tmpSS = checkSym(tmpSS);
      tmpS = getInv(tmpSS);
      S = tmpS * sigma2;
      S = checkSym(S);
      M = tmpS * (tXDX * Betahat);
      
    } else if ((method == "LASSO") || (method == "ALASSO")) {
      tmpSS = tXDX + inveta2mat;
      tmpSS = checkSym(tmpSS);
      tmpS = getInv(tmpSS);
      S = tmpS * sigma2;
      S = checkSym(S);
      M = tmpS * (tXDX * Betahat);
    }
    
    arma::colvec Beta = rmvnorm(M, S);
  
  return(Beta);
}


double updateSigma2(arma::mat resi, arma::colvec phi, arma::colvec beta, 
                    arma::mat invphieta2mat, arma::mat invbetaeta2mat, 
                    int T, int phiq, arma::mat phiA, int betap, arma::mat betaA, 
                    double a, double b, Rcpp::String method, int Xflg) {
  
  int p = beta.n_elem;
  int q = phi.n_elem;
  double sigma2 = 0.0;
  
  arma::mat tResiResi = resi.t() * resi;
  double tmptResiResi = tResiResi(0);
  
  //Rcpp::Rcout << "tmptResiResi:" << tmptResiResi << std::endl;
  
  arma::mat tCoefVarCoef; 
  double tmptCoefVarCoef;
  
  double tmpa = (T - q) / 2.0 + a;
  double tmpb = tmptResiResi / 2.0 + b;
  
  
  arma::mat tmpcoef; 
  
  // update phi contribution
  if (method == "regression") {
    tCoefVarCoef = phi.t() * getInv(phiA) * phi;
    tmptCoefVarCoef = tCoefVarCoef(0);
    tmpa = tmpa + q / 2.0;
    tmpb = tmpb + tmptCoefVarCoef / 2.0;
  } else if ((method == "LASSO") || (method == "ALASSO")) {
    tCoefVarCoef = phi.t() * invphieta2mat * phi;
    tmptCoefVarCoef = tCoefVarCoef(0);
    tmpa = tmpa + q / 2.0;
    tmpb = tmpb + tmptCoefVarCoef / 2.0;
  }
  
  // update beta contribution
  if (Xflg == 1) {
    if (method == "regression") {
      tCoefVarCoef = beta.t() * getInv(betaA) * beta;
      tmptCoefVarCoef = tCoefVarCoef(0);
      tmpa = tmpa + p / 2.0;
      tmpb = tmpb + tmptCoefVarCoef / 2.0;
    } else if ((method == "LASSO") || (method == "ALASSO")) {
      tCoefVarCoef = beta.t() * invbetaeta2mat * beta;
      tmptCoefVarCoef = tCoefVarCoef(0);
      tmpa = tmpa + p / 2.0;
      tmpb = tmpb + tmptCoefVarCoef / 2.0;
    }
  }
  
  sigma2 = 1.0 / R::rgamma(tmpa, 1.0 / (tmpb));
  
  return(sigma2);
}

arma::mat updateinveta2(arma::colvec coef, double sigma2, arma::mat lambda2, int q, double tol) {
 arma::mat coefm =arma::conv_to<arma::mat>::from(coef); 
 arma::mat inveta2(q, 1);
 arma::mat muPrime =arma::sqrt(sigma2 * (lambda2 %arma::pow(coefm, -2)));
 arma::colvec tmp; 
  int gg;
  
  double tmpmuPrime;
  double tmplambda2;
  
  for (gg = 0; gg < q; gg++) {
    tmpmuPrime = muPrime(gg, 0);
    tmplambda2 = lambda2(gg, 0);
    //tmp = rrinvgauss(1, tmpmuPrime, tmplambda2);
    if ((-tol < coefm(gg, 0)) && (coefm(gg, 0) < tol)) {
      inveta2(gg, 0) = 1 / R::rgamma(1.0 / 2.0, 2.0 / tmplambda2);
    } else {
      tmp = rinvgaussiancpp(1, tmpmuPrime, tmplambda2);
      inveta2(gg, 0) = tmp(0);
    }
    
    if (inveta2(gg, 0) < tol) {
      inveta2(gg, 0) = tol;
    }
    
  }
  return(inveta2);
}

arma::mat updatelambda2(arma::mat eta2, int q, 
                        double lambda2alpha, double lambda2beta, Rcpp::String method){
 arma::mat lambda2(q, 1); 
  double shape;
  double scale;
 arma::vec tmpsum; 
  int gg;
 // update lambda2
  if (method == "LASSO") {
    shape = lambda2alpha + q;
    tmpsum = arma::sum(eta2);
    scale = lambda2beta + tmpsum(0)/2.0;
    lambda2.fill(R::rgamma(shape, scale));
  } else if (method == "ALASSO") {
    shape = lambda2alpha + 1;
    for (gg = 0; gg < q; gg++) {
      scale = lambda2beta + eta2(gg)/2.0;
      lambda2(gg) = R::rgamma(shape, scale);
    }
  }
  return(lambda2);
}

Rcpp::List updateTauGamma(arma::colvec Y, arma::mat X, arma::colvec Beta,
                          arma::mat Tau, arma::mat Gamma, 
                          double mu0, double sigma2, double pho, double xi2,
                          int T, int q,arma::mat D, arma::mat H_, int Xflg, int m){
  
 arma::mat Ht(T, 1); 
 arma::mat DHt(T, m);
 arma::mat tHtDHt(m, 1); 
 arma::mat Hnot(T, m - 1);
 arma::mat tmpDHt;
 arma::mat tmptHtDHt; 
  
 arma::mat Taunot(m - 1, 1);
 arma::mat Taut(1, 1);
 arma::mat Gammanot(m - 1, 1);
 arma::mat Gammat(1, 1);
  
 arma::mat zetanot(T, 1);
 arma::mat zetat(T, 1);
  
 arma::mat pvec(m, 1); 
  
 arma::mat muGamma = Gamma;
 arma::mat sigma2Gamma = Gamma;
  
  double tmpzetanot;
  double tmpzetat;
 arma::mat tmp;
  double p;
  
 arma::mat Gammathat(m, 1);
  double st;
  double mt; 
  arma::mat tmpresi = Y - mu0;
  
  if (Xflg == 1) {
   tmpresi = tmpresi - X * Beta;
  }
  
  
  int jj;
  
    
   // get DHt and tHtDHt
    
    //for (jj = 0; jj < m; jj++){
   //  Ht = H_.col(jj);
   //  DHt.col(jj) = D * Ht;
   //  tmp = Ht.t() * DHt.col(jj);
   //  tHtDHt(jj) = tmp(0);
    //}
    
   // update Tau and Gamma
    for (jj = 0; jj < m; jj++) {
      Hnot = removeCol(H_, jj);
      Ht = H_.col(jj);
      
      //tmpDHt = DHt.col(jj);
      //tmptHtDHt = tHtDHt(jj);
      
      //############
      Taunot = removeRow(Tau, jj);
      Taut = Tau.row(jj);
      
      Gammanot = removeRow(Gamma, jj);
      Gammat = Gamma.row(jj);
      
      //update Tau
      zetanot = tmpresi - Hnot * (Taunot % Gammanot);
      zetat = zetanot - Ht * Gammat;
      
      tmp =arma::exp(-1.0 / 2.0 / sigma2 * zetanot.t() * D * zetanot);
      tmpzetanot = tmp(0);
      
      tmp =arma::exp(-1.0 / 2.0 / sigma2 * zetat.t() * D * zetat);
      tmpzetat = tmp(0);
      
      p = pho * tmpzetat / (pho * tmpzetat + (1 - pho) * tmpzetanot);
      
      Tau(jj) = R::rbinom(1, p);
      pvec(jj) = p;
      //############
      if (Tau(jj) == 1) {
        
        DHt.col(jj) = D * Ht;
        tmp = Ht.t() * DHt.col(jj);
        tHtDHt(jj) = tmp(0);
        
        tmpDHt = DHt.col(jj);
        tmptHtDHt = tHtDHt(jj);
        
        //#update Gamma
        tmp = tmpDHt.t() * zetanot / tmptHtDHt;
        Gammathat = tmp(0);
        tmp = 1.0 / (tmptHtDHt / sigma2 + 1 / xi2);
        st = tmp(0);
        tmp = st * (tmptHtDHt * Gammathat) / sigma2;
        mt = tmp(0);
      } else {
        mt = 0;
        st = xi2;
      }
      
      Gamma(jj) = R::rnorm(mt, sqrt(st));
      muGamma(jj) = mt;
      sigma2Gamma(jj) = st;
    }
 
  
  Rcpp::List out = Rcpp::List::create(
    Rcpp::_["Tau"] = Tau,
    Rcpp::_["Gamma"] = Gamma,
    Rcpp::_["p"] = pvec,
    Rcpp::_["muGamma"] = muGamma,
    Rcpp::_["sigma2Gamma"] = sigma2Gamma
  );
  return(out);
}

arma::colvec getresidual(arma::colvec Y, arma::colvec Mu, arma::colvec Phi, int q) {
  int T = Y.n_elem;
  int phiq = q;
  arma::colvec V_ = Y - Mu;
  arma::colvec V = V_.rows(phiq, T - 1);
  arma::mat Vas_ = getV(V_, phiq);
  arma::mat Vas = Vas_.rows(phiq, T - 1);
  arma::colvec tmpresi = V - Vas * Phi;
  return(tmpresi);
}

// [[Rcpp::export]]
Rcpp::List updatepars(arma::colvec Y, Rcpp::List Bset, Rcpp::List oldpars, 
                      Rcpp::Nullable<Rcpp::NumericMatrix> X = R_NilValue, 
                      Rcpp::Nullable<Rcpp::NumericMatrix> H = R_NilValue,
                      Rcpp::Nullable<Rcpp::NumericMatrix> G = R_NilValue) {
  
  //////////////////////////////////////
  // initialize
  
  int T = Y.n_rows;
  
  double tol = Bset["tol"];
  Rcpp::String method = Bset["method"];
  int phimono = Bset["phimono"];
  int phiq = Bset["phiq"];
  arma::mat phiA = Rcpp::as<arma::mat>(Bset["phiA"]);
  double phibound0 = Bset["phibound0"];
  double phiboundqplus1 = Bset["phiboundqplus1"];
  arma::mat betaA = Rcpp::as<arma::mat>(Bset["betaA"]);
  double gammaxi2 = Bset["gammaxi2"];
  double tautheta1 = Bset["tautheta1"];
  double tautheta2 = Bset["tautheta2"];
  double sigma2a = Bset["sigma2a"];
  double sigma2b = Bset["sigma2b"];
  int updatepenalty = Bset["updatepenalty"];
  double lambda2alpha = Bset["lambda2alpha"];
  double lambda2beta = Bset["lambda2beta"];
  
  arma::mat oldPhi = Rcpp::as<arma::mat>(oldpars["Phi"]);
  arma::mat oldBeta = Rcpp::as<arma::mat>(oldpars["Beta"]);
  double oldsigma2 = oldpars["sigma2"];
  arma::mat oldTau = Rcpp::as<arma::mat>(oldpars["Tau"]);
  arma::mat oldGamma = Rcpp::as<arma::mat>(oldpars["Gamma"]);
  double oldmu0 = oldpars["mu0"];
  arma::mat oldMu = Rcpp::as<arma::mat>(oldpars["Mu"]);
  double oldpho = oldpars["pho"];
  arma::mat oldeta2 = Rcpp::as<arma::mat>(oldpars["eta2"]);
  arma::mat oldlambda2 = Rcpp::as<arma::mat>(oldpars["lambda2"]);
  
  
  
  //////////////////////////////////////
  
  arma::mat X_;
  int betap = 0;
  int Xflg = 0;
  
  if (X.isNotNull()) {
    X_ = Rcpp::as<arma::mat>(X);
    betap = X_.n_cols;
    Xflg = 1;
  }
  
  arma::mat H_;
  int m = 0;
  int Hflg = 0;
  
  if (H.isNotNull()) {
    H_ = Rcpp::as<arma::mat>(H);
    m = H_.n_cols;
    Hflg = 1;
  }
  
  arma::mat G_;
  
  if (G.isNotNull()) {
    G_ = Rcpp::as<arma::mat>(G);
  } else {
    G_ = getGMat(T, phiq);
  }
  
  arma::mat oldphieta2 = oldeta2.rows(0, phiq - 1);
  arma::mat oldinvphieta2 = arma::pow(oldphieta2, -1);
  arma::mat oldinvphieta2mat(phiq, phiq);
  oldinvphieta2mat.diag() = oldinvphieta2;
  
  arma::mat oldbetaeta2 = oldeta2.rows(phiq, phiq + betap - 1);
  arma::mat oldinvbetaeta2 = arma::pow(oldbetaeta2, -1);
  arma::mat oldinvbetaeta2mat(betap, betap);
  oldinvbetaeta2mat.diag() = oldinvbetaeta2;
  
  //////////////////////////////////////
  
  int dm = phiq + betap;
  arma::colvec coef(dm);
  
  arma::mat V_;
  arma::mat V;
  arma::mat Vas_;
  arma::mat Vas;
  
  arma::mat PhiMat;
  arma::mat C;
  arma::mat D;
  
  arma::mat Phi;
  arma::mat Beta;
  
  Rcpp::List TauGamma;
  arma::mat Tau;
  arma::mat Gamma;
  
  arma::colvec One(T);
  One.ones();
  
  arma::mat tOneD;
  arma::mat tOneDOne;
  arma::mat invtOneDOne;
  
  arma::mat tmpresi;
  arma::mat tmpmu0hat;
  arma::mat tmpmu0sigma2;
  double mu0;
  
  arma::mat Mu;
  
  double sigma2;
  
  arma::mat eta2;
  arma::mat inveta2;
  
  arma::mat lambda2;
  double pho;
  //////////////////////////////////////
  // update Phi
  V_ = Y - oldMu;
  V = V_.rows(phiq, T - 1);
  Vas_ = getV(V_, phiq);
  Vas = Vas_.rows(phiq, T - 1);
  
  Phi = updatePhi(V, Vas, phiA, 
                    oldPhi, oldsigma2, oldinvphieta2mat, 
                    phibound0, phiboundqplus1,
                    phimono, method);
  
  //////////////////////////////////////
  //Calculate Phi Matrix 
  PhiMat = getPhiMat(Phi, T);
  
  //Calculate C Matrix
  C = G_ - PhiMat;
  
  //Calculate D Matrix
  D = C.t() * C;
  
  //////////////////////////////////////
  // update Beta
  if (Xflg == 1) {
    Beta = updateBeta(Y, X_, betaA, 
                    H_, oldTau, oldGamma,
                    oldmu0, oldsigma2, oldinvbetaeta2mat, 
                    method, D, Hflg);
  }
  
  //////////////////////////////////////
  // update Tau and Gamma
  if (Hflg == 1) {
    TauGamma = updateTauGamma(Y, X_, Beta,
                          oldTau, oldGamma, 
                          oldmu0, oldsigma2, oldpho, gammaxi2,
                          T, phiq, D, H_, Xflg, m);
    
    Tau = Rcpp::as<arma::mat>(TauGamma["Tau"]);
    Gamma = Rcpp::as<arma::mat>(TauGamma["Gamma"]);
  }
  
  //////////////////////////////////////
  // update mu0
  tmpresi = Y;
  if (Xflg == 1) {
    tmpresi = tmpresi - X_ * Beta;
  }
  if (Hflg == 1) {
    tmpresi = tmpresi - H_ * (Tau % Gamma);
  }
  tOneD = One.t() * D;
  tOneDOne = tOneD * One;
  invtOneDOne = arma::pow(tOneDOne, -1);
  
  tmpmu0hat = invtOneDOne * tOneD * tmpresi;
  tmpmu0sigma2 = oldsigma2 * invtOneDOne;
  
  mu0 = R::rnorm(tmpmu0hat(0), sqrt(tmpmu0sigma2(0)));
  
  //////////////////////////////////////
  // update Mu
  Mu = mu0 * One;
  if (Xflg == 1) {
    Mu = Mu + X_ * Beta;
  }
  if (Hflg == 1) {
    Mu = Mu + H_ * (Tau % Gamma);
  }

  //////////////////////////////////////
  // update sigma2
  V_ = Y - Mu;
  V = V_.rows(phiq, T - 1);
  Vas_ = getV(V_, phiq);
  Vas = Vas_.rows(phiq, T - 1);
  tmpresi = V - Vas * Phi;
  
  sigma2 = updateSigma2(tmpresi, Phi, Beta, 
                    oldinvphieta2mat, oldinvbetaeta2mat, 
                    T, phiq, phiA, betap, betaA, 
                    sigma2a, sigma2b, method, Xflg);
  
  //////////////////////////////////////
  // update inveta2
  coef = arma::join_cols(Phi, Beta);
  inveta2 = updateinveta2(coef, sigma2, oldlambda2, phiq + betap, tol);
  eta2 = arma::pow(inveta2, -1);
  
  //////////////////////////////////////
  // update lambda2
  if (updatepenalty == 1) {
    lambda2 = updatelambda2(eta2, phiq + betap, 
                        lambda2alpha, lambda2beta, method);
  }
  
  //////////////////////////////////////
  // update pho
  if (Hflg == 1) {
    pho = R::rbeta(tautheta1 + arma::accu(Tau), tautheta2 + m - arma::accu(Tau));
  }
  
  //////////////////////////////////////
  // update output
  Rcpp::List out;
  
  out = Rcpp::List::create(
    Rcpp::_["Phi"] = Phi,
    Rcpp::_["Beta"] = Beta,
    Rcpp::_["sigma2"] = sigma2,
    Rcpp::_["Tau"] = Tau,
    Rcpp::_["Gamma"] = Gamma,
    Rcpp::_["mu0"] = mu0,
    Rcpp::_["Mu"] = Mu,
    Rcpp::_["pho"] = pho,
    Rcpp::_["eta2"] = eta2,
    Rcpp::_["lambda2"] = lambda2,
    Rcpp::_["residuals"] = tmpresi
  );
  
  return(out);
}

double llf(arma::mat resi, double sigma2) {
  int T = resi.n_rows;
  
  
  double ll = 0.0;
  int i;
  for (i = 0; i < T; i++) {
    ll = ll + R::dnorm4(resi(i), 0, sqrt(sigma2), 1);
  }
  return(ll);
} 

double llYJf(arma::mat Y, arma::mat resi, double sigma2, double theta, int q) {
  arma::mat tmpY = Y;
  tmpY.shed_rows(0, q - 1);
  int T = tmpY.n_rows;
  double ll = llf(resi, sigma2);
  int i;
  for (i = 0; i < T; i++) {
    ll = ll + log(pow(abs(tmpY(i)) + 1, (theta - 1) * sign(tmpY(i))));
  }
  return(ll);
} 

//' Absolute-value-constrained normal distribution
//' 
//' gets a sample from a normal distribution whose absolute observations are constrained.
//'
//' @param n is sample size.
//' @export
//' @examples
//' rtwosegnorm(10, 1, 2, 0, 1)
// [[Rcpp::export]]
arma::colvec invyeojohnsontr(arma::colvec Yyj, double theta, double eps) {
  int T = Yyj.n_elem;
 arma::colvec ans(T);
  
  for (int i = 0; i < T; i++) {
    if ((Yyj(i) >= 0.0) && (abs(theta) > eps)) {
      ans(i) = pow(Yyj(i) * theta + 1.0, 1.0 / theta) - 1.0;
    }
    else if ((Yyj(i) >= 0.0) && (abs(theta) <= eps)) {
      ans(i) = exp(Yyj(i));
    }
    else if ((Yyj(i) < 0.0) && (abs(theta - 2.0) > eps)) {
      ans(i) = 1.0 - pow(-(2.0 - theta) * Yyj(i) + 1.0, 1.0 / (2.0 - theta));
    }
    else if ((Yyj(i) < 0.0) && (abs(theta - 2.0) <= eps)) {
      ans(i) = -exp(-Yyj(i));
    }
  }
  return(ans);
} 

//' Absolute-value-constrained normal distribution
//' 
//' gets a sample from a normal distribution whose absolute observations are constrained.
//'
//' @param n is sample size.
//' @export
//' @examples
//' rtwosegnorm(10, 1, 2, 0, 1)
// [[Rcpp::export]]
arma::colvec yeojohnsontr(arma::colvec Y, double theta, double eps) {
  int T = Y.n_elem;
 arma::colvec ans(T);
  
  for (int i = 0; i < T; i++) {
    if ((Y(i) >= 0) && (abs(theta) > eps)) {
      ans(i) = (pow(Y(i) + 1, theta) - 1) / theta;
    }
    else if ((Y(i) >= 0) && (abs(theta) <= eps)) {
      ans(i) = log(Y(i));
    }
    else if ((Y(i) < 0) && (abs(theta - 2) > eps)) {
      ans(i) = -((pow(-Y(i) + 1, 2 - theta) - 1) / (2 - theta));
    }
    else if ((Y(i) < 0) && (abs(theta - 2) <= eps)) {
      ans(i) = -log(-Y(i));
    }
  }
  
  return(ans);
} 

// [[Rcpp::export]]
arma::mat updateThetaYJ(arma::colvec Y, arma::mat Phi, arma::mat Mu, double sigma2, 
                        double oldtheta, int burnin, int nsim, double tol) {
  
  //double pi = 3.14159265359;
  
  //int T = Y.n_elem;
  int q = Phi.n_rows;
  
  double u;
  double A;
  
  double oldtheta_ = oldtheta;
  arma::colvec oldresi = yeojohnsontr(Y, oldtheta_, tol);
  oldresi = getresidual(oldresi, Mu, Phi, q);
  double oldll = llYJf(Y, oldresi, sigma2, oldtheta_, q);
  
  double thetaas = oldtheta;
  arma::colvec newresi = oldresi;
  double newll = oldll;
    
  arma::mat thetaout(nsim, 1); 
  
  double tmp; 
  double pd;
  
  int i;
  int j = 0;
  
  for (i = 0; i < (nsim + burnin); i++) {
    u = R::runif(0.0, 1.0);
    thetaas = R::rnorm(oldtheta_, 0.1);
    newresi = yeojohnsontr(Y, thetaas, tol);
    newresi = getresidual(newresi, Mu, Phi, q);
    newll = llYJf(Y, newresi, sigma2, thetaas, q);

    tmp = newll - oldll + (R::dnorm4(thetaas, 1, 0.1, 1) - R::dnorm4(oldtheta_, 1, 0.1, 1));
    
    pd = exp(tmp);
    
    A = std::min(1.0, pd);
    
    if (u < A) {
      oldtheta_ = thetaas;
      oldll = newll;
    } 
    
    if (i >= burnin) {
      thetaout(j, 0) = oldtheta_;
      j = j + 1;
    }
  }
  
  return(thetaout);
  
}

//' Absolute-value-constrained normal distribution
//' 
//' gets a sample from a normal distribution whose absolute observations are constrained.
//'
//' @param n is sample size.
//' @export
//' @examples
//' rtwosegnorm(10, 1, 2, 0, 1)
// [[Rcpp::export]]
arma::colvec rtrnorm(int n, double mean, double sd, double lower, double upper) {
 arma::colvec out(n);
 arma::colvec U =arma::randu(n);
  
  double slower = (lower - mean) / sd;
  double supper = (upper - mean) / sd;
  
  double Z = R::pnorm5(supper, 0.0, 1.0, 1, 0) - R::pnorm5(slower, 0.0, 1.0, 1, 0);
  
  int i;
  
  for (i = 0; i < n; i++) {
    out(i) = R::qnorm5(U(i) * Z + R::pnorm5(slower, 0.0, 1.0, 1, 0), 0.0, 1.0, 1, 0) * sd + mean;
  }
  
  return(out);
  
}

double dtrnorm(double x, double mean, double sd, double lower, double upper) {
  double beta = (upper - mean) / sd;
  double alpha = (lower - mean) / sd;
  double z = R::pnorm5(beta, 0.0, 1.0, 1, 0) - R::pnorm5(alpha, 0.0, 1.0, 1, 0);
  double xi = (x - mean) / sd;
  double out = 0.0;
  if ((alpha <= xi) && (xi <= beta)) {
    out = R::dnorm4(xi, 0, 1, 0) / sd / z;
  } 
  
  return(out);
}

// [[Rcpp::export]]
double updateZt(arma::colvec Y, arma::mat Phi, arma::mat Mu, double sigma2, 
                        double theta, int leftcensoring, int rounding, arma::colvec oldZ,
                        int t, double tol, double lb, double ub) {
  
  //double pi = 3.14159265359;
  
  //int T = Y.n_elem;
  
  int q = Phi.n_rows;
  
  double u;
  double A;
  
  double oldZt = oldZ(t);
  arma::colvec oldZ_ = oldZ;
  arma::colvec oldYZ = Y + oldZ_;
  arma::colvec oldresi = yeojohnsontr(oldYZ, theta, tol);
  oldresi = getresidual(oldresi, Mu, Phi, q);
  double oldll = llYJf(Y, oldresi, sigma2, theta, q);
  
  double newZt = oldZt;
  arma::colvec newZ = oldZ;
  arma::colvec newYZ = oldYZ;
  arma::colvec newresi = oldresi;
  double newll = oldll;
    
  arma::mat tmp; 
  double pd;
  
  //int i;
  //int j = 0;
  
    u = R::runif(0.0, 1.0);
    tmp = rtrnorm(1, oldZt, 0.1, lb, ub);
    newZt = tmp(0);
    newZ(t) = newZt;
    
    newYZ = Y + newZ;
    newresi = yeojohnsontr(newYZ, theta, tol);
    newresi = getresidual(newresi, Mu, Phi, q);
    newll = llYJf(Y, newresi, sigma2, theta, q);
      
    pd = newll + log(dtrnorm(newZt, 0.0, 0.1, lb, ub)) - 
      (oldll + log(dtrnorm(oldZt, 0.0, 0.1, lb, ub))) -
      (log(dtrnorm(newZt, oldZt, 0.1, lb, ub)) - log(dtrnorm(oldZt, newZt, 0.1, lb, ub)));
    
    pd = exp(pd);
    
    A = std::min(1.0, pd);
    
    if (u < A) {
      oldZt = newZt;
      oldll = newll;
    } 
    

  return(oldZt);
  
}

// [[Rcpp::export]]
arma::mat updateZ(arma::colvec Y, arma::mat Phi, arma::mat Mu, double sigma2, 
               double theta, int leftcensoring, int rounding, arma::colvec oldZ,
               int burnin, int nsim, double tol) {
  
  //double pi = 3.14159265359;
  
  int T = Y.n_elem;
  
  //int q = Phi.n_rows;
  
  arma::mat Zout(T, nsim);
  arma::mat tmpZ = oldZ;
  double tmpZt;
  
  int flgl;
  int flgr;
  
  double lbl = 0.0;
  double ubl = 0.0;
  
  double lbr = 0.0;
  double ubr = 0.0;
  
  double lb = 0.0;
  double ub = 0.0;
  
  int i;
  int j = 0;
  int t;
  
  if ((leftcensoring == 1) || (rounding == 1)) {
    for (i = 0; i < (burnin + nsim); i++) {
      for (t = 0; t < T; t++) {
        
        flgl = 0;
        flgr = 0;
        
        if (leftcensoring == 1) {
          if (Y(t) <= 0.0) {
            lbl = (-1.0) * arma::datum::inf;
            ubl = 0.0;
            flgl = 1;
          }
        }
        if (rounding == 1) {
          lbr = -0.5;
          ubr = 0.5;
          flgr = 1;
        }
        if ((flgl == 0) && (flgr == 0)) {
          tmpZt = 0.0;
        } else {
          if ((flgl == 1) && (flgr == 0)) {
            lb = lbl;
            ub = ubl;
          } else if ((flgl == 0) && (flgr == 1)) {
            lb = lbr;
            ub = ubr;
          } else if ((flgl == 1) && (flgr == 1)) {
            lb = lbl;
            ub = ubr;
          }
          tmpZt = updateZt(Y, Phi, Mu, sigma2, 
               theta, leftcensoring, rounding, tmpZ,
               t, tol, lb, ub);
        }
        tmpZ(t) = tmpZt;
      }
      if (i >= burnin) {
        Zout.col(j) = tmpZ;
        j = j + 1;
      }
    }
  }
  
  return(Zout);
  
}