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
                       int MonoFlg, Rcpp::String method) {
  
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
 
  if (MonoFlg == 0) {
    Phi = rmvnorm(M, S);
  } else {
    //if ((method == "MonoLASSO") || (method == "MonoALASSO")) {
      //tmpSS = tVasVas + inveta2mat;
      //tmpSS = checkSym(tmpSS);
      //tmpS = getInv(tmpSS);
      //S = tmpS * sigma2;
      //S = checkSym(S);
      //M = tmpS * (tVasVas * Phihat);
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
        tmp = rtwosegnorm(1, bound2, bound1, mj, sqrt(sj));
        Phi(gg) = tmp(0);
      }
    //}
  }
  return(Phi);
}


double updateSigma2(arma::mat resi,arma::colvec Phi,arma::mat inveta2mat, int T, int q, 
                   arma::mat A, double a, double b, Rcpp::String method) {
  double sigma2 = 0.0;
  
 arma::mat tResiResi = resi.t() * resi;
  double tmptResiResi = tResiResi(0);
  
 arma::mat tPhiVarPhi; 
  double tmptPhiVarPhi;
 // update sigma2
  if (method == "MT") {
    sigma2 = 1.0 / R::rgamma((T - q) / 2.0 + a, 1.0 / (tmptResiResi / 2.0 + b));
  } else if (method == "regression") {
    tPhiVarPhi = Phi.t() * A * Phi;
    tmptPhiVarPhi = tPhiVarPhi(0);
    sigma2 = 1.0 / R::rgamma(T / 2.0 + a, 1.0 / (tmptResiResi / 2.0 + tmptPhiVarPhi / 2.0 + b));
  } else if ((method == "LASSO") || (method == "ALASSO") || (method == "MonoLASSO") || (method == "MonoALASSO")) {
    tPhiVarPhi = Phi.t() * inveta2mat * Phi;
    tmptPhiVarPhi = tPhiVarPhi(0);
    sigma2 = 1.0 / R::rgamma(T / 2.0 + a, 1.0 / (tmptResiResi / 2.0 + tmptPhiVarPhi / 2.0 + b));
  }
  
  return(sigma2);
}

arma::mat updateinveta2(arma::colvec Phi, double sigma2,arma::mat lambda2, int q, double tol) {
 arma::mat Phim =arma::conv_to<arma::mat>::from(Phi); 
 arma::mat inveta2(q, 1);
 arma::mat muPrime =arma::sqrt(sigma2 * (lambda2 %arma::pow(Phim, -2)));
 arma::colvec tmp; 
  int gg;
  
  double tmpmuPrime;
  double tmplambda2;
  
  for (gg = 0; gg < q; gg++) {
    tmpmuPrime = muPrime(gg, 0);
    tmplambda2 = lambda2(gg, 0);
    //tmp = rrinvgauss(1, tmpmuPrime, tmplambda2);
    if ((-tol < Phim(gg, 0)) && (Phim(gg, 0) < tol)) {
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

arma::mat updatelambda2(arma::mat eta2, int q, double alpha, double beta, Rcpp::String method){
 arma::mat lambda2(q, 1); 
  double shape;
  double scale;
 arma::vec tmpsum; 
  int gg;
 // update lambda2
  if ((method == "LASSO") || (method == "MonoLASSO")) {
    shape = alpha + q;
    tmpsum =arma::sum(eta2);
    scale = beta + tmpsum(0)/2.0;
    lambda2.fill(R::rgamma(shape, scale));
  } else if ((method == "ALASSO") || (method == "MonoALASSO")) {
    shape = alpha + 1;
    for (gg = 0; gg < q; gg++) {
      scale = beta + eta2(gg)/2.0;
      lambda2(gg) = R::rgamma(shape, scale);
    }
  }
  return(lambda2);
}

Rcpp::List updateTauGamma(arma::colvec tmpY,arma::colvec Phi,arma::mat Tau,arma::mat Gamma, 
                          double mu0, double sigma2, double pho, double xi2,
                          int T, int q,arma::mat D,arma::mat H_, int Hflg, int m){
  
  //Rcpp::Rcout << "m:" << m << std::endl;
  //Rcpp::Rcout << "Hflg:" << Hflg << std::endl;
  
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
  
  int jj;
  
  //Rcpp::Rcout << 2.1 << std::endl;
  
  if (Hflg == 1) {
    
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
      zetanot = tmpY - Hnot * (Taunot % Gammanot);
      zetat = zetanot - Ht * Gammat;
      
      tmp =arma::exp(-1.0 / 2.0 / sigma2 * zetanot.t() * D * zetanot);
      tmpzetanot = tmp(0);
      tmpzetanot = log(tmpzetanot);
      
      tmp =arma::exp(-1.0 / 2.0 / sigma2 * zetat.t() * D * zetat);
      tmpzetat = tmp(0);
      tmpzetat = log(tmpzetat);
      
      //p = pho * tmpzetat / (pho * tmpzetat + (1 - pho) * tmpzetanot
      p = pho / (pho + (1 - pho) * exp(tmpzetanot - tmpzetat));
      
      Tau(jj) = R::rbinom(1, p);
      pvec(jj) = p;
      
      //Rcpp::Rcout << 2.2 << std::endl;
      
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
      
      //Rcpp::Rcout << 2.3 << std::endl;
      
      Gamma(jj) = R::rnorm(mt, sqrt(st));
      muGamma(jj) = mt;
      sigma2Gamma(jj) = st;
      
      //Rcpp::Rcout << 2.4 << std::endl;
    }
    
  }
  
  //Rcpp::Rcout << 2.5 << std::endl;
  
  Rcpp::List out = Rcpp::List::create(
    Rcpp::_["Tau"] = Tau,
    Rcpp::_["Gamma"] = Gamma,
    Rcpp::_["p"] = pvec,
    Rcpp::_["muGamma"] = muGamma,
    Rcpp::_["sigma2Gamma"] = sigma2Gamma
  );
  return(out);
}


Rcpp::List updateZetaBeta(arma::colvec tmpY,arma::colvec Phi,arma::mat Zeta,arma::mat Beta, 
                          double mu0, double sigma2, double pho, double xi2,
                          int T, int q,arma::mat D,arma::mat X_, int Xflg, int m){
  
  //Rcpp::Rcout << "m:" << m << std::endl;
  //Rcpp::Rcout << "Hflg:" << Hflg << std::endl;
  
 arma::mat Ht(T, 1); 
 arma::mat DHt(T, m);
 arma::mat tHtDHt(m, 1); 
 arma::mat Hnot(T, m - 1);
 arma::mat tmpDHt;
 arma::mat tmptHtDHt; 
  
 arma::mat Zetanot(m - 1, 1);
 arma::mat Zetat(1, 1);
 arma::mat Betanot(m - 1, 1);
 arma::mat Betat(1, 1);
  
 arma::mat zetanot(T, 1);
 arma::mat zetat(T, 1);
  
 arma::mat pvec(m, 1); 
  
 arma::mat muBeta = Beta;
 arma::mat sigma2Beta = Beta;
  
  double tmpzetanot;
  double tmpzetat;
 arma::mat tmp;
  double p;
  
 arma::mat Betathat(m, 1);
  double st;
  double mt; 
  
  int jj;
  
  //Rcpp::Rcout << 2.1 << std::endl;
  
  if (Xflg == 1) {
    
   // get DHt and tHtDHt
    
    //for (jj = 0; jj < m; jj++){
   //  Ht = H_.col(jj);
   //  DHt.col(jj) = D * Ht;
   //  tmp = Ht.t() * DHt.col(jj);
   //  tHtDHt(jj) = tmp(0);
    //}
    
   // update Tau and Gamma
    for (jj = 0; jj < m; jj++) {
      Hnot = removeCol(X_, jj);
      Ht = X_.col(jj);
      
      //tmpDHt = DHt.col(jj);
      //tmptHtDHt = tHtDHt(jj);
      
      //############
      Zetanot = removeRow(Zeta, jj);
      Zetat = Zeta.row(jj);
      
      Betanot = removeRow(Beta, jj);
      Betat = Beta.row(jj);
      
      //update Zeta
      zetanot = tmpY - Hnot * (Zetanot % Betanot);
      zetat = zetanot - Ht * Betat;
      
      tmp =arma::exp(-1.0 / 2.0 / sigma2 * zetanot.t() * D * zetanot);
      tmpzetanot = tmp(0);
      tmpzetanot = log(tmpzetanot);
      
      tmp =arma::exp(-1.0 / 2.0 / sigma2 * zetat.t() * D * zetat);
      tmpzetat = tmp(0);
      tmpzetat = log(tmpzetat);
      
      //p = pho * tmpzetat / (pho * tmpzetat + (1 - pho) * tmpzetanot);
      p = pho / (pho + (1 - pho) * exp(tmpzetanot - tmpzetat));
      
      
      Zeta(jj) = R::rbinom(1, p);
      pvec(jj) = p;
      
      //Rcpp::Rcout << 2.2 << std::endl;
      
      //############
      if (Zeta(jj) == 1) {
        
        DHt.col(jj) = D * Ht;
        tmp = Ht.t() * DHt.col(jj);
        tHtDHt(jj) = tmp(0);
        
        tmpDHt = DHt.col(jj);
        tmptHtDHt = tHtDHt(jj);
        
        //#update Gamma
        tmp = tmpDHt.t() * zetanot / tmptHtDHt;
        Betathat = tmp(0);
        tmp = 1.0 / (tmptHtDHt / sigma2 + 1 / xi2);
        st = tmp(0);
        tmp = st * (tmptHtDHt * Betathat) / sigma2;
        mt = tmp(0);
      } else {
        mt = 0;
        st = xi2;
      }
      
      //Rcpp::Rcout << 2.3 << std::endl;
      
      Beta(jj) = R::rnorm(mt, sqrt(st));
      muBeta(jj) = mt;
      sigma2Beta(jj) = st;
      
      //Rcpp::Rcout << 2.4 << std::endl;
    }
    
  }
  
  //Rcpp::Rcout << 2.5 << std::endl;
  
  Rcpp::List out = Rcpp::List::create(
    Rcpp::_["Zeta"] = Zeta,
    Rcpp::_["Beta"] = Beta,
    Rcpp::_["p"] = pvec,
    Rcpp::_["muBeta"] = muBeta,
    Rcpp::_["sigma2Beta"] = sigma2Beta
  );
  return(out);
}


Rcpp::List updatemu0Mu(arma::colvec Y,arma::mat Tau,arma::mat Gamma, double sigma2,
                      arma::mat One,arma::mat D,arma::mat H_, int Hflg, int T, double tol){
  
 arma::colvec zeta;
 arma::mat HTauGamma;
 arma::mat tOneD = One.t() * D;
 arma::mat tmp = tOneD * One;
  double tOneDOne = tmp(0);
  double mu0hat;
  double sq;
  double mu0;
 arma::mat Mu(T, 1);
  
  if (tOneDOne < tol) {
    tOneDOne = tol;
  }
  
  if (Hflg == 0) {
    zeta = Y;
  } else {
    HTauGamma = H_ * (Tau % Gamma);
    zeta = Y - HTauGamma;
  }
  
  //#cat("tOneDOne:", tOneDOne, "\n")
  tmp = tOneD * zeta / tOneDOne;
  mu0hat = tmp(0);
  sq = sigma2 / tOneDOne;
  mu0 = R::rnorm(mu0hat, sqrt(sq));
  
  if (Hflg == 0) {
    Mu.fill(mu0);
  } else {
    Mu = mu0 + HTauGamma;
  }
  
  Rcpp::List out = Rcpp::List::create(
    Rcpp::_["mu0"] = mu0,
    Rcpp::_["Mu"] = Mu
  ); 
  
  return(out);
}

Rcpp::List updatemu0MuX(arma::colvec Y, arma::mat Zeta, arma::mat Beta, 
                        arma::mat Tau,arma::mat Gamma, 
                        double sigma2, arma::mat One,arma::mat D,
                      arma::mat H_, int Hflg, arma::mat X_, int Xflg, 
                      int T, double tol){
  
 arma::colvec tmpY;
 arma::mat HTauGamma;
 arma::mat XZetaBeta;
 arma::mat tOneD = One.t() * D;
 arma::mat tmp = tOneD * One;
  double tOneDOne = tmp(0);
  double mu0hat;
  double sq;
  double mu0;
 arma::mat Mu(T, 1);
  
  if (tOneDOne < tol) {
    tOneDOne = tol;
  }
  
  //Rcpp::Rcout << 1 << std::endl;
  
  tmpY = Y;
  
  //Rcpp::Rcout << 2 << std::endl;
  
  if (Hflg == 1) {
    HTauGamma = H_ * (Tau % Gamma);
    tmpY = tmpY - HTauGamma;
  }
  
 // Rcpp::Rcout << 3 << std::endl;
  
  //Rcpp::Rcout << Zeta << std::endl;
  //Rcpp::Rcout << Beta << std::endl;
  
  if (Xflg == 1) {
    XZetaBeta = X_ * (Zeta % Beta);
    tmpY = tmpY - XZetaBeta;
  }
  
  //Rcpp::Rcout << 4 << std::endl;
  
  //#cat("tOneDOne:", tOneDOne, "\n")
  tmp = tOneD * tmpY / tOneDOne;
  mu0hat = tmp(0);
  sq = sigma2 / tOneDOne;
  mu0 = R::rnorm(mu0hat, sqrt(sq));
  
  Mu.fill(mu0);
  
  if (Hflg == 1) {
    Mu = Mu + HTauGamma;
  }
  
  if (Xflg == 1) {
    Mu = Mu + XZetaBeta;
  }
  
  Rcpp::List out = Rcpp::List::create(
    Rcpp::_["mu0"] = mu0,
    Rcpp::_["Mu"] = Mu
  ); 
  
  return(out);
}


// [[Rcpp::export]]
Rcpp::List GibbsRFLSMcpp(arma::colvec& Y,int& q, 
                        arma::mat& A, double& a, double& b, double& alpha, double& beta, 
                         double& theta1, double& theta2, double& xi2,
                         Rcpp::String& method, double& bound0, double& boundqplus1,
                         int& nsim, int& by, int& burnin,
                         double& tol, Rcpp::Nullable<Rcpp::NumericMatrix> H = R_NilValue) {
  
  auto start = std::chrono::system_clock::now();
  std::time_t start_time = std::chrono::system_clock::to_time_t(start);
  
  Rcpp::Rcout << "Start training using " << method.get_cstring() << " at " << std::ctime(&start_time) <<  std::endl;
  
  /////////////////////////////////
 arma::mat H_;
  
 // Calculate H
  int Hflg = 1;
  int m = 1;
  if (H.isNotNull()) {
    H_ = Rcpp::as<arma::mat>(H);
    m = H_.n_cols;
  } else {
    Hflg = 0;
  }
  
  int T = Y.n_elem;
  
 // Calculate G
 arma::mat G = getGMat(T, q);
  
 // Initialize ones
 arma::mat One(T, 1);
  One.ones();
  
 // Initialize the output
 arma::mat Phiout(q, nsim);
  Phiout.zeros();
  
 arma::mat sigma2out(1, nsim);
  sigma2out.zeros();
  
 arma::mat Tauout(m, nsim);
  Tauout.zeros();
  
 arma::mat Gammaout(m, nsim);
  Gammaout.zeros();
  
 arma::mat muGammaout(m, nsim);
  muGammaout.zeros();
 arma::mat sigma2Gammaout(m, nsim);
  sigma2Gammaout.zeros();
  
 arma::mat pout(m, nsim);
  pout.zeros();
  
 arma::mat mu0out(1, nsim);
  mu0out.zeros();
  
 arma::mat Muout(T, nsim);
  Muout.zeros();
  
 arma::mat phoout(1, nsim);
  phoout.zeros();
  
 arma::mat eta2out(q, nsim);
  eta2out.zeros();
  
 arma::mat lambda2out(q, nsim);
  lambda2out.zeros();
  
 // Is it mono?
  int MonoFlg = 0;
  if ((method == "MonoLASSO") || (method == "MonoALASSO")) {
    MonoFlg = 1;
  }
  
 // Initialize the learning
  Rcpp::List model0 = arimacpp(Y, q);
  
  Rcpp::NumericVector coef = model0["coef"];
  
  double mu0 = coef[q];
  
 arma::mat Phihat(q, 1);
  Rcpp::NumericMatrix varcoef = model0["var.coef"];
  double tmpphi;
  int ii;
  for (ii = 0; ii < q; ii++) {
    tmpphi = coef[ii];
    Phihat(ii) = tmpphi;
  }
  
 arma::mat Phi = Phihat;
  
  
  double bound1;
  double bound2;
  
  double tmpPhi;
  double tmpPhiVar;
  
 arma::colvec tmp;
  
  int gg;
  if (MonoFlg == 1) {
    for (gg = 0; gg < q; gg++) {
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
      tmpPhi = Phi(gg);
      tmpPhiVar = varcoef(gg, gg);
      if (!((bound2 <= abs(tmpPhi)) && 
          (abs(tmpPhi) <= bound1))) {
        tmp = rtwosegnorm(1, boundqplus1, bound1, 
                          tmpPhi, sqrt(tmpPhiVar));
        Phi(gg) = tmp(0);
      }
    }
  }
  
 arma::mat Mu(T, 1);
  Mu.fill(mu0);
  double sigma2 = model0["sigma2"];
  
 arma::mat Tau(m, 1);
  Tau.zeros();
  
 arma::mat Gamma(m, 1);
  Gamma.zeros();
  
 arma::mat muGamma(m, 1);
  muGamma.zeros();
 arma::mat sigma2Gamma(m, 1);
  sigma2Gamma.zeros();
  
 arma::mat pvec(m, 1);
  pvec.zeros();
  
  double pho = R::rbeta(theta1, theta2);
  
 arma::mat eta2(q, 1);
  eta2.zeros();
  
 arma::mat lambda2(q, 1);
  lambda2.zeros();
  
 arma::mat inveta2(q, 1);
  
  eta2 =arma::pow(Phi, 2);
  inveta2 =arma::pow(eta2, -1);
  
 arma::mat inveta2mat(q, q);
  inveta2mat.diag() = inveta2;
  
  if ((method == "LASSO") || (method == "MonoLASSO")) {
    lambda2.fill(pow(q * sqrt(sigma2) /arma::accu(arma::abs(Phi)), 2));
  } else if ((method == "ALASSO") || (method == "MonoALASSO")) {
    for (gg = 0; gg < q; gg++) {
      lambda2(gg) = pow((sqrt(sigma2) / abs(Phi(gg))), 2);
    }
  }
  
 arma::mat DHt(T, T);
  DHt.zeros();
 arma::mat tHtDHt(T, 1);
  tHtDHt.zeros();
  
  int rr = 0;
  
  int TotalSim = nsim * by + burnin;
  
  //outputseq <- seq(burnin + 1, TotalSim, step)
  
 arma::mat V_(T, 1);
  V_.zeros();
  
 arma::mat V(T - q, 1);
  V.zeros();
  
 arma::mat Vas_(T, q);
  Vas_.zeros();
 arma::mat Vas(T - q, q);
  Vas.zeros();
  
  
 arma::mat VasPhi(T - q, 1);
 arma::mat resi(T - q, 1);
  
  
 arma::mat PhiMat(T - q, T);
 arma::mat C;
 arma::mat D;
  
  Rcpp::List TauGamma; 
  Rcpp::List mu0Mu; 
  
 arma::mat tmpSumTau; 
  
  for (ii = 0; ii < TotalSim; ii++) {
    
    if (ii % 100 == 0) {
      Rcpp::Rcout <<"Training: " << ((ii + 0.0) / (TotalSim + 0.0) * 100.0) << '%' << std::endl;
    }
    
    //update V
    V_ = Y - Mu;
    V = V_.rows(q, T - 1);
    Vas_ = getV(V_, q);
    Vas = Vas_.rows(q, T - 1);
    
    //Rcpp::Rcout << Mu << std::endl;
    //Rcpp::Rcout << V << std::endl;
    
   // update Phi
    Phi = updatePhi(V, Vas, A, 
                    Phi, sigma2, inveta2mat, 
                    bound0, boundqplus1,
                    MonoFlg, method);
    
   // Get residuals
    VasPhi = Vas * Phi;
    resi = V - VasPhi;
    
   // update sigma2
    sigma2 = updateSigma2(resi, Phi, inveta2mat, T, q, 
                          A, a, b, method);
    
   // update eta2
    inveta2 = updateinveta2(Phi, sigma2, lambda2, q, tol);
    eta2 =arma::pow(inveta2, -1);
    inveta2mat.diag() = inveta2;
    
   // update lambda2
    lambda2 = updatelambda2(eta2, q, alpha, beta, method);
    
    ///////////////////////////////////////////////////
    //update the random level shift model
    ///////////////////////////////////////////////////
    
    //Calculate Phi Matrix 
    PhiMat = getPhiMat(Phi, T);
    
    //Calculate C Matrix
    C = G - PhiMat;
    
    //Calculate D Matrix
    D = C.t() * C;
    
    //#update Tau and Gamma
    
    TauGamma = updateTauGamma(Y, Phi, Tau, Gamma, 
                              mu0, sigma2, pho, xi2,
                              T, q, D, H_, Hflg, m);
    
    Tau = Rcpp::as<arma::mat>(TauGamma["Tau"]);
    Gamma = Rcpp::as<arma::mat>(TauGamma["Gamma"]);
    pvec = Rcpp::as<arma::mat>(TauGamma["p"]);
    muGamma = Rcpp::as<arma::mat>(TauGamma["muGamma"]);
    sigma2Gamma = Rcpp::as<arma::mat>(TauGamma["sigma2Gamma"]);
    
    //#update mu0 and Mu
    
    mu0Mu = updatemu0Mu(Y, Tau, Gamma, sigma2,
                        One, D, H_, Hflg, T, tol);
    
    mu0 = mu0Mu["mu0"];
    //Rcpp::Rcout << mu0 << std::endl;
    
    Mu = Rcpp::as<arma::mat>(mu0Mu["Mu"]);
    
    //#update pho
    if (Hflg == 1) {
      tmpSumTau =arma::sum(Tau);
      pho = R::rbeta(theta1 + tmpSumTau(0), theta2 + m - tmpSumTau(0));
    }
    
    if (ii >= burnin) {
      if (ii % by == 0) {
        Phiout.col(rr) = Phi;
        sigma2out.col(rr) = sigma2;
        Tauout.col(rr) = Tau;
        Gammaout.col(rr) = Gamma;
        muGammaout.col(rr) = muGamma;
        sigma2Gammaout.col(rr) = sigma2Gamma;
        pout.col(rr) = pvec;
        mu0out.col(rr) = mu0;
        Muout.col(rr) = Mu;
        phoout.col(rr) = pho;
        eta2out.col(rr) = eta2;
        lambda2out.col(rr) = lambda2;
        rr = rr + 1;
      }
    }
    
  }
  
  /////////////////////////////////
  
  Rcpp::Rcout <<"Training: 100%" << std::endl;
  
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  std::time_t end_time = std::chrono::system_clock::to_time_t(end);
  
  Rcpp::Rcout << "Finished training at " << std::ctime(&end_time)
              << "Elapsed time: " << elapsed_seconds.count() << "s"
              << std::endl;
  
  /////////////////////////////////
  
  Rcpp::List out = Rcpp::List::create(
    _["Phi"] = Phiout,
    _["sigma2"] = sigma2out,
    _["Tau"] = Tauout,
    _["Gamma"] = Gammaout,
    _["muGamma"] = muGammaout,
    _["sigma2Gamma"] = sigma2Gammaout,
    _["p"] = pout,
    _["mu0"] = mu0out,
    _["Mu"] = Muout,
    //_["pho"] = phoout,
    //_["eta2"] = eta2out,
    _["lambda2"] = lambda2out
  );
  
  return(out);
  
}


// [[Rcpp::export]]
Rcpp::List GibbsRFLSMUpdatecpp(arma::colvec Y,int q, 
                        arma::mat A, double a, double b, double alpha, double beta, 
                         double theta1, double theta2, double xi2,
                         Rcpp::String method, int monophi, double bound0, double boundqplus1,
                         int nsim, int by, int burnin,
                         double tol, 
                         Rcpp::Nullable<Rcpp::NumericMatrix> G = R_NilValue,
                        Rcpp::Nullable<Rcpp::List> oldpars = R_NilValue,
                         Rcpp::Nullable<Rcpp::NumericMatrix> H = R_NilValue) {
  
  //auto start = std::chrono::system_clock::now();
  //std::time_t start_time = std::chrono::system_clock::to_time_t(start);
  
  //Rcpp::Rcout << "Start training using " << method.get_cstring() << " at " << std::ctime(&start_time) <<  std::endl;
  
  /////////////////////////////////
 arma::mat H_;
  
 // Calculate H
  int Hflg = 1;
  int m = 1;
  if (H.isNotNull()) {
    H_ = Rcpp::as<arma::mat>(H);
    m = H_.n_cols;
  } else {
    Hflg = 0;
  }
  
  int T = Y.n_elem;
  
 // Calculate G
 arma::mat G_;
  
  if (G.isNotNull()) {
    G_ = Rcpp::as<arma::mat>(G);
  } else {
    G_ = getGMat(T, q);
  }
  
 // Initialize ones
 arma::mat One(T, 1);
  One.ones();
  
 // Initialize the output
 arma::mat Phiout(q, nsim);
  Phiout.zeros();
  
 arma::mat sigma2out(1, nsim);
  sigma2out.zeros();
  
 arma::mat Tauout(m, nsim);
  Tauout.zeros();
  
 arma::mat Gammaout(m, nsim);
  Gammaout.zeros();
  
 arma::mat muGammaout(m, nsim);
  muGammaout.zeros();
 arma::mat sigma2Gammaout(m, nsim);
  sigma2Gammaout.zeros();
  
 arma::mat pout(m, nsim);
  pout.zeros();
  
 arma::mat mu0out(1, nsim);
  mu0out.zeros();
  
 arma::mat Muout(T, nsim);
  Muout.zeros();
  
 arma::mat phoout(1, nsim);
  phoout.zeros();
  
 arma::mat eta2out(q, nsim);
  eta2out.zeros();
  
 arma::mat lambda2out(q, nsim);
  lambda2out.zeros();
  
 // Is it mono?
  int MonoFlg = monophi;
  //if ((method == "MonoLASSO") || (method == "MonoALASSO")) {
  //  MonoFlg = 1;
  //}
  
 // Initialize the learning
  
  Rcpp::List model0;
  Rcpp::NumericVector coef;
  double mu0;
 arma::mat Phihat(q, 1);
  Rcpp::NumericMatrix varcoef;
  double tmpphi;
  int ii;
 arma::mat Phi;
  double bound1;
  double bound2;
  
  double tmpPhi;
  double tmpPhiVar;
  
 arma::colvec tmp;
  int gg;
  
 arma::mat Mu(T, 1);
  double sigma2;
  
 arma::mat Tau(m, 1);
  Tau.zeros();
  
 arma::mat Gamma(m, 1);
  Gamma.zeros();
  
 arma::mat muGamma(m, 1);
  muGamma.zeros();
 arma::mat sigma2Gamma(m, 1);
  sigma2Gamma.zeros();
  
 arma::mat pvec(m, 1);
  pvec.zeros();
  
  double pho;
  
 arma::mat eta2(q, 1);
  eta2.zeros();
  
 arma::mat lambda2(q, 1);
  lambda2.zeros();
  
 arma::mat inveta2(q, 1);
 arma::mat inveta2mat(q, q);
  
  Rcpp::List oldpars_ = Rcpp::as<Rcpp::List>(oldpars);
  
  //Rcpp::Rcout << 1 << std::endl;
  
  if (oldpars.isNotNull()) {
    Phi = Rcpp::as<arma::mat>(oldpars_["Phi"]);
    Mu = Rcpp::as<arma::mat>(oldpars_["Mu"]);
    sigma2 = oldpars_["sigma2"];
    Tau = Rcpp::as<arma::mat>(oldpars_["Tau"]);
    Gamma = Rcpp::as<arma::mat>(oldpars_["Gamma"]);
    pho = oldpars_["pho"];
    
    //Rcpp::Rcout << 2 << std::endl;
    
    if ((method == "LASSO") || (method == "ALASSO") || (method == "MonoLASSO") || (method == "MonoALASSO")) {
      eta2 = Rcpp::as<arma::mat>(oldpars_["eta2"]);
      inveta2 =arma::pow(eta2, -1);
      inveta2mat.diag() = inveta2;
      lambda2 = Rcpp::as<arma::mat>(oldpars_["lambda2"]);
    }
    
    mu0 = oldpars_["mu0"];
  } else {
    model0 = arimacpp(Y, q);
    
    coef = model0["coef"];
    
    mu0 = coef[q];
    
    varcoef = Rcpp::as<Rcpp::NumericMatrix>(model0["var.coef"]);
    for (ii = 0; ii < q; ii++) {
      tmpphi = coef[ii];
      Phihat(ii) = tmpphi;
    }
    
    Phi = Phihat;
    
    if (MonoFlg == 1) {
      for (gg = 0; gg < q; gg++) {
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
        tmpPhi = Phi(gg);
        tmpPhiVar = varcoef(gg, gg);
        if (!((bound2 <= abs(tmpPhi)) && 
            (abs(tmpPhi) <= bound1))) {
          tmp = rtwosegnorm(1, boundqplus1, bound1, 
                            tmpPhi, sqrt(tmpPhiVar));
          Phi(gg) = tmp(0);
        }
      }
    }
    
    Mu.fill(mu0);
    sigma2 = model0["sigma2"];
    
    pho = R::rbeta(theta1, theta2);
    
    eta2 =arma::pow(Phi, 2);
    inveta2 =arma::pow(eta2, -1);
    
    inveta2mat.diag() = inveta2;
    
    if ((method == "LASSO") || (method == "MonoLASSO")) {
      lambda2.fill(pow(q * sqrt(sigma2) /arma::accu(arma::abs(Phi)), 2));
    } else if ((method == "ALASSO") || (method == "MonoALASSO")) {
      for (gg = 0; gg < q; gg++) {
        lambda2(gg) = pow((sqrt(sigma2) / abs(Phi(gg))), 2);
      }
    }
    
  }
  
 arma::mat DHt(T, T);
  DHt.zeros();
 arma::mat tHtDHt(T, 1);
  tHtDHt.zeros();
  
  int rr = 0;
  
  int TotalSim = nsim * by + burnin;
  
  //outputseq <- seq(burnin + 1, TotalSim, step)
  
 arma::mat V_(T, 1);
  V_.zeros();
  
 arma::mat V(T - q, 1);
  V.zeros();
  
 arma::mat Vas_(T, q);
  Vas_.zeros();
 arma::mat Vas(T - q, q);
  Vas.zeros();
  
  
 arma::mat VasPhi(T - q, 1);
 arma::mat resi(T - q, 1);
  
  
 arma::mat PhiMat(T - q, T);
 arma::mat C;
 arma::mat D;
  
  Rcpp::List TauGamma; 
  Rcpp::List mu0Mu; 
  
 arma::mat tmpSumTau; 
  
  for (ii = 0; ii < TotalSim; ii++) {
    
    //if (ii % 100 == 0) {
   //  Rcpp::Rcout <<"Training: " << ((ii + 0.0) / (TotalSim + 0.0) * 100.0) << '%' << std::endl;
    //}
    
    //update V
    V_ = Y - Mu;
    V = V_.rows(q, T - 1);
    Vas_ = getV(V_, q);
    Vas = Vas_.rows(q, T - 1);
    
    //Rcpp::Rcout << Mu << std::endl;
    //Rcpp::Rcout << V << std::endl;
    
   // update Phi
    Phi = updatePhi(V, Vas, A, 
                    Phi, sigma2, inveta2mat, 
                    bound0, boundqplus1,
                    MonoFlg, method);
    
   // Get residuals
    VasPhi = Vas * Phi;
    resi = V - VasPhi;
    
   // update sigma2
    sigma2 = updateSigma2(resi, Phi, inveta2mat, T, q, 
                          A, a, b, method);
    
   // update eta2
    inveta2 = updateinveta2(Phi, sigma2, lambda2, q, tol);
    eta2 =arma::pow(inveta2, -1);
    inveta2mat.diag() = inveta2;
    
   // update lambda2
    lambda2 = updatelambda2(eta2, q, alpha, beta, method);
    
    ///////////////////////////////////////////////////
    //update the random level shift model
    ///////////////////////////////////////////////////
    
    //Calculate Phi Matrix 
    PhiMat = getPhiMat(Phi, T);
    
    //Calculate C Matrix
    C = G_ - PhiMat;
    
    //Calculate D Matrix
    D = C.t() * C;
    
    //#update Tau and Gamma
    
    TauGamma = updateTauGamma(Y, Phi, Tau, Gamma, 
                              mu0, sigma2, pho, xi2,
                              T, q, D, H_, Hflg, m);
    
    Tau = Rcpp::as<arma::mat>(TauGamma["Tau"]);
    Gamma = Rcpp::as<arma::mat>(TauGamma["Gamma"]);
    pvec = Rcpp::as<arma::mat>(TauGamma["p"]);
    muGamma = Rcpp::as<arma::mat>(TauGamma["muGamma"]);
    sigma2Gamma = Rcpp::as<arma::mat>(TauGamma["sigma2Gamma"]);
    
    //#update mu0 and Mu
    
    mu0Mu = updatemu0Mu(Y, Tau, Gamma, sigma2,
                        One, D, H_, Hflg, T, tol);
    
    mu0 = mu0Mu["mu0"];
    //Rcpp::Rcout << mu0 << std::endl;
    
    Mu = Rcpp::as<arma::mat>(mu0Mu["Mu"]);
    
    //#update pho
    if (Hflg == 1) {
      tmpSumTau =arma::sum(Tau);
      pho = R::rbeta(theta1 + tmpSumTau(0), theta2 + m - tmpSumTau(0));
    }
    
    if (ii >= burnin) {
      if (ii % by == 0) {
        Phiout.col(rr) = Phi;
        sigma2out.col(rr) = sigma2;
        Tauout.col(rr) = Tau;
        Gammaout.col(rr) = Gamma;
        muGammaout.col(rr) = muGamma;
        sigma2Gammaout.col(rr) = sigma2Gamma;
        pout.col(rr) = pvec;
        mu0out.col(rr) = mu0;
        Muout.col(rr) = Mu;
        phoout.col(rr) = pho;
        eta2out.col(rr) = eta2;
        lambda2out.col(rr) = lambda2;
        rr = rr + 1;
      }
    }
    
  }
  
  /////////////////////////////////
  
  //Rcpp::Rcout <<"Training: 100%" << std::endl;
  //
  //auto end = std::chrono::system_clock::now();
  //std::chrono::duration<double> elapsed_seconds = end-start;
  //std::time_t end_time = std::chrono::system_clock::to_time_t(end);
  //
  //Rcpp::Rcout << "Finished training at " << std::ctime(&end_time)
 //            << "Elapsed time: " << elapsed_seconds.count() << "s"
 //            << std::endl;
  
  /////////////////////////////////
  
  Rcpp::List out = Rcpp::List::create(
    _["Phi"] = Phiout,
    _["sigma2"] = sigma2out,
    _["Tau"] = Tauout,
    _["Gamma"] = Gammaout,
    //_["muGamma"] = muGammaout,
    //_["sigma2Gamma"] = sigma2Gammaout,
    //_["p"] = pout,
    _["mu0"] = mu0out,
    _["Mu"] = Muout,
    _["pho"] = phoout,
    _["eta2"] = eta2out,
    _["lambda2"] = lambda2out
  );
  
  return(out);
  
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
arma::mat lhf(arma::colvec Y,arma::mat Phi,arma::mat Mu, double sigma2) {
  
  double pi = 3.14159265359;
  
  int q = Phi.n_rows;
  //Rcpp::Rcout << q << std::endl;
  int T = Y.n_elem;
  //Rcpp::Rcout << T << std::endl;
 arma::mat V_;
 arma::mat V; 
 arma::mat Vas_;
 arma::mat Vas;
 arma::mat VasPhi;
 arma::mat resi; 
  
  V_ = Y - Mu;
  V = V_.rows(q, T - 1);
  Vas_ = getV(V_, q);
  Vas = Vas_.rows(q, T - 1);
  
  //Rcpp::Rcout << V << std::endl;
  //Rcpp::Rcout << Vas << std::endl;
  
  VasPhi = Vas * Phi;
  
  //Rcpp::Rcout << VasPhi << std::endl;
  
  resi = V - VasPhi;
  
  //Rcpp::Rcout << resi << std::endl;
  
 arma::mat lh = sqrt(1.0 / 2.0 / pi / sigma2) * exp(- 1.0 / 2.0 *arma::pow(resi, 2.0) / sigma2);
  return(lh);
  
}

// [[Rcpp::export]]
arma::colvec boxcoxtr(arma::colvec Y, double theta) {
  
 arma::colvec Ybc = (arma::pow(Y, theta) - 1.0) / theta;
  return(Ybc);
  
}

// [[Rcpp::export]]
arma::colvec invboxcoxtr(arma::colvec Ybc, double theta) {
  
 arma::colvec Y =arma::pow(Ybc * theta + 1.0, 1.0 / theta);
  return(Y);
  
}

// [[Rcpp::export]]
arma::colvec invyeojohnsontr(arma::colvec Yyj, double theta, double eps) {
  int T = Yyj.n_elem;
 arma::colvec ans(T);
  
  for (int i = 0; i < T; i++) {
    if (Yyj(i) >= 0.0 && abs(theta) > eps) {
      ans(i) = pow(Yyj(i) * theta + 1.0, 1.0 / theta) - 1.0;
    }
    else if (Yyj(i) >= 0.0 && abs(theta) <= eps) {
      ans(i) = exp(Yyj(i));
    }
    else if (Yyj(i) < 0.0 && abs(theta - 2.0) > eps) {
      ans(i) = 1.0 - pow(-(2.0 - theta) * Yyj(i) + 1.0, 1.0 / (2.0 - theta));
    }
    else if (Yyj(i) < 0.0 && abs(theta - 2.0) <= eps) {
      ans(i) = -exp(-Yyj(i));
    }
  }
  return(ans);
} 


// [[Rcpp::export]]
arma::colvec yeojohnsontr(arma::colvec Y, double theta, double eps) {
  int T = Y.n_elem;
 arma::colvec ans(T);
  
  for (int i = 0; i < T; i++) {
    if (Y(i) >= 0 && abs(theta) > eps) {
      ans(i) = (pow(Y(i) + 1, theta) - 1) / theta;
    }
    else if (Y(i) >= 0 && abs(theta) <= eps) {
      ans(i) = log(Y(i));
    }
    else if (Y(i) < 0 && abs(theta - 2) > eps) {
      ans(i) = -((pow(-Y(i) + 1, 2 - theta) - 1) / (2 - theta));
    }
    else if (Y(i) < 0 && abs(theta - 2) <= eps) {
      ans(i) = -log(-Y(i));
    }
  }
  
  return(ans);
} 


// [[Rcpp::export]]
arma::mat lhBCf(arma::colvec Y,arma::mat Phi,arma::mat Mu, double sigma2, double theta) {
  
  int T = Y.n_elem;
  int q = Phi.n_rows;
  
 arma::colvec Ybc = boxcoxtr(Y, theta);
 arma::mat lh = lhf(Ybc, Phi, Mu, sigma2);
 arma::mat lhBC = lh %arma::pow(Y.rows(q, T - 1), theta - 1.0);
  
  return(lhBC);
  
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
arma::mat lhYJf(arma::colvec Y,arma::mat Phi,arma::mat Mu, double sigma2, double theta) {
  
  int T = Y.n_elem;
  int q = Phi.n_rows;
  
 arma::colvec Yyj = yeojohnsontr(Y, theta, 0.000001);
 arma::mat lh = lhf(Yyj, Phi, Mu, sigma2);
  
 arma::mat lhYJ = lh;
  
  for (int i = 0; i < (T - q); i++) {
    lhYJ(i) = lhYJ(i) * pow(abs(Y(i + q)) + 1, (theta - 1) * sign(Y(i + q)));
  }
  
  return(lhYJ);
  
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
arma::mat llhYJf(arma::colvec Y,arma::mat Phi,arma::mat Mu, double sigma2, double theta, double eps) {
  
  int T = Y.n_elem;
  int q = Phi.n_rows;
  
 arma::colvec Yyj = yeojohnsontr(Y, theta, eps);
 arma::mat lh = lhf(Yyj, Phi, Mu, sigma2);
  
 arma::mat llhYJ = arma::log(lh);
  
  for (int i = 0; i < (T - q); i++) {
    llhYJ(i) = llhYJ(i) + log(pow(abs(Y(i + q)) + 1, (theta - 1) * sign(Y(i + q))));
  }
  
  return(llhYJ);
  
}


// [[Rcpp::export]]
arma::mat thetaBoxCoxMH(arma::colvec Y,arma::mat Phi,arma::mat Mu, double sigma2, 
                  double oldtheta, int burnin, int nsim, double tol) {
  
  double pi = 3.14159265359;
  
  int T = Y.n_elem;
  int q = Phi.n_rows;
  
  double u;
  double oldtheta_ = oldtheta;
  double thetaas;
  double A;
  
 arma::mat thetaout(nsim, 1); 
  
 arma::mat oldlhBC = lhBCf(Y, Phi, Mu, sigma2, oldtheta_);
  //double sumoldllhBC =arma::accu(oldllhBC);
 arma::uvec ind0 =arma::find(oldlhBC <= tol); 
  oldlhBC(ind0).fill(tol);
  
 arma::uvec ind1; 
  
 arma::mat newlhBC; 
  //double sumnewllhBC;
  
 arma::mat tmp; 
  double pd;
  
  int i;
  int j = 0;
  
  for (i = 0; i < (nsim + burnin); i++) {
    u = R::runif(0.0, 1.0);
    thetaas = R::rnorm(oldtheta_, 0.1);
    newlhBC = lhBCf(Y, Phi, Mu, sigma2, thetaas);
    ind1 =arma::find(newlhBC <= tol);
    newlhBC(ind1).fill(tol);
    
    tmp =arma::cumprod(newlhBC %arma::pow(oldlhBC, -1));
    tmp = tmp * (1 / sqrt(1.0 / 2.0 / pi) * exp(- 1.0 / 2.0 * thetaas * thetaas)) /
      (1 / sqrt(1.0 / 2.0 / pi) * exp(- 1.0 / 2.0 * oldtheta_ * oldtheta_));
      
    //Rcpp::Rcout << tmp << std::endl;
    pd = tmp(T - q - 1);
    //Rcpp::Rcout << pd << std::endl;
    A = std::min(1.0, pd);
    //Rcpp::Rcout << tmp(T - q - 1) << std::endl;
    //Rcpp::Rcout << A << std::endl;
    
    if (u < A) {
      oldtheta_ = thetaas;
      oldlhBC = newlhBC;
    } 
    
    //Rcpp::Rcout << oldtheta_ << std::endl;
    
    if (i >= burnin) {
      thetaout(j, 0) = oldtheta_;
      j = j + 1;
    }
  }
  
  return(thetaout);
  
}

// [[Rcpp::export]]
arma::mat thetaYeoJohnsonMH(arma::colvec Y,arma::mat Phi,arma::mat Mu, double sigma2, 
                        double oldtheta, int burnin, int nsim, double tol) {
  
  double pi = 3.14159265359;
  
  int T = Y.n_elem;
  int q = Phi.n_rows;
  
  double u;
  double oldtheta_ = oldtheta;
  double thetaas;
  double A;
  
 arma::mat thetaout(nsim, 1); 
  
 arma::mat oldlhYJ = lhYJf(Y, Phi, Mu, sigma2, oldtheta_);
  //double sumoldllhBC =arma::accu(oldllhBC);
 arma::uvec ind0 =arma::find(oldlhYJ <= tol); 
  oldlhYJ(ind0).fill(tol);
  
 arma::uvec ind1; 
  
 arma::mat newlhYJ; 
  //double sumnewllhBC;
  
 arma::mat tmp; 
  double pd;
  
  int i;
  int j = 0;
  
  for (i = 0; i < (nsim + burnin); i++) {
    u = R::runif(0.0, 1.0);
    thetaas = R::rnorm(oldtheta_, 0.1);
    newlhYJ = lhYJf(Y, Phi, Mu, sigma2, thetaas);
    ind1 =arma::find(newlhYJ <= tol);
    newlhYJ(ind1).fill(tol);
    
    tmp =arma::cumprod(newlhYJ %arma::pow(oldlhYJ, -1));
    tmp = tmp * (1 / sqrt(1.0 / 2.0 / pi) * exp(- 1.0 / 2.0 * thetaas * thetaas)) /
      (1 / sqrt(1.0 / 2.0 / pi) * exp(- 1.0 / 2.0 * oldtheta_ * oldtheta_));
    
    //Rcpp::Rcout << tmp << std::endl;
    pd = tmp(T - q - 1);
    //Rcpp::Rcout << pd << std::endl;
    A = std::min(1.0, pd);
    //Rcpp::Rcout << tmp(T - q - 1) << std::endl;
    //Rcpp::Rcout << A << std::endl;
    
    if (u < A) {
      oldtheta_ = thetaas;
      oldlhYJ = newlhYJ;
    } 
    
    //Rcpp::Rcout << oldtheta_ << std::endl;
    
    if (i >= burnin) {
      thetaout(j, 0) = oldtheta_;
      j = j + 1;
    }
  }
  
  return(thetaout);
  
}

// [[Rcpp::export]]
arma::mat updatethetaYJMH(arma::colvec Y,arma::mat Phi,arma::mat Mu, double sigma2, 
                        double oldtheta, int burnin, int nsim, double tol) {
  
  double pi = 3.14159265359;
  
  int T = Y.n_elem;
  int q = Phi.n_rows;
  
  double u;
  double oldtheta_ = oldtheta;
  double thetaas;
  double A;
  
 arma::mat thetaout(nsim, 1); 
  
 arma::mat oldllhYJ = llhYJf(Y, Phi, Mu, sigma2, oldtheta_, tol);

 arma::mat newllhYJ = oldllhYJ; 

 arma::mat tmp; 
  double pd;
  
  int i;
  int j = 0;
  
  for (i = 0; i < (nsim + burnin); i++) {
    u = R::runif(0.0, 1.0);
    thetaas = R::rnorm(oldtheta_, 0.1);
    newllhYJ = llhYJf(Y, Phi, Mu, sigma2, thetaas, tol);
    tmp = arma::accu(newllhYJ - oldllhYJ);
    tmp = tmp + R::dnorm4(thetaas, 1, 0.025, 1) - R::dnorm4(oldtheta_, 1, 0.025, 1);
      //(log(1 / sqrt(1.0 / 2.0 / pi) * exp(- 1.0 / 2.0 * thetaas * thetaas)) - 
      //log(1 / sqrt(1.0 / 2.0 / pi) * exp(- 1.0 / 2.0 * oldtheta_ * oldtheta_)));
    pd = exp(tmp(0));
    
    A = std::min(1.0, pd);
    
    if (u < A) {
      oldtheta_ = thetaas;
      oldllhYJ = newllhYJ;
    } 
    
    if (i >= burnin) {
      thetaout(j, 0) = oldtheta_;
      j = j + 1;
    }
  }
  
  return(thetaout);
  
}

// [[Rcpp::export]]
Rcpp::List updateKappaTheta(arma::colvec Y,arma::mat Phi,arma::mat Mu, double sigma2, 
                        double oldkappa, double oldtheta, int burnin, int nsim, double tol) {
  
  double pi = 3.14159265359;
  
  int T = Y.n_elem;
  int q = Phi.n_rows;
  
  double u;
  double oldtheta_ = oldtheta;
  double thetaas;
  arma::mat tmpmat;
  double A;
  
 arma::mat thetaout(nsim, 1); 
  
 arma::mat oldllhYJ = llhYJf(Y, Phi, Mu, sigma2, oldtheta_, tol);
 arma::mat oldllhYJ1 = llhYJf(Y, Phi, Mu, sigma2, 1.0, tol);
 
 double tmp = arma::accu(oldllhYJ1 - oldllhYJ);
 double p = 0.5 / (0.5 + (1 - 0.5) * exp(tmp));
 
 double Kappa = R::rbinom(1, p);
 
 if (Kappa == 1.0) {
   tmpmat = updatethetaYJMH(Y, Phi, Mu, sigma2, 
                         oldtheta, 0, 1, tol);
   thetaas = tmpmat(0);
 } else {
   thetaas = R::rnorm(1, 0.1);
 }
 
 Rcpp::List out;
 out = Rcpp::List::create(
   Rcpp::_["Kappa"] = Kappa,
   Rcpp::_["Theta"] = thetaas
 );
  
 return(out);
  
}

// [[Rcpp::export]]
Rcpp::List GibbsRFLSMBoxCoxcpp(arma::colvec& Y,int& q, 
                        arma::mat& A, double& a, double& b, double& alpha, double& beta, 
                         double& theta1, double& theta2, double& xi2,
                         Rcpp::String& method, double& bound0, double& boundqplus1,
                         int updateBC, double& theta,
                         int& nsim, int& by, int& burnin,
                         double& tol, 
                         Rcpp::Nullable<Rcpp::NumericMatrix> G = R_NilValue,
                         Rcpp::Nullable<Rcpp::List> oldpars = R_NilValue,
                         Rcpp::Nullable<Rcpp::NumericMatrix> H = R_NilValue) {
  
 arma::mat H_;
  
 // Calculate H
  int m = 1;
  if (H.isNotNull()) {
    m = H_.n_cols;
  } 
  
  int T = Y.n_elem;
  
  int TotalSim = nsim * by + burnin;
  
  Rcpp::List GibbsRFLSMModel; 
 arma::colvec Ybc;
  
  Rcpp::List oldpars_;
  
 arma::mat Phi(q, 1);
  Phi.zeros();
  
 arma::mat Mu(T - q, 1);
  Mu.zeros();
  
  double sigma2 = 1.0;
  double theta_ = theta;
 arma::mat tmp; 
  
  Ybc = boxcoxtr(Y, theta_);
  
  if (oldpars.isNotNull()) {
    oldpars_ = Rcpp::as<Rcpp::List>(oldpars);
    Phi = Rcpp::as<arma::mat>(oldpars_["Phi"]);
    Mu = Rcpp::as<arma::mat>(oldpars_["Mu"]);
    sigma2 = oldpars_["sigma2"];
  } else {
    
    oldpars_ = GibbsRFLSMUpdatecpp(Ybc, q, 
                             A, a, b, alpha, beta, 
                             theta1, theta2, xi2, 
                             method, 0, bound0, boundqplus1, 
                             1, 1, TotalSim / 10, tol,
                             G, oldpars, H);
    Phi = Rcpp::as<arma::mat>(oldpars_["Phi"]);
    Mu = Rcpp::as<arma::mat>(oldpars_["Mu"]);
    sigma2 = oldpars_["sigma2"];
  }
  
  //Rcpp::Rcout << 1 << std::endl;
  
  Rcpp::NumericMatrix GG;
 arma::mat G_;
  
  if (G.isNotNull()) {
    GG = Rcpp::as<Rcpp::NumericMatrix>(G);
  } else {
    G_ = getGMat(T, q);
    GG = Rcpp::as<Rcpp::NumericMatrix>(Rcpp::wrap(G_));
  }
  
  //Rcpp::Rcout << 2 << std::endl;
  
  
  
  //////////////////////////
  
 arma::mat Phiout(q, nsim); 
 arma::mat sigma2out(1, nsim);
 arma::mat Tauout(m, nsim);
 arma::mat Gammaout(m, nsim);
 arma::mat mu0out(1, nsim);
 arma::mat Muout(T, nsim);
 arma::mat phoout(1, nsim);
 arma::mat eta2out(q, nsim);
 arma::mat lambda2out(q, nsim);
 arma::mat thetaout(1, nsim);
  
  //////////////////////////
  
  int i;
  int rr = 0;
  
  for (i = 0; i < TotalSim; i++) {
    
    if (i % 100 == 0) {
      Rcpp::Rcout <<"Training: " << ((i + 0.0) / (TotalSim + 0.0) * 100.0) << '%' << std::endl;
    }
    
    if (updateBC == 1) {
      tmp = thetaBoxCoxMH(Y, Phi, Mu, sigma2, 
                    theta_, 1, 1, tol);
      theta_ = tmp(0);
      //Rcpp::Rcout << theta_ << std::endl;
      Ybc = boxcoxtr(Y, theta_);
    }
    
 //  Rcpp::Rcout << 4 << std::endl;
 //  
    oldpars_ = GibbsRFLSMUpdatecpp(Ybc, q, 
                             A, a, b, alpha, beta, 
                             theta1, theta2, xi2, 
                             method, 0, bound0, boundqplus1, 
                             1, 1, 0, tol,
                             GG, oldpars_, H);
 //  
 //  Rcpp::Rcout << 5 << std::endl;
    
    if (i >= burnin) {
      if (i % by == 0) {
        Phiout.col(rr) = Rcpp::as<arma::mat>(oldpars_["Phi"]);
        sigma2out.col(rr) = Rcpp::as<arma::mat>(oldpars_["sigma2"]);
        Tauout.col(rr) = Rcpp::as<arma::mat>(oldpars_["Tau"]);
        Gammaout.col(rr) = Rcpp::as<arma::mat>(oldpars_["Gamma"]);
        mu0out.col(rr) = Rcpp::as<arma::mat>(oldpars_["mu0"]);
        Muout.col(rr) = Rcpp::as<arma::mat>(oldpars_["Mu"]);
        phoout.col(rr) = Rcpp::as<arma::mat>(oldpars_["pho"]);
        eta2out.col(rr) = Rcpp::as<arma::mat>(oldpars_["eta2"]);
        lambda2out.col(rr) = Rcpp::as<arma::mat>(oldpars_["lambda2"]);
        thetaout.col(rr) = theta_;
        rr = rr + 1;
      }
    }
    
  }
  
  Rcpp::List out; 
  out = Rcpp::List::create(
    _["Phi"] = Phiout,
    _["sigma2"] = sigma2out,
    _["Tau"] = Tauout,
    _["Gamma"] = Gammaout,
    _["mu0"] = mu0out,
    _["Mu"] = Muout,
    _["pho"] = phoout,
    _["eta2"] = eta2out,
    _["lambda2"] = lambda2out,
    _["theta"] = thetaout
  );
  
  return(out);
  
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


arma::colvec getucY(arma::colvec Yyj,arma::colvec Y,arma::mat Phi,arma::mat Mu, double sigma2, double theta, double eps) {
  
  int q = Phi.n_rows;
  //Rcpp::Rcout << q << std::endl;
  int T = Y.n_elem;
  //Rcpp::Rcout << T << std::endl;
 arma::mat V(T, 1); 
 arma::mat Vas(q, 1);
 arma::mat VasPhi;
  
 arma::colvec ucY = Y;
 arma::colvec ucYyj = Yyj; 
  
 arma::mat fit(T - q, 1);
 arma::colvec tmp;
  
  //Rcpp::Rcout << 1 << std::endl;
  
  for (int i = 0; i < T; i++) {
    V(i, 0) = ucYyj(i) - Mu(i, 0);
    
    //Rcpp::Rcout << i << std::endl;
    
    if (i >= q) {
      for (int j = 0; j < q; j++) {
        Vas(j, 0) = V(i - 1 - j, 0);
      }
      VasPhi = Vas.t() * Phi;
      fit(i - q, 0) = Mu(i, 0) + VasPhi(0);
      
      if (Y(i) == 0) { 
        //Rcpp::Rcout << 3.01 << std::endl;
        tmp = rtrnorm(1, fit(i - q, 0), sqrt(sigma2), (-1.0) *arma::math::inf(), 0.0);
        //Rcpp::Rcout << 3.1 << std::endl;
        ucYyj(i) = tmp(0); 
        //Rcpp::Rcout << 3.2 << std::endl;
      }
      //Rcpp::Rcout << ucYyj(i)  << std::endl;
    }
    
    //Rcpp::Rcout << 3 << std::endl;
    //Rcpp::Rcout << ucYyj(i)  << std::endl;
    
    
    //Rcpp::Rcout << 4 << std::endl;
  }
  
  //Rcpp::Rcout << ucYyj << std::endl;
  
  //ucY = invyeojohnsontr(ucYyj, theta, eps);
  
  //Rcpp::Rcout << 5 << std::endl;
  
  return(ucYyj);
  
}



// [[Rcpp::export]]
Rcpp::List GibbsRFLSMYeoJohnsonZcpp(arma::colvec& Y,int& q, 
                              arma::mat& A, double& a, double& b, double& alpha, double& beta, 
                               double& theta1, double& theta2, double& xi2,
                               Rcpp::String& method, double& bound0, double& boundqplus1,
                               int updateYJ, double& theta,
                               int updateZ, double eps,
                               int& nsim, int& by, int& burnin,
                               double& tol, 
                               Rcpp::Nullable<Rcpp::NumericMatrix> G = R_NilValue,
                               Rcpp::Nullable<Rcpp::List> oldpars = R_NilValue,
                               Rcpp::Nullable<Rcpp::NumericMatrix> H = R_NilValue) {
  
  
  auto start = std::chrono::system_clock::now();
  std::time_t start_time = std::chrono::system_clock::to_time_t(start);
  
  Rcpp::Rcout << "Start training using " << method.get_cstring() << " at " << std::ctime(&start_time) <<  std::endl;
  
  ///////////////////////////////////
  
 arma::mat H_;
  
 // Calculate H
  int m = 1;
  if (H.isNotNull()) {
    H_ = Rcpp::as<arma::mat>(H);
    m = H_.n_cols;
  } 
  
  int T = Y.n_elem;
  
  int TotalSim = nsim * by + burnin;
  
  Rcpp::List GibbsRFLSMModel; 
 arma::colvec Yyj = Y;
  
  Rcpp::List oldpars_;
  
 arma::mat Phi(q, 1);
  Phi.zeros();
  
 arma::mat Mu(T - q, 1);
  Mu.zeros();
  
  double sigma2 = 1.0;
  double theta_ = theta;
 arma::mat tmp; 
  
  if (updateYJ == 1) {
    Yyj = yeojohnsontr(Y, theta_, eps);
  }
  
  
  if (oldpars.isNotNull()) {
    oldpars_ = Rcpp::as<Rcpp::List>(oldpars);
    Phi = Rcpp::as<arma::mat>(oldpars_["Phi"]);
    Mu = Rcpp::as<arma::mat>(oldpars_["Mu"]);
    sigma2 = oldpars_["sigma2"];
  } else {
    
    oldpars_ = GibbsRFLSMUpdatecpp(Yyj, q, 
                                   A, a, b, alpha, beta, 
                                   theta1, theta2, xi2, 
                                   method, bound0, boundqplus1, 
                                   1, 1, 0, TotalSim / 10, tol,
                                   G, oldpars, H);
    Phi = Rcpp::as<arma::mat>(oldpars_["Phi"]);
    Mu = Rcpp::as<arma::mat>(oldpars_["Mu"]);
    sigma2 = oldpars_["sigma2"];
  }
  
  if (updateZ == 1) {
    Yyj = getucY(Yyj, Y, Phi, Mu, sigma2, theta_, eps);
  }
  
  Rcpp::NumericMatrix GG;
 arma::mat G_;
  
  if (G.isNotNull()) {
    GG = Rcpp::as<Rcpp::NumericMatrix>(G);
  } else {
    G_ = getGMat(T, q);
    GG = Rcpp::as<Rcpp::NumericMatrix>(Rcpp::wrap(G_));
  }
  
  //Rcpp::Rcout << 2 << std::endl;
  
  
  
  //////////////////////////
  
 arma::mat Phiout(q, nsim); 
 arma::mat sigma2out(1, nsim);
 arma::mat Tauout(m, nsim);
 arma::mat Gammaout(m, nsim);
 arma::mat mu0out(1, nsim);
 arma::mat Muout(T, nsim);
 arma::mat phoout(1, nsim);
 arma::mat eta2out(q, nsim);
 arma::mat lambda2out(q, nsim);
 arma::mat thetaout(1, nsim);
 arma::mat Yyjout(T, nsim);
  
  //////////////////////////
  
  
  //////////////////////////
  
  int i;
  int rr = 0;
  
  for (i = 0; i < TotalSim; i++) {
    
    if (i % 100 == 0) {
      Rcpp::Rcout <<"Training: " << ((i + 0.0) / (TotalSim + 0.0) * 100.0) << '%' << std::endl;
    }
    
   //  Rcpp::Rcout << 4 << std::endl;
   //  
    oldpars_ = GibbsRFLSMUpdatecpp(Yyj, q, 
                                   A, a, b, alpha, beta, 
                                   theta1, theta2, xi2, 
                                   method, 0, bound0, boundqplus1, 
                                   1, 1, 0, tol,
                                   GG, oldpars_, H);
   //  
   //  Rcpp::Rcout << 5 << std::endl;
    
    if (updateYJ == 1) {
      tmp = thetaYeoJohnsonMH(Y, Phi, Mu, sigma2, 
                              theta_, 1, 1, tol);
      
      theta_ = tmp(0);
      //Rcpp::Rcout << theta_ << std::endl;
      Yyj = yeojohnsontr(Y, theta_, eps);
    }
    
    if (updateZ == 1) {
      Yyj = getucY(Yyj, Y, Phi, Mu, sigma2, theta_, eps);
    }
    
    if (i >= burnin) {
      if (i % by == 0) {
        Phiout.col(rr) = Rcpp::as<arma::mat>(oldpars_["Phi"]);
        sigma2out.col(rr) = Rcpp::as<arma::mat>(oldpars_["sigma2"]);
        Tauout.col(rr) = Rcpp::as<arma::mat>(oldpars_["Tau"]);
        Gammaout.col(rr) = Rcpp::as<arma::mat>(oldpars_["Gamma"]);
        mu0out.col(rr) = Rcpp::as<arma::mat>(oldpars_["mu0"]);
        Muout.col(rr) = Rcpp::as<arma::mat>(oldpars_["Mu"]);
        phoout.col(rr) = Rcpp::as<arma::mat>(oldpars_["pho"]);
        eta2out.col(rr) = Rcpp::as<arma::mat>(oldpars_["eta2"]);
        lambda2out.col(rr) = Rcpp::as<arma::mat>(oldpars_["lambda2"]);
        thetaout.col(rr) = theta_;
        Yyjout.col(rr) = Yyj;
        rr = rr + 1;
      }
    }
    
  }
  
  Rcpp::Rcout <<"Training: 100%" << std::endl;
  
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  std::time_t end_time = std::chrono::system_clock::to_time_t(end);
  
  Rcpp::Rcout << "Finished training at " << std::ctime(&end_time)
              << "Elapsed time: " << elapsed_seconds.count() << "s"
              << std::endl;
  
  
  Rcpp::List out; 
  out = Rcpp::List::create(
    _["Phi"] = Phiout,
    _["sigma2"] = sigma2out,
    _["Tau"] = Tauout,
    _["Gamma"] = Gammaout,
    _["mu0"] = mu0out,
    _["Mu"] = Muout,
    _["pho"] = phoout,
    _["eta2"] = eta2out,
    _["lambda2"] = lambda2out,
    _["theta"] = thetaout,
    _["Yyj"] = Yyjout
  );
  
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
arma::mat updateZt(arma::colvec Y, arma::colvec Z, arma::mat Phi,arma::mat Mu, double sigma2, 
                        double theta, int t, int leftcensoring, int rounding, int burnin, double tol) {
  
  double pi = 3.14159265359;
  
  int T = Y.n_elem;
  int q = Phi.n_rows;
  
  double u;
  double A;
  
  arma::mat tmp;
  
  double oldZt = Z(t);
  arma::colvec oldZ = Z;
  double newZt = oldZt;
  arma::colvec newZ = oldZ;
  
  arma::colvec oldYZ = Y + oldZ;
  arma::colvec newYZ = Y + newZ;
  
  arma::mat oldllhYJ = llhYJf(oldYZ, Phi, Mu, sigma2, theta, tol);
  arma::mat newllhYJ = oldllhYJ; 
 
  double lbl = -(1.0) * arma::datum::inf;
  double ubl = 0.0;
  
  double lbr = -0.5;
  double ubr = 0.5;
  
  double lb = 0.0;
  double ub = 0.0;
  
  int flgl = 0;
  int flgr = 0;
  
  if (rounding == 1) {
    lb = lbr;
    ub = ubr;
    flgl = 1;
  }
  
  if (leftcensoring == 1) {
    if (Y(t) <= 0) {
      lb = lbl;
      ub = ubl;
      flgr = 1;
    }
  }
  
  if ((leftcensoring == 1) && ((rounding == 1))) {
    if (Y(t) <= 0) {
      lb = lbl;
      ub = ubr;
      flgl = 1;
      flgr = 1;
    }
  }
  
  arma::mat Ztout(1, 1);
  
  double pd;
  
  int i;
  int j = 0;
  
  if ((flgl == 1) || (flgr == 1)) {
    for (i = 0; i < (1 + burnin); i++) {
      u = R::runif(0.0, 1.0);
      tmp = rtrnorm(1, oldZt, 0.1, lb, ub);
      newZt = tmp(0);
      newZ(t) = newZt;
      newYZ = Y + newZ;
      
      newllhYJ = llhYJf(newYZ, Phi, Mu, sigma2, theta, tol);
      
      tmp = arma::accu(newllhYJ - oldllhYJ);
      tmp = tmp + log(dtrnorm(newZt, 0.0, 0.1, lb, ub)) - log(dtrnorm(oldZt, 0.0, 0.1, lb, ub)) -
        (log(dtrnorm(newZt, oldZt, 0.1, lb, ub)) - log(dtrnorm(oldZt, newZt, 0.1, lb, ub)));
      //tmp = tmp - (log(dtrnorm(newZt, oldZt, 0.1, lb, ub)) - log(dtrnorm(oldZt, newZt, 0.1, lb, ub)));
      pd = exp(tmp(0));
      
      A = std::min(1.0, pd);
      
      if (u < A) {
        oldZt = newZt;
        oldllhYJ = newllhYJ;
      } 
      
      if (i >= burnin) {
        Ztout(j, 0) = oldZt;
        j = j + 1;
      }
    }
  } else {
    Ztout.zeros();
  }
  
  
  
  return(Ztout);
  
}


// [[Rcpp::export]]
arma::mat updateZtMD(arma::colvec Y, arma::colvec Z, arma::mat Phi,arma::mat Mu, double sigma2, 
                        double theta, int t, arma::colvec missingdata, double missingdatalb, double missingdataub, 
                        int burnin, double tol) {
  
  double pi = 3.14159265359;
  
  int T = Y.n_elem;
  int q = Phi.n_rows;
  
  double u;
  double A;
  
  arma::mat tmp;
  
  double oldZt = Z(t);
  arma::colvec oldZ = Z;
  double newZt = oldZt;
  arma::colvec newZ = oldZ;
  
  arma::colvec oldYZ = Y + oldZ;
  arma::colvec newYZ = Y + newZ;
  
  arma::mat oldllhYJ = llhYJf(oldYZ, Phi, Mu, sigma2, theta, tol);
  arma::mat newllhYJ = oldllhYJ; 
 
  arma::mat Ztout(1, 1);
  
  double pd;
  
  int missingdatat = missingdata(t);

  Rcpp::Rcout << "t:" << t << std::endl;
  Rcpp::Rcout << "missingdatat:" << missingdatat << std::endl;
  
  double lb = missingdatalb;
  double ub = missingdataub;
  
  int i;
  int j = 0;
  
  if (missingdatat == 1) {
    for (i = 0; i < (1 + burnin); i++) {
      u = R::runif(0.0, 1.0);
      tmp = rtrnorm(1, oldZt, 0.1, lb, ub);
      newZt = tmp(0);
      newZ(t) = newZt;
      newYZ = Y + newZ;
      
      newllhYJ = llhYJf(newYZ, Phi, Mu, sigma2, theta, tol);
      
      tmp = arma::accu(newllhYJ - oldllhYJ);
      tmp = tmp + log(dtrnorm(newZt, 0.0, 0.1, lb, ub)) - log(dtrnorm(oldZt, 0.0, 0.1, lb, ub)) -
        (log(dtrnorm(newZt, oldZt, 0.1, lb, ub)) - log(dtrnorm(oldZt, newZt, 0.1, lb, ub)));
      //tmp = tmp - (log(dtrnorm(newZt, oldZt, 0.1, lb, ub)) - log(dtrnorm(oldZt, newZt, 0.1, lb, ub)));
      pd = exp(tmp(0));
      
      A = std::min(1.0, pd);
      
      if (u < A) {
        oldZt = newZt;
        oldllhYJ = newllhYJ;
      } 
      
      if (i >= burnin) {
        Ztout(j, 0) = oldZt;
        j = j + 1;
      }
    }
  } else {
    Ztout.zeros();
  }
  
  return(Ztout);
  
}

// [[Rcpp::export]]
arma::mat updateZZ(arma::colvec Y, arma::colvec Z, arma::mat Phi,arma::mat Mu, double sigma2, 
                        double theta, int leftcensoring, int rounding, int burnin, int nsim, double tol) {
  
  int T = Y.n_elem;
  int q = Phi.n_rows;
  
  arma::mat tmpZ = Z;
  arma::mat Zout(T, nsim);
  //Rcpp::Rcout << "Zout:" << Zout << std::endl;
  arma::mat tmpZt(1, 1);
  tmpZt.zeros();
  
  int i = 0;
  int t = 0;
  
  for (i = 0; i < nsim; i++) {
    for (t = 0; t < T; t++) {
      
      tmpZ.row(t) = updateZt(Y, tmpZ, Phi, Mu, sigma2, 
                        theta, t, leftcensoring, rounding, burnin, tol);
    }
    //Rcpp::Rcout << "tmpZ:" << tmpZ << std::endl;
    Zout.col(i) = tmpZ;
    //Rcpp::Rcout << "Zout:" << Zout << std::endl;
  }
  
  return(Zout);
  
}


// [[Rcpp::export]]
arma::mat updateZZMD(arma::colvec Y, arma::colvec Z, arma::mat Phi,arma::mat Mu, double sigma2, 
                        double theta, arma::colvec missingdata, double missingdatalb, double missingdataub, 
                        int burnin, int nsim, double tol) {
  
  int T = Y.n_elem;
  int q = Phi.n_rows;
  
  arma::mat tmpZ = Z;
  arma::mat Zout(T, nsim);
  //Rcpp::Rcout << "Zout:" << Zout << std::endl;
  arma::mat tmpZt(1, 1);
  tmpZt.zeros();
  
  int i = 0;
  int t = 0;
  
  for (i = 0; i < nsim; i++) {
    for (t = 0; t < T; t++) {
      
      tmpZ.row(t) = updateZtMD(Y, tmpZ, Phi, Mu, sigma2, 
                        theta, t, missingdata, missingdatalb, missingdataub, 
                        burnin, tol);
    }
    //Rcpp::Rcout << "tmpZ:" << tmpZ << std::endl;
    Zout.col(i) = tmpZ;
    //Rcpp::Rcout << "Zout:" << Zout << std::endl;
  }
  
  return(Zout);
  
}



// [[Rcpp::export]]
arma::mat updateZ0(arma::colvec Y, arma::colvec Z, arma::mat Phi,arma::mat Mu, double sigma2, 
                        double theta, int leftcensoring, int rounding, int burnin, int nsim, double tol) {
  
  int T = Y.n_elem;
  int q = Phi.n_rows;
  
  arma::mat tmpZ = Z;
  arma::mat Zout(T, nsim);
  //Rcpp::Rcout << "Zout:" << Zout << std::endl;
  arma::mat tmpZt(1, 1);
  tmpZt.zeros();
  
  int i = 0;
  int t = 0;
  
  for (i = 0; i < nsim; i++) {
    for (t = 0; t < q; t++) {
      
      tmpZ.row(t) = updateZt(Y, tmpZ, Phi, Mu, sigma2, 
                        theta, t, leftcensoring, rounding, burnin, tol);
    }
    //Rcpp::Rcout << "tmpZ:" << tmpZ << std::endl;
    Zout.col(i) = tmpZ;
    //Rcpp::Rcout << "Zout:" << Zout << std::endl;
  }
  
  //Zout = Zout.rows(0, q);
  
  return(Zout);
  
}

// [[Rcpp::export]]
Rcpp::List GibbsRFLSMYeoJohnsonZcpp1(arma::colvec& Y,int& q, 
                              arma::mat& A, double& a, double& b, double& alpha, double& beta, 
                               double& theta1, double& theta2, double& xi2,
                               Rcpp::String& method, int monophi, double& bound0, double& boundqplus1,
                               int updateYJ, double& theta,
                               int leftcensoring, int rounding, double eps,
                               int& nsim, int& by, int& burnin,
                               double& tol, 
                               Rcpp::Nullable<Rcpp::NumericMatrix> G = R_NilValue,
                               Rcpp::Nullable<Rcpp::List> oldpars = R_NilValue,
                               Rcpp::Nullable<Rcpp::NumericMatrix> H = R_NilValue) {
  
  
  auto start = std::chrono::system_clock::now();
  std::time_t start_time = std::chrono::system_clock::to_time_t(start);
  
  Rcpp::Rcout << "Start training using " << method.get_cstring() << " at " << std::ctime(&start_time) <<  std::endl;
  
  ///////////////////////////////////
  
 arma::mat H_;
  
 // Calculate H
  int m = 1;
  if (H.isNotNull()) {
    H_ = Rcpp::as<arma::mat>(H);
    m = H_.n_cols;
  } 
  
  int T = Y.n_elem;
  
  int TotalSim = nsim * by + burnin;
  
  Rcpp::List GibbsRFLSMModel; 
 arma::colvec Yyj = Y;
  
  Rcpp::List oldpars_;
  
 arma::mat Phi(q, 1);
  Phi.zeros();
  
 arma::mat Mu(T - q, 1);
  Mu.zeros();
  
  double sigma2 = 1.0;
  double theta_ = theta;
 arma::mat tmp; 
  
  if (updateYJ == 1) {
    Yyj = yeojohnsontr(Y, theta_, eps);
  }
  
  
  if (oldpars.isNotNull()) {
    oldpars_ = Rcpp::as<Rcpp::List>(oldpars);
    Phi = Rcpp::as<arma::mat>(oldpars_["Phi"]);
    Mu = Rcpp::as<arma::mat>(oldpars_["Mu"]);
    sigma2 = oldpars_["sigma2"];
  } else {
    
    oldpars_ = GibbsRFLSMUpdatecpp(Yyj, q, 
                                   A, a, b, alpha, beta, 
                                   theta1, theta2, xi2, 
                                   method, monophi, bound0, boundqplus1, 
                                   1, 1, 0, tol,
                                   G, oldpars, H);
    Phi = Rcpp::as<arma::mat>(oldpars_["Phi"]);
    Mu = Rcpp::as<arma::mat>(oldpars_["Mu"]);
    sigma2 = oldpars_["sigma2"];
  }
  
  //if (updateZ == 1) {
  //  Yyj = getucY(Yyj, Y, Phi, Mu, sigma2, theta_, eps);
  //}
  
  Rcpp::NumericMatrix GG;
 arma::mat G_;
  
  if (G.isNotNull()) {
    GG = Rcpp::as<Rcpp::NumericMatrix>(G);
  } else {
    G_ = getGMat(T, q);
    GG = Rcpp::as<Rcpp::NumericMatrix>(Rcpp::wrap(G_));
  }
  
  //Rcpp::Rcout << 2 << std::endl;
  
  arma::mat Z(T, 1);
  Z.zeros();
  
  //////////////////////////
  
 arma::mat Phiout(q, nsim); 
 arma::mat sigma2out(1, nsim);
 arma::mat Tauout(m, nsim);
 arma::mat Gammaout(m, nsim);
 arma::mat mu0out(1, nsim);
 arma::mat Muout(T, nsim);
 arma::mat phoout(1, nsim);
 
 arma::mat eta2out;
 arma::mat lambda2out;
 if ((method == "LASSO") || (method == "ALASSO")) {
   eta2out.set_size(q, nsim);
   lambda2out.set_size(q, nsim);
 }
 
 arma::mat thetaout;
 if (updateYJ == 1) {
   thetaout.set_size(1, nsim);
 }
 
 arma::mat Zout;
 if ((leftcensoring == 1) || (rounding == 1)) {
   Zout.set_size(T, nsim);
 }
  
  //////////////////////////
  
  
  //////////////////////////
  
  int i;
  int rr = 0;
  
  for (i = 0; i < TotalSim; i++) {
    
    if (i % 100 == 0) {
      Rcpp::Rcout <<"Training: " << ((i + 0.0) / (TotalSim + 0.0) * 100.0) << '%' << std::endl;
    }
    
    //if (updateZ == 1) {
    if ((leftcensoring == 1) || (rounding == 1)) {
      Z = updateZZ(Y, Z, Phi, Mu, sigma2, 
                  theta_, leftcensoring, rounding, 0, 1, tol);
      Yyj = Y + Z;
    }
    
    if (updateYJ == 1) {
      //tmp = thetaYeoJohnsonMH(Y, Phi, Mu, sigma2, 
      //                        theta_, 1, 1, tol);
      tmp = updatethetaYJMH(Yyj, Phi, Mu, sigma2, 
                              theta_, 0, 1, tol);
      theta_ = tmp(0);
      //Rcpp::Rcout << theta_ << std::endl;
      Yyj = yeojohnsontr(Yyj, theta_, eps);
    }
    
   //  Rcpp::Rcout << 4 << std::endl;
   //  
    oldpars_ = GibbsRFLSMUpdatecpp(Yyj, q, 
                                   A, a, b, alpha, beta, 
                                   theta1, theta2, xi2, 
                                   method, monophi, bound0, boundqplus1, 
                                   1, 1, 0, tol,
                                   GG, oldpars_, H);
    
    Phi = Rcpp::as<arma::mat>(oldpars_["Phi"]);
    Mu = Rcpp::as<arma::mat>(oldpars_["Mu"]);
    sigma2 = oldpars_["sigma2"];
   //  
   //  Rcpp::Rcout << 5 << std::endl;
    
    
    
    //if (updateZ == 1) {
    //  Yyj = getucY(Yyj, Y, Phi, Mu, sigma2, theta_, eps);
    //}
    
    if (i >= burnin) {
      if (i % by == 0) {
        Phiout.col(rr) = Rcpp::as<arma::mat>(oldpars_["Phi"]);
        sigma2out.col(rr) = Rcpp::as<arma::mat>(oldpars_["sigma2"]);
        Tauout.col(rr) = Rcpp::as<arma::mat>(oldpars_["Tau"]);
        Gammaout.col(rr) = Rcpp::as<arma::mat>(oldpars_["Gamma"]);
        mu0out.col(rr) = Rcpp::as<arma::mat>(oldpars_["mu0"]);
        Muout.col(rr) = Rcpp::as<arma::mat>(oldpars_["Mu"]);
        phoout.col(rr) = Rcpp::as<arma::mat>(oldpars_["pho"]);
        
        
        if ((method == "LASSO") || (method == "ALASSO")) {
          eta2out.col(rr) = Rcpp::as<arma::mat>(oldpars_["eta2"]);
          lambda2out.col(rr) = Rcpp::as<arma::mat>(oldpars_["lambda2"]);
        }
        
        if (updateYJ == 1) {
          thetaout.col(rr) = theta_;
        }
        //Yyjout.col(rr) = Yyj;
        
        if ((leftcensoring == 1) || (rounding == 1)) {
          Zout.col(rr) = Z;
        }
        rr = rr + 1;
      }
    }
    
  }
  
  Rcpp::Rcout <<"Training: 100%" << std::endl;
  
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  std::time_t end_time = std::chrono::system_clock::to_time_t(end);
  
  Rcpp::Rcout << "Finished training at " << std::ctime(&end_time)
              << "Elapsed time: " << elapsed_seconds.count() << "s"
              << std::endl;
  
  
  Rcpp::List out; 
  out = Rcpp::List::create(
    _["Phi"] = Phiout,
    _["sigma2"] = sigma2out,
    _["Tau"] = Tauout,
    _["Gamma"] = Gammaout,
    _["mu0"] = mu0out,
    _["Mu"] = Muout,
    _["pho"] = phoout,
    _["eta2"] = eta2out,
    _["lambda2"] = lambda2out,
    _["theta"] = thetaout,
    _["Z"] = Zout
    //_["Yyj"] = Yyjout
  );
  
  return(out);
  
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
arma::mat updateZSim(arma::colvec Y, arma::colvec oldZ, arma::mat Phi,arma::mat Mu, double sigma2, 
                  double theta, double eps, int leftcensoring, int rounding) {
  
  double pi = 3.14159265359;
  
  arma::colvec newYZ = Y + oldZ; 
  arma::colvec newZ = oldZ;
  int q = Phi.n_rows;
  //Rcpp::Rcout << q << std::endl;
  int T = Y.n_elem;
  //Rcpp::Rcout << T << std::endl;
 arma::mat V(T, 1); 
 arma::mat Vas(1, q);
 arma::mat MuX(1, 1); 
 arma::mat MuinvhX(1, 1); 
 arma::mat varinhX(1, 1); 
 
 double lb;
 double ub;
 
 int flg;
 arma::mat tmp(1, 1);
 
 int t;
 int j;
 
 //Rcpp::Rcout << 1 << std::endl;
 
 if ((leftcensoring == 1) || (rounding == 1)) {
   for (t = 0; t < T; t++) {
     flg = 0;
     
     V.row(t) = yeojohnsontr(newYZ.row(t), theta, eps) - Mu(t);
     
     //Rcpp::Rcout << 2 << std::endl;
     
     if (t >= q) {
       for (j = 0; j < q; j++) {
         Vas(j) = V(t - j - 1);
       }
       
       MuX = Mu(t) + Vas * Phi;
       
       if (Y(t) <= 0) {
         if ((leftcensoring == 1) && (rounding == 1)) {
            lb = -(1.0) * arma::datum::inf;
           
            tmp.fill(0.5);
            tmp = yeojohnsontr(tmp, theta, eps);
            ub = tmp(0);
            
            flg = 1;
         } else if ((leftcensoring == 0) && (rounding == 1)) {
            tmp.fill(-0.5);
            tmp = yeojohnsontr(tmp, theta, eps);
            lb = tmp(0);
            
            tmp.fill(0.5);
            tmp = yeojohnsontr(tmp, theta, eps);
            ub = tmp(0);
            
            flg = 1;
         } else if ((leftcensoring == 1) && (rounding == 0)) {
            lb = -(1.0) * arma::datum::inf;
            ub = 0.0;
            flg = 1;
         }
       } else {
         if (rounding == 1) {
           tmp.fill(Y(t) - 0.5);
           tmp = yeojohnsontr(tmp, theta, eps);
           lb = tmp(0);
            
           tmp.fill(Y(t) + 0.5);
           tmp = yeojohnsontr(tmp, theta, eps);
           ub = tmp(0);
            
           flg = 1;
         }
       }
       
       if (flg == 1) {
         tmp = rtrnorm(1, MuX(0), sqrt(sigma2), lb, ub);
         tmp = invyeojohnsontr(tmp, theta, eps);
         newYZ(t) = tmp(0);
       } else {
         newYZ(t) = newYZ(t);
       }
       
       newZ(t) = newYZ(t) - Y(t);
     }
   }
   
   
 } 
 
  return(newZ);
  
}



// [[Rcpp::export]]
Rcpp::List GibbsRFLSMXUpdatecpp(arma::colvec Y,int q, 
                        arma::mat A, double a, double b, double alpha, double beta, 
                         double theta1, double theta2, double xi2,
                         Rcpp::String method, int monophi, double bound0, double boundqplus1,
                         int nsim, int by, int burnin,
                         double tol, 
                         Rcpp::Nullable<Rcpp::NumericMatrix> G = R_NilValue,
                        Rcpp::Nullable<Rcpp::List> oldpars = R_NilValue,
                        Rcpp::Nullable<Rcpp::NumericMatrix> X = R_NilValue,
                         Rcpp::Nullable<Rcpp::NumericMatrix> H = R_NilValue) {
  
  //auto start = std::chrono::system_clock::now();
  //std::time_t start_time = std::chrono::system_clock::to_time_t(start);
  
  //Rcpp::Rcout << "Start training using " << method.get_cstring() << " at " << std::ctime(&start_time) <<  std::endl;
  
  /////////////////////////////////
 arma::mat X_;
  
 // Calculate H
  int Xflg = 1;
  int betap = 0;
  if (X.isNotNull()) {
    X_ = Rcpp::as<arma::mat>(X);
    betap = X_.n_cols;
  } else {
    Xflg = 0;
  }
 
 
 arma::mat H_;
  
 // Calculate H
  int Hflg = 1;
  int m = 0;
  if (H.isNotNull()) {
    H_ = Rcpp::as<arma::mat>(H);
    m = H_.n_cols;
  } else {
    Hflg = 0;
  }
  
  int T = Y.n_elem;
  
 // Calculate G
 arma::mat G_;
  
  if (G.isNotNull()) {
    G_ = Rcpp::as<arma::mat>(G);
  } else {
    G_ = getGMat(T, q);
  }
  
  //Rcpp::Rcout << 1 << std::endl;
  
 // Initialize ones
 arma::mat One(T, 1);
  One.ones();
  
 // Initialize the output
 arma::mat Phiout(q, nsim);
  Phiout.zeros();
  
 arma::mat sigma2out(1, nsim);
  sigma2out.zeros();
  
  arma::mat Tauout;
  arma::mat Gammaout;
  if (Hflg == 1) {
    Tauout.set_size(m, nsim);
    Tauout.zeros();
  
    Gammaout.set_size(m, nsim);
    Gammaout.zeros();
  }
 
  arma::mat Zetaout;
  arma::mat Betaout;
  if (Xflg == 1) {
    Zetaout.set_size(betap, nsim);
    Zetaout.zeros();
  
    Betaout.set_size(betap, nsim);
    Betaout.zeros();
  }
  
  //Rcpp::Rcout << 1 << std::endl;
  
 //arma::mat muGammaout(m, nsim);
 // muGammaout.zeros();
 //arma::mat sigma2Gammaout(m, nsim);
 // sigma2Gammaout.zeros();
 // 
 //arma::mat pout(m, nsim);
 // pout.zeros();
  
 arma::mat mu0out(1, nsim);
  mu0out.zeros();
  
 arma::mat Muout(T, nsim);
  Muout.zeros();
  
  arma::mat phoout;
  if (Hflg == 1) {
    phoout.set_size(1, nsim);
    phoout.zeros();
  }
  
  arma::mat phobetaout;
  if (Xflg == 1) {
    phobetaout.set_size(1, nsim);
    phobetaout.zeros();
  }
 
 arma::mat eta2out;
  arma::mat lambda2out;
  if ((method == "LASSO") || (method == "ALASSO")) {
    eta2out.set_size(q, nsim);
    eta2out.zeros();
  
    lambda2out.set_size(q, nsim);
    lambda2out.zeros();
  }
  
  
 
  
 // Is it mono?
  int MonoFlg = monophi;
  //if ((method == "MonoLASSO") || (method == "MonoALASSO")) {
  //  MonoFlg = 1;
  //}
  
  
  //Rcpp::Rcout << 1 << std::endl;
  
 // Initialize the learning
  
  Rcpp::List model0;
  Rcpp::NumericVector coef;
  double mu0;
 arma::mat Phihat(q, 1);
 arma::mat Betahat(betap, 1);
  Rcpp::NumericMatrix varcoef;
  double tmpphi;
  double tmpbeta;
  int ii;
 arma::mat Phi;
  double bound1;
  double bound2;
  
  double tmpPhi;
  double tmpPhiVar;
  
 arma::colvec tmp;
  int gg;
  
 arma::mat Mu(T, 1);
  double sigma2;
  
  arma::mat Tau;
  arma::mat Gamma;
  if (Hflg == 1) {
    Tau.set_size(m, 1);
    Tau.zeros();
  
    Gamma.set_size(m, 1);
    Gamma.zeros();
  }
  
  arma::mat Zeta;
  arma::mat Beta;
  if (Xflg == 1) {
    Zeta.set_size(betap, 1);
    Zeta.zeros();
  
    Beta.set_size(betap, 1);
    Beta.zeros();
  }
 
  
 //arma::mat muGamma(m, 1);
 // muGamma.zeros();
 //arma::mat sigma2Gamma(m, 1);
 // sigma2Gamma.zeros();
 // 
 //arma::mat pvec(m, 1);
 // pvec.zeros();
  
  double pho = 0.0;
  double phobeta = 0.0;
  
 arma::mat eta2(q, 1);
  eta2.zeros();
  
 arma::mat lambda2(q, 1);
  lambda2.zeros();
  
 arma::mat inveta2(q, 1);
 arma::mat inveta2mat(q, q);
  
  Rcpp::List oldpars_ = Rcpp::as<Rcpp::List>(oldpars);
  

  
  //Rcpp::Rcout << 1 << std::endl;
  
  if (oldpars.isNotNull()) {
    Phi = Rcpp::as<arma::mat>(oldpars_["Phi"]);
    Mu = Rcpp::as<arma::mat>(oldpars_["Mu"]);
    sigma2 = oldpars_["sigma2"];
    
    if (Hflg == 1) {
      Tau = Rcpp::as<arma::mat>(oldpars_["Tau"]);
      Gamma = Rcpp::as<arma::mat>(oldpars_["Gamma"]);
      pho = oldpars_["pho"];
    }
    
    if (Xflg == 1) {
      Zeta = Rcpp::as<arma::mat>(oldpars_["Zeta"]);
      Beta = Rcpp::as<arma::mat>(oldpars_["Beta"]);
      phobeta = oldpars_["phobeta"];
    }
    
    //Rcpp::Rcout << 2 << std::endl;
    
    if ((method == "LASSO") || (method == "ALASSO") || (method == "MonoLASSO") || (method == "MonoALASSO")) {
      eta2 = Rcpp::as<arma::mat>(oldpars_["eta2"]);
      inveta2 =arma::pow(eta2, -1);
      inveta2mat.diag() = inveta2;
      lambda2 = Rcpp::as<arma::mat>(oldpars_["lambda2"]);
    }
    
    mu0 = oldpars_["mu0"];
  } else {
 
    //Rcpp::Rcout << 3 << std::endl;
 
      if (Xflg == 1) {
        model0 = arimaxcpp(Y, q, X_);
      } else {
        model0 = arimacpp(Y, q);
      }
      
    //Rcpp::Rcout << 1 << std::endl;
    
      coef = model0["coef"];
      
      mu0 = coef[q];
      
      //Rcpp::Rcout << "coef:" << coef << std::endl;
      //Rcpp::Rcout << "mu0:" << mu0 << std::endl;
      
      varcoef = Rcpp::as<Rcpp::NumericMatrix>(model0["var.coef"]);
      for (ii = 0; ii < q; ii++) {
        tmpphi = coef[ii];
        Phihat(ii) = tmpphi;
      }
      
      Phi = Phihat;
      
      
      if (MonoFlg == 1) {
        for (gg = 0; gg < q; gg++) {
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
          tmpPhi = Phi(gg);
          tmpPhiVar = varcoef(gg, gg);
          if (!((bound2 <= abs(tmpPhi)) && 
              (abs(tmpPhi) <= bound1))) {
            tmp = rtwosegnorm(1, boundqplus1, bound1, 
                              tmpPhi, sqrt(tmpPhiVar));
            Phi(gg) = tmp(0);
          }
        }
      }
      
      //Rcpp::Rcout << "Phi:" << Phi << std::endl;
      
      if (Xflg == 1) {
        for (ii = q + 1; ii < (q + betap + 1); ii++) {
          //Rcpp::Rcout << "ii:" << ii << std::endl;
          
          tmpbeta = coef[ii];
          Betahat(ii - (q + 1)) = tmpbeta;
          
          //Rcpp::Rcout << "tmpbeta:" << tmpbeta << std::endl;
        }
      
        Beta = Betahat;
        Zeta.ones();
      }
      
      //Rcpp::Rcout << "betap:" << betap << std::endl;
      //Rcpp::Rcout << "Beta:" << Beta << std::endl;
      //Rcpp::Rcout << "Zeta:" << Zeta << std::endl;
      
      Mu.fill(mu0);
      
      if (Xflg == 1) {
        Mu = Mu + X_ * (Zeta % Beta);
        phobeta = R::rbeta(theta1, theta2);
      }
      
      //Rcpp::Rcout << 1 << std::endl;
      
      sigma2 = model0["sigma2"];
      
      pho = R::rbeta(theta1, theta2);
      
      eta2 =arma::pow(Phi, 2);
      inveta2 =arma::pow(eta2, -1);
      
      inveta2mat.diag() = inveta2;
      
      if ((method == "LASSO") || (method == "MonoLASSO")) {
        lambda2.fill(pow(q * sqrt(sigma2) /arma::accu(arma::abs(Phi)), 2));
      } else if ((method == "ALASSO") || (method == "MonoALASSO")) {
        for (gg = 0; gg < q; gg++) {
          lambda2(gg) = pow((sqrt(sigma2) / abs(Phi(gg))), 2);
        }
      }
    
    
    
  }
  
  
  //Rcpp::Rcout << 4 << std::endl;
  
 arma::mat DHt(T, T);
  DHt.zeros();
 arma::mat tHtDHt(T, 1);
  tHtDHt.zeros();
  
  int rr = 0;
  
  int TotalSim = nsim * by + burnin;
  
  //outputseq <- seq(burnin + 1, TotalSim, step)
  
 arma::mat V_(T, 1);
  V_.zeros();
  
 arma::mat V(T - q, 1);
  V.zeros();
  
 arma::mat Vas_(T, q);
  Vas_.zeros();
 arma::mat Vas(T - q, q);
  Vas.zeros();
  
  
 arma::mat VasPhi(T - q, 1);
 arma::mat resi(T - q, 1);
  
  
 arma::mat PhiMat(T - q, T);
 arma::mat C;
 arma::mat D;
  
  Rcpp::List TauGamma; 
  Rcpp::List ZetaBeta; 
  Rcpp::List mu0Mu; 
  
  arma::mat tmpY;
  
 arma::mat tmpSumTau; 
  
  //Rcpp::Rcout << 1 << std::endl;
  
  for (ii = 0; ii < TotalSim; ii++) {
    
    //if (ii % 100 == 0) {
   //  Rcpp::Rcout <<"Training: " << ((ii + 0.0) / (TotalSim + 0.0) * 100.0) << '%' << std::endl;
    //}
    
    //Rcpp::Rcout << "ii:" << ii << std::endl;
    
    //update V
    V_ = Y - Mu;
    V = V_.rows(q, T - 1);
    Vas_ = getV(V_, q);
    Vas = Vas_.rows(q, T - 1);
    
    //Rcpp::Rcout << Mu << std::endl;
    //Rcpp::Rcout << V << std::endl;
    
   // update Phi
    Phi = updatePhi(V, Vas, A, 
                    Phi, sigma2, inveta2mat, 
                    bound0, boundqplus1,
                    MonoFlg, method);
    
    //Rcpp::Rcout << Phi << std::endl;
    
   // Get residuals
    VasPhi = Vas * Phi;
    resi = V - VasPhi;
    
   // update sigma2
    sigma2 = updateSigma2(resi, Phi, inveta2mat, T, q, 
                          A, a, b, method);
    
    //Rcpp::Rcout << sigma2 << std::endl;
    
   // update eta2
    inveta2 = updateinveta2(Phi, sigma2, lambda2, q, tol);
    eta2 =arma::pow(inveta2, -1);
    inveta2mat.diag() = inveta2;
    
    //Rcpp::Rcout << inveta2 << std::endl;
    
   // update lambda2
    lambda2 = updatelambda2(eta2, q, alpha, beta, method);
    
    //Rcpp::Rcout << lambda2 << std::endl;
    
    ///////////////////////////////////////////////////
    //update the random level shift model
    ///////////////////////////////////////////////////
    
    //Calculate Phi Matrix 
    PhiMat = getPhiMat(Phi, T);
    
    //Calculate C Matrix
    C = G_ - PhiMat;
    
    //Calculate D Matrix
    D = C.t() * C;
    
    //#update Tau and Gamma
    
    //Rcpp::Rcout << 1 << std::endl;
    
    if (Hflg == 1) {
      tmpY = Y - mu0;
    
      if (Xflg == 1) {
        tmpY = tmpY - X_ * (Zeta % Beta);
      }
      
      TauGamma = updateTauGamma(tmpY, Phi, Tau, Gamma, 
                                mu0, sigma2, pho, xi2,
                                T, q, D, H_, Hflg, m);
      
      Tau = Rcpp::as<arma::mat>(TauGamma["Tau"]);
      Gamma = Rcpp::as<arma::mat>(TauGamma["Gamma"]);
    }
    
    //Rcpp::Rcout << 2 << std::endl;
    
    //Rcpp::Rcout << Tau << std::endl;
    //Rcpp::Rcout << Gamma << std::endl;
    
    //pvec = Rcpp::as<arma::mat>(TauGamma["p"]);
    //muGamma = Rcpp::as<arma::mat>(TauGamma["muGamma"]);
    //sigma2Gamma = Rcpp::as<arma::mat>(TauGamma["sigma2Gamma"]);
    
    //#update Zeta and Beta
    if (Xflg == 1) {
      tmpY = Y - mu0;
    
      if (Hflg == 1) {
        tmpY = tmpY - H_ * (Tau % Gamma);
      }
    
      ZetaBeta = updateZetaBeta(tmpY, Phi, Zeta, Beta, 
                              mu0, sigma2, phobeta, xi2,
                              T, q, D, X_, Xflg, betap);
    
      Zeta = Rcpp::as<arma::mat>(ZetaBeta["Zeta"]);
      Beta = Rcpp::as<arma::mat>(ZetaBeta["Beta"]);
    }
    
    //Rcpp::Rcout << 3 << std::endl;
    
    //Rcpp::Rcout << Zeta << std::endl;
    //Rcpp::Rcout << Beta << std::endl;
    
    //#update mu0 and Mu
    
    mu0Mu = updatemu0MuX(Y, Zeta, Beta, 
                        Tau, Gamma, 
                        sigma2, One, D,
                        H_, Hflg, X_, Xflg, 
                        T, tol);
    
    mu0 = mu0Mu["mu0"];
    //Rcpp::Rcout << mu0 << std::endl;
    
    Mu = Rcpp::as<arma::mat>(mu0Mu["Mu"]);
    
    //Rcpp::Rcout << 4 << std::endl;
    
    //Rcpp::Rcout << mu0 << std::endl;
    //Rcpp::Rcout << Mu << std::endl;
    
    //#update pho
    if (Hflg == 1) {
      tmpSumTau =arma::sum(Tau);
      pho = R::rbeta(theta1 + tmpSumTau(0), theta2 + m - tmpSumTau(0));
    }
    
    if (Xflg == 1) {
      tmpSumTau =arma::sum(Zeta);
      phobeta = R::rbeta(theta1 + tmpSumTau(0), theta2 + betap - tmpSumTau(0));
    }
    
    //Rcpp::Rcout << pho << std::endl;
    //Rcpp::Rcout << phobeta << std::endl;
    
    if (ii >= burnin) {
      if (ii % by == 0) {
        Phiout.col(rr) = Phi;
        sigma2out.col(rr) = sigma2;
        
        if (Hflg == 1) {
          Tauout.col(rr) = Tau;
          Gammaout.col(rr) = Gamma;
          phoout.col(rr) = pho;
        }
        if (Xflg == 1) {
          Zetaout.col(rr) = Zeta;
          Betaout.col(rr) = Beta;
          phobetaout.col(rr) = phobeta;
        }
        
        //muGammaout.col(rr) = muGamma;
        //sigma2Gammaout.col(rr) = sigma2Gamma;
        //pout.col(rr) = pvec;
        mu0out.col(rr) = mu0;
        Muout.col(rr) = Mu;
        
        
        if ((method == "LASSO") || (method == "ALASSO")) {
          eta2out.col(rr) = eta2;
          lambda2out.col(rr) = lambda2;
        }
        
        rr = rr + 1;
      }
    }
    
  }
  
  //Rcpp::Rcout << 1 << std::endl;
  
  /////////////////////////////////
  
  //Rcpp::Rcout <<"Training: 100%" << std::endl;
  //
  //auto end = std::chrono::system_clock::now();
  //std::chrono::duration<double> elapsed_seconds = end-start;
  //std::time_t end_time = std::chrono::system_clock::to_time_t(end);
  //
  //Rcpp::Rcout << "Finished training at " << std::ctime(&end_time)
 //            << "Elapsed time: " << elapsed_seconds.count() << "s"
 //            << std::endl;
  
  /////////////////////////////////
  
  Rcpp::List out = Rcpp::List::create(
    _["Phi"] = Phiout,
    _["sigma2"] = sigma2out,
    _["Tau"] = Tauout,
    _["Gamma"] = Gammaout,
    _["Zeta"] = Zetaout,
    _["Beta"] = Betaout,
    //_["muGamma"] = muGammaout,
    //_["sigma2Gamma"] = sigma2Gammaout,
    //_["p"] = pout,
    _["mu0"] = mu0out,
    _["Mu"] = Muout,
    _["pho"] = phoout,
    _["phobeta"] = phobetaout,
    _["eta2"] = eta2out,
    _["lambda2"] = lambda2out
  );
  
  //Rcpp::Rcout << 1 << std::endl;
  
  return(out);
  
}


// [[Rcpp::export]]
Rcpp::List GibbsRFLSMXYJZcpp(arma::colvec& Y,int& q, 
                              arma::mat& A, double& a, double& b, double& alpha, double& beta, 
                               double& theta1, double& theta2, double& xi2,
                               Rcpp::String& method, int monophi, double& bound0, double& boundqplus1,
                               int updateYJ, double& theta,
                               int leftcensoring, int rounding, double eps,
                               arma::colvec missingdata, double missingdatalb, double missingdataub,
                               int& nsim, int& by, int& burnin,
                               double& tol, 
                               Rcpp::Nullable<Rcpp::NumericMatrix> G = R_NilValue,
                               Rcpp::Nullable<Rcpp::List> oldpars = R_NilValue,
                               Rcpp::Nullable<Rcpp::NumericMatrix> X = R_NilValue,
                               Rcpp::Nullable<Rcpp::NumericMatrix> H = R_NilValue) {
  
  
  auto start = std::chrono::system_clock::now();
  std::time_t start_time = std::chrono::system_clock::to_time_t(start);
  
  if (monophi == 1) {
    Rcpp::Rcout << "Start training using " << method.get_cstring() << " with monotonicity constraints" << " at " << std::ctime(&start_time) <<  std::endl;
  
  } else {
    Rcpp::Rcout << "Start training using " << method.get_cstring() << " without monotonicity constraints" << " at " << std::ctime(&start_time) <<  std::endl;
  
  }
  
  
  ///////////////////////////////////
  
   arma::mat X_;
  
 // Calculate X
  int betap = 0;
  int Xflg = 0;
  if (X.isNotNull()) {
    X_ = Rcpp::as<arma::mat>(X);
    betap = X_.n_cols;
    Xflg = 1;
  } 
  
 arma::mat H_;
  
 // Calculate H
  int m = 0;
  int Hflg = 0;
  if (H.isNotNull()) {
    H_ = Rcpp::as<arma::mat>(H);
    m = H_.n_cols;
    Hflg = 1;
  } 
  
  
  int T = Y.n_elem;
  
  int TotalSim = nsim * by + burnin;
  
  Rcpp::List GibbsRFLSMModel; 
 arma::colvec Yyj = Y;
  
  Rcpp::List oldpars_;
  
 arma::mat Phi(q, 1);
  Phi.zeros();
  
 arma::mat Mu(T - q, 1);
  Mu.zeros();
  
  double sigma2 = 1.0;
  double theta_ = theta;
  //double theta0 = theta;
  //double kappa = 0.0;
 arma::mat tmp; 
  
  //Rcpp::Rcout << 1 << std::endl;
  
  if (updateYJ == 1) {
    Yyj = yeojohnsontr(Y, theta_, eps);
  }
  
  //Rcpp::Rcout << 2 << std::endl;
  
  if (oldpars.isNotNull()) {
    oldpars_ = Rcpp::as<Rcpp::List>(oldpars);
    Phi = Rcpp::as<arma::mat>(oldpars_["Phi"]);
    Mu = Rcpp::as<arma::mat>(oldpars_["Mu"]);
    sigma2 = oldpars_["sigma2"];
  } else {
    
    //Rcpp::Rcout << 3 << std::endl;
    
    oldpars_ = GibbsRFLSMXUpdatecpp(Yyj, q, 
                                   A, a, b, alpha, beta, 
                                   theta1, theta2, xi2, 
                                   method, monophi, bound0, boundqplus1, 
                                   1, 1, 0, tol,
                                   G, oldpars, X, H);
    Phi = Rcpp::as<arma::mat>(oldpars_["Phi"]);
    Mu = Rcpp::as<arma::mat>(oldpars_["Mu"]);
    sigma2 = oldpars_["sigma2"];
  }
  
  //Rcpp::Rcout << 21 << std::endl;
  
  //if (updateZ == 1) {
  //  Yyj = getucY(Yyj, Y, Phi, Mu, sigma2, theta_, eps);
  //}
  
  Rcpp::NumericMatrix GG;
 arma::mat G_;
  
  if (G.isNotNull()) {
    GG = Rcpp::as<Rcpp::NumericMatrix>(G);
  } else {
    G_ = getGMat(T, q);
    GG = Rcpp::as<Rcpp::NumericMatrix>(Rcpp::wrap(G_));
  }
  
  //Rcpp::Rcout << 2 << std::endl;
  
  arma::mat Z(T, 1);
  Z.zeros();
  
  //////////////////////////
  
 arma::mat Phiout(q, nsim); 
 arma::mat sigma2out(1, nsim);
 
 arma::mat Tauout;
 arma::mat Gammaout;
 arma::mat phoout;
 if (Hflg == 1) {
   Tauout.set_size(m, nsim);
   Gammaout.set_size(m, nsim);
   phoout.set_size(1, nsim);
 }
 
 arma::mat Zetaout;
 arma::mat Betaout;
 arma::mat phobetaout;
 if (Xflg == 1) {
   Zetaout.set_size(betap, nsim);
   Betaout.set_size(betap, nsim);
   phobetaout.set_size(1, nsim);
 }

 arma::mat mu0out(1, nsim);
 arma::mat Muout(T, nsim);
 
 
 arma::mat eta2out;
 arma::mat lambda2out;
 if ((method == "LASSO") || (method == "ALASSO")) {
   eta2out.set_size(q, nsim);
   lambda2out.set_size(q, nsim);
 }
 
 arma::mat thetaout;
 //arma::mat kappaout;
 if (updateYJ == 1) {
   thetaout.set_size(1, nsim);
   //kappaout.set_size(1, nsim);
 }
 
 Rcpp::List tmpList;
 
 arma::mat Zout;
 if ((leftcensoring == 1) || (rounding == 1) || (arma::accu(missingdata) > 0)) {
   Zout.set_size(T, nsim);
 }
  
  //////////////////////////
  
  
  //////////////////////////
  
  int i;
  int rr = 0;
  
  for (i = 0; i < TotalSim; i++) {
    
    //Rcpp::Rcout << "i:" << i << std::endl;
    
    if (i % 100 == 0) {
      Rcpp::Rcout <<"Training: " << ((i + 0.0) / (TotalSim + 0.0) * 100.0) << '%' << std::endl;
    }
    
    //if (updateZ == 1) {
    if ((leftcensoring == 1) || (rounding == 1) || (arma::accu(missingdata) > 0)) {
      Z = updateZZ(Y, Z, Phi, Mu, sigma2, 
                  theta_, leftcensoring, rounding, 0, 1, tol);
      
      Z = updateZZMD(Y, Z, Phi, Mu, sigma2, 
                  theta_, missingdata, missingdatalb, missingdataub, 0, 1, tol);
      
      Yyj = Y + Z;
    } else {
      Yyj = Y;
    }
    
    //Rcpp::Rcout << 2 << std::endl;
    
    if (updateYJ == 1) {
      //tmp = thetaYeoJohnsonMH(Y, Phi, Mu, sigma2, 
      //                        theta_, 1, 1, tol);
      tmp = updatethetaYJMH(Yyj, Phi, Mu, sigma2, 
                              theta_, 0, 1, tol);
      
      //tmpList = updateKappaTheta(Yyj, Phi, Mu, sigma2, 
      //                        kappa, theta0, 0, 1, tol);
      
      //Rcpp::Rcout << 3 << std::endl;
      
      //kappa = tmpList["Kappa"];
      //theta0 = tmpList["Theta"];
      //theta_ = pow(theta0, kappa);
      theta_ = tmp(0);
      //Rcpp::Rcout << theta_ << std::endl;
      Yyj = yeojohnsontr(Yyj, theta_, eps);
    } 
    
    //Rcpp::Rcout << 4 << std::endl;
   // 
   
   //Rcpp::Rcout << Yyj << std::endl;
   // Rcpp::Rcout << q << std::endl;
   // Rcpp::Rcout << A << std::endl;
   // Rcpp::Rcout << a << std::endl;
   // Rcpp::Rcout << b << std::endl;
   //Rcpp::Rcout << alpha << std::endl;
   //Rcpp::Rcout << beta << std::endl;
   //Rcpp::Rcout << theta1 << std::endl;
   //Rcpp::Rcout << theta2 << std::endl;
   //Rcpp::Rcout << xi2 << std::endl;
   //Rcpp::Rcout << monophi << std::endl;
   //Rcpp::Rcout << bound0 << std::endl;
   //Rcpp::Rcout << boundqplus1 << std::endl;
   //Rcpp::Rcout << tol << std::endl;
   
    //Phi = Rcpp::as<arma::mat>(oldpars_["Phi"]);
    //Mu = Rcpp::as<arma::mat>(oldpars_["Mu"]);
    //sigma2 = oldpars_["sigma2"];
    
    //Rcpp::Rcout << "Phi:" << Rcpp::as<arma::mat>(oldpars_["Phi"]) << std::endl;
    //Rcpp::Rcout << "Mu:" << Rcpp::as<arma::mat>(oldpars_["Mu"] )<< std::endl;
    //Rcpp::Rcout << "sigma2:" << Rcpp::as<double>(oldpars_["sigma2"]) << std::endl;
    //Rcpp::Rcout << "Tau:" << Rcpp::as<arma::mat>(oldpars_["Tau"]) << std::endl;
    //Rcpp::Rcout << "Gamma:" << Rcpp::as<arma::mat>(oldpars_["Gamma"]) << std::endl;
    //Rcpp::Rcout << "Zeta:" << Rcpp::as<arma::mat>(oldpars_["Zeta"]) << std::endl;
    //Rcpp::Rcout << "Beta:" << Rcpp::as<arma::mat>(oldpars_["Beta"]) << std::endl;
    //Rcpp::Rcout << "mu0:" << Rcpp::as<double>(oldpars_["mu0"]) << std::endl;
    //Rcpp::Rcout << "Mu:" << Rcpp::as<arma::mat>(oldpars_["Mu"]) << std::endl;
    //Rcpp::Rcout << "pho:" << Rcpp::as<double>(oldpars_["pho"]) << std::endl;
    //Rcpp::Rcout << "phobeta:" << Rcpp::as<double>(oldpars_["phobeta"]) << std::endl;
    //Rcpp::Rcout << "eta2:" << Rcpp::as<arma::mat>(oldpars_["eta2"]) << std::endl;
    //Rcpp::Rcout << "lambda2:" << Rcpp::as<arma::mat>(oldpars_["lambda2"]) << std::endl;
   
    oldpars_ = GibbsRFLSMXUpdatecpp(Yyj, q, 
                                   A, a, b, alpha, beta, 
                                   theta1, theta2, xi2, 
                                   method, monophi, bound0, boundqplus1, 
                                   1, 1, 0, tol,
                                   GG, oldpars_, X, H);
    
    Phi = Rcpp::as<arma::mat>(oldpars_["Phi"]);
    Mu = Rcpp::as<arma::mat>(oldpars_["Mu"]);
    sigma2 = oldpars_["sigma2"];
   //  
    // Rcpp::Rcout << 5 << std::endl;
    
    
    
    //if (updateZ == 1) {
    //  Yyj = getucY(Yyj, Y, Phi, Mu, sigma2, theta_, eps);
    //}
    
    if (i >= burnin) {
      if (i % by == 0) {
        Phiout.col(rr) = Rcpp::as<arma::mat>(oldpars_["Phi"]);
        sigma2out.col(rr) = Rcpp::as<arma::mat>(oldpars_["sigma2"]);
        
        if (Hflg == 1) {
          Tauout.col(rr) = Rcpp::as<arma::mat>(oldpars_["Tau"]);
          Gammaout.col(rr) = Rcpp::as<arma::mat>(oldpars_["Gamma"]);
        }
        
        if (Xflg == 1) {
          Zetaout.col(rr) = Rcpp::as<arma::mat>(oldpars_["Zeta"]);
          Betaout.col(rr) = Rcpp::as<arma::mat>(oldpars_["Beta"]);
        }
        
        mu0out.col(rr) = Rcpp::as<arma::mat>(oldpars_["mu0"]);
        Muout.col(rr) = Rcpp::as<arma::mat>(oldpars_["Mu"]);
        
        
        
        if ((method == "LASSO") || (method == "ALASSO")) {
          eta2out.col(rr) = Rcpp::as<arma::mat>(oldpars_["eta2"]);
          lambda2out.col(rr) = Rcpp::as<arma::mat>(oldpars_["lambda2"]);
        }
        
        if (updateYJ == 1) {
          thetaout.col(rr) = theta_;
          //kappaout.col(rr) = kappa;
        }
        //Yyjout.col(rr) = Yyj;
        
        if ((leftcensoring == 1) || (rounding == 1)) {
          Zout.col(rr) = Z;
        }
        rr = rr + 1;
      }
    }
    
  }
  
  Rcpp::Rcout <<"Training: 100%" << std::endl;
  
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  std::time_t end_time = std::chrono::system_clock::to_time_t(end);
  
  Rcpp::Rcout << "Finished training at " << std::ctime(&end_time)
              << "Elapsed time: " << elapsed_seconds.count() << "s"
              << std::endl;
  
  
  Rcpp::List out; 
  out = Rcpp::List::create(
    _["Phi"] = Phiout,
    _["sigma2"] = sigma2out,
    _["Tau"] = Tauout,
    _["Gamma"] = Gammaout,
    _["Zeta"] = Zetaout,
    _["Beta"] = Betaout,
    _["mu0"] = mu0out,
    _["Mu"] = Muout,
    //_["pho"] = phoout,
    _["eta2"] = eta2out,
    _["lambda2"] = lambda2out,
    _["theta"] = thetaout,
    //_["kappa"] = kappaout,
    _["Z"] = Zout
    //_["Yyj"] = Yyjout
  );
  
  return(out);
  
}


// [[Rcpp::export]]
Rcpp::List GibbsRFLSMXYJZcpp1(arma::colvec& Y,int& q, 
                              arma::mat& A, double& a, double& b, double& alpha, double& beta, 
                               double& theta1, double& theta2, double& xi2,
                               Rcpp::String& method, int monophi, double& bound0, double& boundqplus1,
                               int updateYJ, double& theta,
                               int leftcensoring, int rounding, double eps,
                               int& nsim, int& by, int& burnin,
                               double& tol, 
                               Rcpp::Nullable<Rcpp::NumericMatrix> G = R_NilValue,
                               Rcpp::Nullable<Rcpp::List> oldpars = R_NilValue,
                               Rcpp::Nullable<Rcpp::NumericMatrix> X = R_NilValue,
                               Rcpp::Nullable<Rcpp::NumericMatrix> H = R_NilValue) {
  
  
  auto start = std::chrono::system_clock::now();
  std::time_t start_time = std::chrono::system_clock::to_time_t(start);
  
  Rcpp::Rcout << "Start training using " << method.get_cstring() << " at " << std::ctime(&start_time) <<  std::endl;
  
  ///////////////////////////////////
  
   arma::mat X_;
  
 // Calculate X
  int betap = 0;
  int Xflg = 0;
  if (X.isNotNull()) {
    X_ = Rcpp::as<arma::mat>(X);
    betap = X_.n_cols;
    Xflg = 1;
  } 
  
 arma::mat H_;
  
 // Calculate H
  int m = 0;
  int Hflg = 0;
  if (H.isNotNull()) {
    H_ = Rcpp::as<arma::mat>(H);
    m = H_.n_cols;
    Hflg = 1;
  } 
  
  
  int T = Y.n_elem;
  
  int TotalSim = nsim * by + burnin;
  
  Rcpp::List GibbsRFLSMModel; 
 arma::colvec Yyj = Y;
  
  Rcpp::List oldpars_;
  
 arma::mat Phi(q, 1);
  Phi.zeros();
  
 arma::mat Mu(T - q, 1);
  Mu.zeros();
  
  double sigma2 = 1.0;
  double theta_ = theta;
 arma::mat tmp; 
  
  //Rcpp::Rcout << 1 << std::endl;
  
  if (updateYJ == 1) {
    Yyj = yeojohnsontr(Y, theta_, eps);
  }
  
  //Rcpp::Rcout << 2 << std::endl;
  
  if (oldpars.isNotNull()) {
    oldpars_ = Rcpp::as<Rcpp::List>(oldpars);
    Phi = Rcpp::as<arma::mat>(oldpars_["Phi"]);
    Mu = Rcpp::as<arma::mat>(oldpars_["Mu"]);
    sigma2 = oldpars_["sigma2"];
  } else {
    
    //Rcpp::Rcout << 3 << std::endl;
    
    oldpars_ = GibbsRFLSMXUpdatecpp(Yyj, q, 
                                   A, a, b, alpha, beta, 
                                   theta1, theta2, xi2, 
                                   method, monophi, bound0, boundqplus1, 
                                   1, 1, 0, tol,
                                   G, oldpars, X, H);
    Phi = Rcpp::as<arma::mat>(oldpars_["Phi"]);
    Mu = Rcpp::as<arma::mat>(oldpars_["Mu"]);
    sigma2 = oldpars_["sigma2"];
  }
  
  //Rcpp::Rcout << 21 << std::endl;
  
  //if (updateZ == 1) {
  //  Yyj = getucY(Yyj, Y, Phi, Mu, sigma2, theta_, eps);
  //}
  
  Rcpp::NumericMatrix GG;
 arma::mat G_;
  
  if (G.isNotNull()) {
    GG = Rcpp::as<Rcpp::NumericMatrix>(G);
  } else {
    G_ = getGMat(T, q);
    GG = Rcpp::as<Rcpp::NumericMatrix>(Rcpp::wrap(G_));
  }
  
  //Rcpp::Rcout << 2 << std::endl;
  
  arma::mat Z(T, 1);
  Z.zeros();
  
  //////////////////////////
  
 arma::mat Phiout(q, nsim); 
 arma::mat sigma2out(1, nsim);
 
 arma::mat Tauout;
 arma::mat Gammaout;
 arma::mat phoout;
 if (Hflg == 1) {
   Tauout.set_size(m, nsim);
   Gammaout.set_size(m, nsim);
   phoout.set_size(1, nsim);
 }
 
 arma::mat Zetaout;
 arma::mat Betaout;
 arma::mat phobetaout;
 if (Xflg == 1) {
   Zetaout.set_size(betap, nsim);
   Betaout.set_size(betap, nsim);
   phobetaout.set_size(1, nsim);
 }

 arma::mat mu0out(1, nsim);
 arma::mat Muout(T, nsim);
 
 
 arma::mat eta2out;
 arma::mat lambda2out;
 if ((method == "LASSO") || (method == "ALASSO")) {
   eta2out.set_size(q, nsim);
   lambda2out.set_size(q, nsim);
 }
 
 arma::mat thetaout;
 if (updateYJ == 1) {
   thetaout.set_size(1, nsim);
 }
 
 arma::mat Zout;
 if ((leftcensoring == 1) || (rounding == 1)) {
   Zout.set_size(T, nsim);
 }
  
  //////////////////////////
  
  
  //////////////////////////
  
  int i;
  int rr = 0;
  
  for (i = 0; i < TotalSim; i++) {
    
    //Rcpp::Rcout << "i:" << i << std::endl;
    
    if (i % 100 == 0) {
      Rcpp::Rcout <<"Training: " << ((i + 0.0) / (TotalSim + 0.0) * 100.0) << '%' << std::endl;
    }
    
    //if (updateZ == 1) {
    //if ((leftcensoring == 1) || (rounding == 1)) {
    //  Z = updateZZ(Y, Z, Phi, Mu, sigma2, 
    //              theta_, leftcensoring, rounding, 0, 1, tol);
    //  Yyj = Y + Z;
    //} else {
    //  Yyj = Y;
    //}
    
    if ((leftcensoring == 1) || (rounding == 1)) {
      Z = updateZ0(Y, Z, Phi, Mu, sigma2, 
                  theta_, leftcensoring, rounding, 0, 1, tol);
      Z = updateZSim(Y, Z, Phi, Mu, sigma2, 
                  theta_, tol, leftcensoring, rounding);
      Yyj = Y + Z;
    } else {
      Yyj = Y;
    }
    
    //Rcpp::Rcout << 2 << std::endl;
    
    if (updateYJ == 1) {
      //tmp = thetaYeoJohnsonMH(Y, Phi, Mu, sigma2, 
      //                        theta_, 1, 1, tol);
      tmp = updatethetaYJMH(Yyj, Phi, Mu, sigma2, 
                              theta_, 0, 1, tol);
      theta_ = tmp(0);
      //Rcpp::Rcout << theta_ << std::endl;
      Yyj = yeojohnsontr(Yyj, theta_, eps);
    } 
    
    //Rcpp::Rcout << 4 << std::endl;
   // 
   
   //Rcpp::Rcout << Yyj << std::endl;
   // Rcpp::Rcout << q << std::endl;
   // Rcpp::Rcout << A << std::endl;
   // Rcpp::Rcout << a << std::endl;
   // Rcpp::Rcout << b << std::endl;
   //Rcpp::Rcout << alpha << std::endl;
   //Rcpp::Rcout << beta << std::endl;
   //Rcpp::Rcout << theta1 << std::endl;
   //Rcpp::Rcout << theta2 << std::endl;
   //Rcpp::Rcout << xi2 << std::endl;
   //Rcpp::Rcout << monophi << std::endl;
   //Rcpp::Rcout << bound0 << std::endl;
   //Rcpp::Rcout << boundqplus1 << std::endl;
   //Rcpp::Rcout << tol << std::endl;
   
    //Phi = Rcpp::as<arma::mat>(oldpars_["Phi"]);
    //Mu = Rcpp::as<arma::mat>(oldpars_["Mu"]);
    //sigma2 = oldpars_["sigma2"];
    
    //Rcpp::Rcout << "Phi:" << Rcpp::as<arma::mat>(oldpars_["Phi"]) << std::endl;
    //Rcpp::Rcout << "Mu:" << Rcpp::as<arma::mat>(oldpars_["Mu"] )<< std::endl;
    //Rcpp::Rcout << "sigma2:" << Rcpp::as<double>(oldpars_["sigma2"]) << std::endl;
    //Rcpp::Rcout << "Tau:" << Rcpp::as<arma::mat>(oldpars_["Tau"]) << std::endl;
    //Rcpp::Rcout << "Gamma:" << Rcpp::as<arma::mat>(oldpars_["Gamma"]) << std::endl;
    //Rcpp::Rcout << "Zeta:" << Rcpp::as<arma::mat>(oldpars_["Zeta"]) << std::endl;
    //Rcpp::Rcout << "Beta:" << Rcpp::as<arma::mat>(oldpars_["Beta"]) << std::endl;
    //Rcpp::Rcout << "mu0:" << Rcpp::as<double>(oldpars_["mu0"]) << std::endl;
    //Rcpp::Rcout << "Mu:" << Rcpp::as<arma::mat>(oldpars_["Mu"]) << std::endl;
    //Rcpp::Rcout << "pho:" << Rcpp::as<double>(oldpars_["pho"]) << std::endl;
    //Rcpp::Rcout << "phobeta:" << Rcpp::as<double>(oldpars_["phobeta"]) << std::endl;
    //Rcpp::Rcout << "eta2:" << Rcpp::as<arma::mat>(oldpars_["eta2"]) << std::endl;
    //Rcpp::Rcout << "lambda2:" << Rcpp::as<arma::mat>(oldpars_["lambda2"]) << std::endl;
   
    oldpars_ = GibbsRFLSMXUpdatecpp(Yyj, q, 
                                   A, a, b, alpha, beta, 
                                   theta1, theta2, xi2, 
                                   method, monophi, bound0, boundqplus1, 
                                   1, 1, 0, tol,
                                   GG, oldpars_, X, H);
    
    Phi = Rcpp::as<arma::mat>(oldpars_["Phi"]);
    Mu = Rcpp::as<arma::mat>(oldpars_["Mu"]);
    sigma2 = oldpars_["sigma2"];
   //  
    // Rcpp::Rcout << 5 << std::endl;
    
    
    
    //if (updateZ == 1) {
    //  Yyj = getucY(Yyj, Y, Phi, Mu, sigma2, theta_, eps);
    //}
    
    if (i >= burnin) {
      if (i % by == 0) {
        Phiout.col(rr) = Rcpp::as<arma::mat>(oldpars_["Phi"]);
        sigma2out.col(rr) = Rcpp::as<arma::mat>(oldpars_["sigma2"]);
        
        if (Hflg == 1) {
          Tauout.col(rr) = Rcpp::as<arma::mat>(oldpars_["Tau"]);
          Gammaout.col(rr) = Rcpp::as<arma::mat>(oldpars_["Gamma"]);
        }
        
        if (Xflg == 1) {
          Zetaout.col(rr) = Rcpp::as<arma::mat>(oldpars_["Zeta"]);
          Betaout.col(rr) = Rcpp::as<arma::mat>(oldpars_["Beta"]);
        }
        
        mu0out.col(rr) = Rcpp::as<arma::mat>(oldpars_["mu0"]);
        Muout.col(rr) = Rcpp::as<arma::mat>(oldpars_["Mu"]);
        
        
        
        if ((method == "LASSO") || (method == "ALASSO")) {
          eta2out.col(rr) = Rcpp::as<arma::mat>(oldpars_["eta2"]);
          lambda2out.col(rr) = Rcpp::as<arma::mat>(oldpars_["lambda2"]);
        }
        
        if (updateYJ == 1) {
          thetaout.col(rr) = theta_;
        }
        //Yyjout.col(rr) = Yyj;
        
        if ((leftcensoring == 1) || (rounding == 1)) {
          Zout.col(rr) = Z;
        }
        rr = rr + 1;
      }
    }
    
  }
  
  Rcpp::Rcout <<"Training: 100%" << std::endl;
  
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  std::time_t end_time = std::chrono::system_clock::to_time_t(end);
  
  Rcpp::Rcout << "Finished training at " << std::ctime(&end_time)
              << "Elapsed time: " << elapsed_seconds.count() << "s"
              << std::endl;
  
  
  Rcpp::List out; 
  out = Rcpp::List::create(
    _["Phi"] = Phiout,
    _["sigma2"] = sigma2out,
    _["Tau"] = Tauout,
    _["Gamma"] = Gammaout,
    _["Zeta"] = Zetaout,
    _["Beta"] = Betaout,
    _["mu0"] = mu0out,
    _["Mu"] = Muout,
    //_["pho"] = phoout,
    _["eta2"] = eta2out,
    _["lambda2"] = lambda2out,
    _["theta"] = thetaout,
    _["Z"] = Zout
    //_["Yyj"] = Yyjout
  );
  
  return(out);
  
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
arma::colvec getfityj(arma::colvec Yyj,arma::mat Phi,arma::mat Mu, double eps) {
  
  int q = Phi.n_rows;
  //Rcpp::Rcout << q << std::endl;
  int T = Yyj.n_elem;
  //Rcpp::Rcout << T << std::endl;
 arma::mat V(T, 1); 
 arma::mat Vas(q, 1);
 arma::mat VasPhi;
  
 arma::mat fit(T - q, 1);
 arma::colvec tmp;
  
  //Rcpp::Rcout << 1 << std::endl;
  
  for (int i = 0; i < T; i++) {
    V(i, 0) = Yyj(i) - Mu(i, 0);
    
    //Rcpp::Rcout << i << std::endl;
    
    if (i >= q) {
      for (int j = 0; j < q; j++) {
        Vas(j, 0) = V(i - 1 - j, 0);
      }
      VasPhi = Vas.t() * Phi;
      fit(i - q, 0) = Mu(i, 0) + VasPhi(0);
    }
    //Rcpp::Rcout << 3 << std::endl;
    //Rcpp::Rcout << ucYyj(i)  << std::endl;
    
    
    //Rcpp::Rcout << 4 << std::endl;
  }
  
  //Rcpp::Rcout << ucYyj << std::endl;
  
  //ucY = invyeojohnsontr(ucYyj, theta, eps);
  
  //Rcpp::Rcout << 5 << std::endl;
  
  return(fit);
  
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
arma::colvec getfit(arma::colvec Yyj,arma::mat Phi,arma::mat Mu, double theta, double eps) {
  
 arma::colvec fittr = getfityj(Yyj, Phi, Mu, eps);
 arma::colvec fit = invyeojohnsontr(fittr, theta, eps);
  return(fit);
  
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
arma::colvec simYyjph1(arma::colvec Yyj,arma::mat Phi,arma::mat Mu, double sigma2) {
  
  int q = Phi.n_rows;
  //Rcpp::Rcout << q << std::endl;
  int T = Yyj.n_elem;
  //Rcpp::Rcout << T << std::endl;
 arma::mat V(T, 1); 
 arma::mat Vas(q, 1);
 arma::mat VasPhi;
  
 arma::mat fit(T - q, 1);
 arma::mat simYyjph1(T - q, 1);
 arma::colvec tmp;
  
  //Rcpp::Rcout << 1 << std::endl;
  
  for (int i = 0; i < T; i++) {
    V(i, 0) = Yyj(i) - Mu(i, 0);
    
    //Rcpp::Rcout << i << std::endl;
    
    if (i >= q) {
      for (int j = 0; j < q; j++) {
        Vas(j, 0) = V(i - 1 - j, 0);
      }
      VasPhi = Vas.t() * Phi;
      fit(i - q, 0) = Mu(i, 0) + VasPhi(0);
      simYyjph1(i - q, 0) = R::rnorm(fit(i - q, 0), sqrt(sigma2));
    }
    
    //Rcpp::Rcout << 3 << std::endl;
    //Rcpp::Rcout << ucYyj(i)  << std::endl;
    
    
    //Rcpp::Rcout << 4 << std::endl;
  }
  
  //Rcpp::Rcout << ucYyj << std::endl;
  
  //ucY = invyeojohnsontr(ucYyj, theta, eps);
  
  //Rcpp::Rcout << 5 << std::endl;
  
  return(simYyjph1);
  
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
arma::colvec simYph1(arma::colvec Yyj,arma::mat Phi,arma::mat Mu, double sigma2, double theta, double eps) {
  
 arma::colvec Yyjph1 = simYyjph1(Yyj, Phi, Mu, sigma2); 
 arma::colvec Yph1 = invyeojohnsontr(Yyjph1, theta, eps); 
 arma::uvec ind0 =arma::find(Yph1 <= 0.0); 
  Yph1(ind0).fill(0.0);
  
  
  return(Yph1);
  
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
arma::colvec simYyjph2(int h, arma::colvec Yyjph1, arma::colvec Phi, arma::colvec Mu, double sigma2) {
  
  int q = Phi.n_rows;
  //Rcpp::Rcout << q << std::endl;
  int T = Yyjph1.n_elem;
  //Rcpp::Rcout << T << std::endl;
 arma::mat V(T + h, 1); 
 arma::mat Vas(q, 1);
 arma::mat VasPhi;
  
 arma::mat fit(T - q + h, 1);
 arma::mat simYyjph2(T - q + h, 1);
 arma::colvec tmp;
  
  //Rcpp::Rcout << 1.1 << std::endl;
  
  //for (int i = 0; i < (T + h); i++) {
  for (int i = 0; i < (T + h); i++) {
    
    
    //Rcpp::Rcout << i << std::endl;
    
    if (i >= q) {
      for (int j = 0; j < q; j++) {
        Vas(j, 0) = V(i - 1 - j, 0);
      }
      VasPhi = Vas.t() * Phi;
      fit(i - q, 0) = Mu(i) + VasPhi(0);
      simYyjph2(i - q, 0) = R::rnorm(fit(i - q, 0), sqrt(sigma2));
    }
    
    //Rcpp::Rcout << 1.2 << std::endl;
    
    if (i < T) {
      V(i, 0) = Yyjph1(i) - Mu(i);
    } else if (i >= T) {
      V(i, 0) = simYyjph2(i - q, 0) - Mu(i);
    }
    
    //Rcpp::Rcout << 1.3 << std::endl;
    
    //Rcpp::Rcout << 3 << std::endl;
    //Rcpp::Rcout << ucYyj(i)  << std::endl;
    
    
    //Rcpp::Rcout << 4 << std::endl;
  }
  
  //Rcpp::Rcout << ucYyj << std::endl;
  
  //ucY = invyeojohnsontr(ucYyj, theta, eps);
  
  //Rcpp::Rcout << 5 << std::endl;
  
  return(simYyjph2);

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
arma::colvec simYph2(int h, arma::colvec Y, arma::colvec Z, arma::colvec Phi,arma::colvec Mu, double sigma2, 
                     int updateYJ, double theta, int leftcensoring, int rounding, double eps, int backtr) {
  
  arma::colvec Yyjph1 = Y;
  if ((leftcensoring == 1) || (rounding == 1)) {
    Yyjph1 = Yyjph1 + Z;
  }
  if (updateYJ == 1) {
    Yyjph1 = yeojohnsontr(Yyjph1, theta, eps);
  }
  
  //Rcpp::Rcout << 1 << std::endl;
  
  arma::colvec Yyjph2 = simYyjph2(h, Yyjph1, Phi, Mu, sigma2); 
  //Rcpp::Rcout << 2 << std::endl;
  
  arma::colvec Yph2 = Yyjph2;
  
  arma::uvec ind0;
  
  if (backtr == 1) { 
  
    if (updateYJ == 1) {
      Yph2 = invyeojohnsontr(Yph2, theta, eps); 
    //Rcpp::Rcout << 3 << std::endl;
    }
   
    if ((leftcensoring == 1)) {
      ind0 =arma::find(Yph2 <= 0.0); 
      Yph2(ind0).fill(0.0);
    }
    
    if ((rounding == 1)) {
      Yph2 = arma::round(Yph2);
    } 
  }
  
  return(Yph2);
  
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
arma::colvec simYyjph2NoY(int h, arma::colvec Yyjph1, arma::colvec Phi, arma::colvec Mu, double sigma2) {
  
  int q = Phi.n_rows;
  //Rcpp::Rcout << q << std::endl;
  int T = Yyjph1.n_elem;
  //Rcpp::Rcout << T << std::endl;
 arma::mat V(q + h, 1); 
 arma::mat Vas(q, 1);
 arma::mat VasPhi;
  
 arma::mat fit(h, 1);
 arma::mat simYyjph2(h, 1);
 arma::colvec tmp;
  
  //Rcpp::Rcout << 1.1 << std::endl;
  
  //for (int i = 0; i < (T + h); i++) {
  for (int i = 0; i < (q + h); i++) {
    
    
    //Rcpp::Rcout << i << std::endl;
    
    if (i >= q) {
      for (int j = 0; j < q; j++) {
        Vas(j, 0) = V(i - 1 - j, 0);
      }
      VasPhi = Vas.t() * Phi;
      fit(i - q, 0) = Mu(i) + VasPhi(0);
      simYyjph2(i - q, 0) = R::rnorm(fit(i - q, 0), sqrt(sigma2));
    }
    
    //Rcpp::Rcout << 1.2 << std::endl;
    
    if (i < q) {
      V(i, 0) = Yyjph1(i) - Mu(i);
    } else if (i >= q) {
      V(i, 0) = simYyjph2(i - q, 0) - Mu(i);
    }
    
    //Rcpp::Rcout << 1.3 << std::endl;
    
    //Rcpp::Rcout << 3 << std::endl;
    //Rcpp::Rcout << ucYyj(i)  << std::endl;
    
    
    //Rcpp::Rcout << 4 << std::endl;
  }
  
  //Rcpp::Rcout << ucYyj << std::endl;
  
  //ucY = invyeojohnsontr(ucYyj, theta, eps);
  
  //Rcpp::Rcout << 5 << std::endl;
  
  return(simYyjph2);

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
arma::colvec simYph2NoY(int h, arma::colvec Y, arma::colvec Z, arma::colvec Phi,arma::colvec Mu, double sigma2, 
                     int updateYJ, double theta, int leftcensoring, int rounding, double eps, int backtr) {
  
  arma::colvec Yyjph1 = Y;
  if ((leftcensoring == 1) || (rounding == 1)) {
    Yyjph1 = Yyjph1 + Z;
  }
  if (updateYJ == 1) {
    Yyjph1 = yeojohnsontr(Yyjph1, theta, eps);
  }
  
  //Rcpp::Rcout << 1 << std::endl;
  
  arma::colvec Yyjph2 = simYyjph2NoY(h, Yyjph1, Phi, Mu, sigma2); 
  //Rcpp::Rcout << 2 << std::endl;
  
  arma::colvec Yph2 = Yyjph2;
  
  arma::uvec ind0;
  
  if (backtr == 1) { 
  
    if (updateYJ == 1) {
      Yph2 = invyeojohnsontr(Yph2, theta, eps); 
    //Rcpp::Rcout << 3 << std::endl;
    }
   
    if ((leftcensoring == 1)) {
      ind0 =arma::find(Yph2 <= 0.0); 
      Yph2(ind0).fill(0.0);
    }
    
    if ((rounding == 1)) {
      Yph2 = arma::round(Yph2);
    } 
  }
  
  return(Yph2);
  
}
