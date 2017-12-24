#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::mat armapmax(arma::mat A, double bounds){
  int n = A.n_rows;
  int p = A.n_cols;
  arma::mat out = arma::mat(n,p);
  for (int i=0; i<n; ++i){
    for (int j=0; j<p; ++j){
      out(i,j) = max(A(i,j), bounds);
    }
  }
  return out;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::mat armapmin(arma::mat A, double bound){
  int n = A.n_rows;
  int p = A.n_cols;
  arma::mat out = arma::mat(n,p);
  for (int i=0; i<n; ++i){
    for (int j=0; j<p; ++j){
      out(i,j) = min(A(i,j), bound);
    }
  }
  return out;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::mat nonnegASC_c(arma::mat B){
  int n = B.n_rows;
  int m = B.n_cols;
  arma::mat A = arma::zeros<arma::mat>(n,m);
  A.each_row() = arma::linspace(1, m, m).t();
  arma::mat B_sort = sort(B,"descend",1);
  arma::mat cum_B = cumsum(B_sort, 1);
  arma::mat sigma = B_sort-(cum_B-1)/A;
  arma::uvec idx = arma::uvec(n);
  for (int i=0; i<n; ++i){
    arma::uvec tempuvec = find(sigma.row(i) > 0);
    idx(i) = tempuvec.size()-1;
  }
  arma::mat tmp = B_sort-sigma;
  sigma = tmp.cols(idx);
  arma::mat X = arma::mat(m,n);
  X.each_col() = sigma.diag();
  arma::mat out = armapmax(B-X, 0);
  return out;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::mat dist_c(arma::mat X){
  int n = X.n_rows;
  //int p = X.n_cols;
  arma::mat out = arma::zeros<arma::mat>(n,n);
  for (int i=0; i<(n-1); ++i){
    for (int j=(i+1); i<n; ++i){
      arma::rowvec diff = X.row(i) - X.row(j);
      out(i,j) = sqrt(accu(diff%diff));
      out(j,i) = out(i,j);
    }
  }
  return out;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::mat constructKernel_c(arma::mat X, double sigma){
  arma::mat dist = dist_c(X) / (2*sigma*sigma);
  return exp(-dist);
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
Rcpp::List RMSC_c(arma::cube T,
               double lambda,
               double mu=1e-3,
               double rho = 1.9,
               int max_iter=100,
               double eps = 1e-9,
               bool verbose = false){
  int m = T.n_rows; int p = T.n_cols; int n = T.n_slices;
  if(m!=p){
    Rcerr << "input matrix T must be square transition matrix\n";
  }
  arma::mat Z = arma::zeros<arma::mat>(m,p);
  arma::cube E = arma::zeros<arma::cube>(m,p,n);
  E.randn();
  arma::cube Y = arma::zeros<arma::cube>(m,p,n);
  arma::cube B = arma::zeros<arma::cube>(m,p,n);
  arma::mat Q = arma::zeros<arma::mat>(m,p);
  arma::mat P = arma::zeros<arma::mat>(m,p);
  arma::mat P_old = arma::zeros<arma::mat>(m,p);
  P_old.randn();
  //arma::vec e = arma::ones<arma::vec>(m);
  int step = 0;
  while(1){
    step += 1;
    double max_inf_norm = -1;
    for (int i=0; i < n; ++i){
      arma::mat Ti = T.slice(i); arma::mat Ei = E.slice(i);
      arma::mat diff = Ti-Ei-P;
      double inf_norm = norm(diff, "inf");
      max_inf_norm = max(max_inf_norm, inf_norm);
    }
    arma::vec s;
    svd(s, P);
    double funV = sum(s) + lambda * accu(abs(E));
    double relChg = norm(P-P_old, "fro")/max(1.0, norm(P_old, "fro"));
    P_old = P;
    arma::mat tmp = P-Q;
    double max_inf_norm2 = norm(tmp, "inf");
    if(verbose && step%10==0){
      cout << "iter" << step <<
        ": \n max_inf_norm = " << max_inf_norm <<
        "\n max_inf_norm2 = " << max_inf_norm2 <<
        "\n relChg = " << relChg <<
        "\n mu = " << mu <<
        "\n funV = " << funV << "\n";
    }
    if(step > 1 & max_inf_norm < eps){
      break;
    }
    if(step > max_iter){
      cout << "reached max iteration \n";
      break;
    }

    //Update P
    arma::mat temp = sum(T-E-Y/mu, 2);
    arma::mat B = (Q- Z/mu + temp)/(n+1);
    P = nonnegASC_c(B);
    for (int i=0; i < m; ++i){
      if(abs(sum(P.row(i)) - 1) > 1e-10){
        Rcerr << "rowsum 1 error";
      }
    }

    //Update Q;
    arma::mat M = P+Z/mu;
    double C = 1/mu;
    arma::mat U; arma::vec Sigma; arma::mat V;
    svd(U,Sigma,V,M);
    arma::uvec svp2 = find(Sigma > C);
    int svp = svp2.size();
    if(svp >=2){
      arma::vec Sigma2 = Sigma.subvec(0,svp-1) - C;
      Q = U.cols(0,svp-1) * diagmat(Sigma2) * V.cols(0,svp-1).t();
    }else if(svp==1){
      double Sigma2 = Sigma(0) - C;
      Q = U.col(0) * V.col(0).t() * Sigma2;
    }else{
      svp = 1;
      Q = arma::zeros<arma::mat>(m,m);
    }

    //Update Ei and Yi
    for (int i=0; i < n; ++i){
      arma::mat C = T.slice(i) - P - Y.slice(i) / mu;
      E.slice(i) = armapmax(C-lambda/mu, 0) + armapmin(C+lambda/mu,0);
      Y.slice(i) = Y.slice(i) + mu * (P+E.slice(i)-T.slice(i));
    }
    //Update Z
    Z = Z + mu*(P-Q);
    mu = min(rho*mu, 1e+10);
  }
  arma::mat U; arma::vec s; arma::mat V;
  svd(U,s,V,P.t());
  arma::vec pi = V.col(0);
  arma::vec Dist = pi / accu(pi);
  arma::vec invDist = 1/Dist;
  arma::mat Distmat = diagmat(Dist);
  arma::mat invDistmat = diagmat(invDist);
  arma::mat sspi = sqrt(invDistmat);
  arma::mat spi = sqrt(Distmat);
  P = (spi * P * sspi + sspi * P.t() * spi)/2;
  return Rcpp::List::create(
    Rcpp::Named("P") = P,
    Rcpp::Named("E") = E
  );
}


/* temporarily give up because of kmeans algorithm
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::vec baseline_spectral_c(arma::mat Phat, int numClust, arma::vec truth){
  arma::mat Utemp;
  arma::vec s;
  arma::mat V;
  svd(Utemp, s, V, Phat);
  arma::mat U = Utemp.rows(0,numClust-1);
  arma::vec Urowsum = sum(U%U,1);
  arma::mat norm_mat(U.n_rows, numClust);
  norm_mat.each_col() = Urowsum;
  arma::uvec errorind = find(U.col(0)==0);
  norm_mat.rows(errorind).fill(1);
  U /= norm_mat;


}

*/

