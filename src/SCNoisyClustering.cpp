#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

Environment quadprog("package:quadprog");
Function solve_QP = quadprog["solve.QP"];

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::mat armapmax(arma::mat A, double bound){
  int n = A.n_rows;
  int p = A.n_cols;
  arma::mat out = arma::mat(n,p);
  for (int i=0; i<n; ++i){
    for (int j=0; j<p; ++j){
      out(i,j) = max(A(i,j), bound);
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
  arma::mat Z  = arma::zeros<arma::mat>(m,p);
  arma::cube E = arma::zeros<arma::cube>(m,p,n);
  E.randn();
  arma::cube Y = arma::zeros<arma::cube>(m,p,n);
  arma::cube B = arma::zeros<arma::cube>(m,p,n);
  arma::mat Q  = arma::zeros<arma::mat>(m,p);
  arma::mat P  = arma::zeros<arma::mat>(m,p);
  arma::mat P_old = arma::zeros<arma::mat>(m,p);
  P_old.randn();
  int step = 0;
  while(1){
    step += 1;
    double max_inf_norm = -1;
    for (int i=0; i < n; ++i){
      arma::mat Ti    = T.slice(i); arma::mat Ei = E.slice(i);
      arma::mat diff  = Ti-Ei-P;
      double inf_norm = norm(diff, "inf");
      max_inf_norm    = max(max_inf_norm, inf_norm);
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
    Rcpp::Named("S") = P,
    Rcpp::Named("E") = E
  );
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
Rcpp::List objective_1_c(const arma::cube P,
                         const int numClust,
                         double mu=1e-3,
                         const double rho = 1.9,
                         const int max_iter=100,
                         const double eps = 1e-9,
                         const bool verbose = false){
  int m = P.n_rows; int p = P.n_cols; int n = P.n_slices;
  arma::vec  funV  = arma::zeros<arma::vec>(max_iter);
  arma::vec  fnorm = arma::zeros<arma::vec>(max_iter);
  arma::vec  w     = arma::ones<arma::vec>(n)* 1/n;
  arma::cube Q     = arma::zeros<arma::cube>(m,p,n);
  arma::cube Y     = arma::zeros<arma::cube>(m,p,n);
  arma::mat  S     = arma::zeros<arma::mat>(m,p);
  arma::mat  S_old = arma::zeros<arma::mat>(m,p);
  S_old.randn();
  int step = 0;
  while(1){
    step += 1;
    if(step > max_iter-1){
      cout << "reached max iteration \n";
      break;
    }
    double max_inf_norm = -1;
    for (int i=0; i < n; ++i){
      arma::mat Qi = Q.slice(i);
      arma::mat diff = S-Qi;
      double inf_norm = norm(diff, "inf");
      max_inf_norm = max(max_inf_norm, inf_norm);
    }
    funV(step-1) = sum(w%w);
    fnorm(step-1) = 0;
    for (int i=0; i < n; ++i){
      funV(step-1) += norm(Q.slice(i)-P.slice(i), "fro") * norm(Q.slice(i)-P.slice(i), "fro") * w(i);
      fnorm(step-1) += norm(Q.slice(i)-P.slice(i), "fro") * norm(Q.slice(i)-P.slice(i), "fro") * w(i);
    }
    double relChg = norm(S-S_old, "fro")/max(1.0, norm(S_old, "fro"));
    S_old = S;
    if(verbose){
      cout << "iter" << step <<
        ": \n max_inf_norm = " << max_inf_norm <<
          "\n relChg = " << relChg <<
            "\n mu = " << mu <<
              "\n funV = " << funV(step-1) << "\n";
    }
    if(step > 1 & max_inf_norm < eps & relChg < eps){
      cout << "converged!\n";
      break;
    }

    //Update Qi;
    for (int i = 0; i < n; ++i){
      arma::mat tmp = mu*S + P.slice(i) * w(i) + Y.slice(i);
      tmp = tmp / (w(i) + mu);
      arma::mat U; arma::vec d; arma::mat V;
      svd(U,d,V,tmp);
      d.subvec(numClust, d.size()-1) = arma::zeros<arma::vec>(d.size()-numClust);
      Q.slice(i) = U * diagmat(d) * V.t();
    }

    //update S;
    arma::mat temp = sum(Q - Y/mu, 2) / n;
    S = nonnegASC_c(temp);

    //update w
    arma::mat Dmat = arma::eye<arma::mat>(n,n);
    arma::vec dvec(n);
    for (int i=0; i<n; ++i){
      dvec(i) = norm(Q.slice(i)-P.slice(i), "fro");
      dvec(i) = -dvec(i) * dvec(i)/2;
    }
    arma::mat Amat = arma::mat(n, n+1);
    Amat.col(0) = arma::ones<arma::vec>(n);
    Amat.submat(0,1,n-1,n) = arma::eye<arma::mat>(n,n);
    arma::vec bvec = arma::zeros<arma::vec>(n+1);
    bvec(0) = 1;
    Rcpp::List qp = solve_QP(Dmat, dvec, Amat, bvec, 1);
    w = as<arma::vec>(qp["solution"]);

    //update Y
    for (int i = 0; i < n; ++i){
      Y.slice(i) = Y.slice(i) + mu * (S - Q.slice(i));
    }

    //update mu
    mu = min(rho*mu, 1e+10);
  }
  arma::mat U; arma::vec s; arma::mat V;
  svd(U,s,V,S);
  arma::vec pi = V.col(0);
  arma::vec Dist = pi / accu(pi) + 1e-10;
  arma::vec invDist = 1/Dist;
  arma::mat Distmat = diagmat(Dist);
  arma::mat invDistmat = diagmat(invDist);
  arma::mat sspi = sqrt(invDistmat);
  arma::mat spi = sqrt(Distmat);
  S = (spi * S * sspi + sspi * S.t() * spi)/2;
  return Rcpp::List::create(
    Rcpp::Named("S") = S,
    Rcpp::Named("f") = funV,
    Rcpp::Named("fnorm") = fnorm,
    Rcpp::Named("w") = w
  );
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
Rcpp::List objective_2_c(const arma::cube P,
                         const int numClust,
                         double mu=1e-3,
                         const double rho = 1.9,
                         const int max_iter=100,
                         const double eps = 1e-9,
                         const bool verbose = false){

  int m = P.n_rows; int p = P.n_cols; int n = P.n_slices;
  arma::vec funV = arma::zeros<arma::vec>(max_iter);
  arma::vec sigma = arma::ones<arma::vec>(n);
  arma::cube Q = arma::zeros<arma::cube>(m,p,n);
  arma::cube Y = arma::zeros<arma::cube>(m,p,n);
  arma::mat S = arma::zeros<arma::mat>(m,p);
  arma::mat S_old = arma::zeros<arma::mat>(m,p);
  S_old.randn();

  int step = 0;
  while(1){
    arma::vec mmsigma = m*m*sigma;
    step += 1;
    if(step > max_iter-1){
      cout << "reached max iteration \n";
      break;
    }
    double max_inf_norm = -1;
    for (int i=0; i < n; ++i){
      arma::mat Qi = Q.slice(i);
      arma::mat diff = S-Qi;
      double inf_norm = norm(diff, "inf");
      max_inf_norm = max(max_inf_norm, inf_norm);
    }
    funV(step-1) = accu(sigma/2);
    for (int i=0; i < n; ++i){
      funV(step-1) += norm(Q.slice(i)-P.slice(i), "fro")/(mmsigma(i));
    }
    double relChg = norm(S-S_old, "fro")/max(1.0, norm(S_old, "fro"));
    S_old = S;
    if(verbose){
      cout << "iter" << step <<
        ": \n max_inf_norm = " << max_inf_norm <<
          "\n relChg = " << relChg <<
            "\n mu = " << mu <<
              "\n funV = " << funV(step-1) << "\n";
    }
    if(step > 1 & max_inf_norm < eps & relChg < eps){
      cout << "converged!\n";
      break;
    }

    //Update Qi;
    for (int i = 0; i < n; ++i){
      arma::mat tmp = (mu*S + P.slice(i)/mmsigma(i) + Y.slice(i));
      tmp = tmp *(mmsigma(i) / (mu*mmsigma(i) + 1));
      arma::mat U; arma::vec d; arma::mat V;
      svd(U,d,V,tmp);
      d.subvec(numClust, d.size()-1) = arma::zeros<arma::vec>(d.size()-numClust);
      Q.slice(i) = U * diagmat(d) * V.t();
      //Q.slice(i) = Q.slice(i).each_col() / sum(Q.slice(i), 1);
    }

    //update S;
    arma::mat tmp = sum(Q - Y/mu, 2) / n;
    S = nonnegASC_c(tmp);

    //update sigma
    for (int i=0; i<n; ++i){
      double tmpnum = norm(Q.slice(i)-P.slice(i), "fro");
      sigma(i) = tmpnum*tmpnum/(m*m);
    }

    //update Y
    for (int i = 0; i < n; ++i){
      Y.slice(i) = Y.slice(i) + mu * (S - Q.slice(i));
    }
    //update mu
    mu = min(rho*mu, 1e+10);
  }
  arma::mat U; arma::vec s; arma::mat V;
  svd(U,s,V,S);
  arma::vec pi = V.col(0);
  arma::vec Dist = pi / accu(pi) + 1e-10;
  arma::vec invDist = 1/Dist;
  arma::mat Distmat = diagmat(Dist);
  arma::mat invDistmat = diagmat(invDist);
  arma::mat sspi = sqrt(invDistmat);
  arma::mat spi = sqrt(Distmat);
  S = (spi * S * sspi + sspi * S.t() * spi)/2;
  return Rcpp::List::create(
    Rcpp::Named("S") = S,
    Rcpp::Named("f") = funV,
    Rcpp::Named("sigma") = sigma
  );
}




// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
Rcpp::List objective_1_c_one_aux(const arma::cube P,
                              const int numClust,
                              double mu=1e-3,
                              const double rho = 1.9,
                              const int max_iter=100,
                              const double eps = 1e-9,
                              const bool verbose = false){
  int m = P.n_rows; int p = P.n_cols; int n = P.n_slices;
  arma::vec funV   = arma::zeros<arma::vec>(max_iter);
  arma::vec w      = arma::ones<arma::vec>(n)* 1/n;
  arma::mat Q      = arma::zeros<arma::mat>(m,p);
  arma::mat Y      = arma::zeros<arma::mat>(m,p);
  arma::mat S      = arma::zeros<arma::mat>(m,p);
  arma::mat S_old  = arma::zeros<arma::mat>(m,p);
  S_old.randn();
  S.randn();
  S = nonnegASC_c(S);
  int step = 0;
  while(1){
    step += 1;
    if(step > max_iter-1){
      cout << "reached max iteration \n";
      break;
    }
    double max_inf_norm = norm(S-Q, "inf");
    funV(step-1) = sum(w%w);
    for (int i=0; i < n; ++i){
      funV(step-1) += norm(Q-P.slice(i), "fro") * w(i);
    }
    double relChg = norm(S-S_old, "fro")/max(1.0, norm(S_old, "fro"));
    S_old = S;
    if(verbose){
      cout << "iter" << step <<
        ": \n max_inf_norm = " << max_inf_norm <<
          "\n relChg = " << relChg <<
            "\n mu = " << mu <<
              "\n funV = " << funV(step-1) << "\n";
    }
    if(step > 1 & max_inf_norm < eps & relChg < eps){
      cout << "converged!\n";
      break;
    }

    //Update Qi;
    arma::mat tmp = S*mu + Y;
    for (int i=0; i < n; ++i){
      tmp = tmp + P.slice(i) * w(i);
    }
    tmp = tmp / (mu + 1);
    arma::mat U; arma::vec d; arma::mat V;
    svd(U,d,V,tmp);
    d.subvec(numClust, d.size()-1) = arma::zeros<arma::vec>(d.size()-numClust);
    Q = U * diagmat(d) * V.t();

    //update S;
    S = nonnegASC_c(Y/mu + Q);

    //update w
    arma::mat Dmat = arma::eye<arma::mat>(n,n);
    arma::vec dvec(n);
    for (int i=0; i<n; ++i){
      dvec(i) = norm(Q-P.slice(i), "fro");
      dvec(i) = dvec(i) * dvec(i) / 2;
    }
    arma::mat Amat = arma::mat(n, n+1);
    Amat.col(0) = arma::ones<arma::vec>(n);
    Amat.submat(0,1,n-1,n) = arma::eye<arma::mat>(n,n);
    arma::vec bvec = arma::zeros<arma::vec>(n+1);
    bvec(0) = 1;
    Rcpp::List qp = solve_QP(Dmat, -dvec, Amat, bvec, 1);
    w = as<arma::vec>(qp["solution"]);

    //update Y
    Y = Y + mu * (S - Q);

    //update mu
    mu = 1;
    //mu = min(rho*mu, 1e+10);
    cout << w.t() << "\n";
  }
  arma::mat U; arma::vec s; arma::mat V;
  svd(U,s,V,S);
  arma::vec pi = V.col(0);
  arma::vec Dist = pi / accu(pi) + 1e-10;
  arma::vec invDist = 1/Dist;
  arma::mat Distmat = diagmat(Dist);
  arma::mat invDistmat = diagmat(invDist);
  arma::mat sspi = sqrt(invDistmat);
  arma::mat spi = sqrt(Distmat);
  S = (spi * S * sspi + sspi * S.t() * spi)/2;
  return Rcpp::List::create(
    Rcpp::Named("Q") = Q,
    Rcpp::Named("S") = S,
    Rcpp::Named("f") = funV,
    Rcpp::Named("w") = w
  );
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
Rcpp::List nuclear_objective_c(const arma::cube P,
                               const double tau,
                               double mu=1e-3,
                               const double rho = 1.9,
                               const int max_iter=100,
                               const double eps = 1e-9,
                               const bool verbose = false){
  int m = P.n_rows; int p = P.n_cols; int n = P.n_slices;
  arma::vec funV   = arma::zeros<arma::vec>(max_iter);
  arma::vec w      = arma::ones<arma::vec>(n)* 1/n;
  arma::mat Q      = arma::zeros<arma::mat>(m,p);
  arma::mat Y      = arma::zeros<arma::mat>(m,p);
  arma::mat S      = sum(P, 2)/n;
  arma::mat S_old  = arma::zeros<arma::mat>(m,p);
  int step = 0;
  while(1){
    step += 1;
    if(step > max_iter-1){
      cout << "reached max iteration \n";
      break;
    }
    double max_inf_norm = norm(S-Q, "inf");
    arma::mat U; arma::vec d; arma::mat V;
    svd(U,d,V,S);
    funV(step-1) = sum(w%w);
    for (int i=0; i < n; ++i){
      funV(step-1) += norm(Q-P.slice(i), "fro") * w(i);
    }
    funV(step-1) += sum(d);
    double relChg = norm(S-S_old, "fro")/max(1.0, norm(S_old, "fro"));
    S_old = S;
    if(verbose){
      cout << "iter" << step <<
        ": \n max_inf_norm = " << max_inf_norm <<
        "\n relChg = " << relChg <<
        "\n mu = " << mu <<
        "\n funV = " << funV(step-1) << "\n";
    }
    if(step > 1 & max_inf_norm < eps & relChg < eps){
      cout << "converged!\n";
      break;
    }

    //Update Qi;
    arma::mat M = mu*S + Y;
    for (int i=0; i<n; ++i){
      M += w(i) * P.slice(i);
    }
    M = M/(mu+2);
    double C = tau/(mu+2);
    svd(U,d,V,M);
    arma::uvec svp2 = find(d > C);
    int svp = svp2.size();
    if(svp >=2){
      arma::vec Sigma2 = d.subvec(0,svp-1) - C;
      Q = U.cols(0,svp-1) * diagmat(Sigma2) * V.cols(0,svp-1).t();
    }else if(svp==1){
      double Sigma2 = d(0) - C;
      Q = U.col(0) * V.col(0).t() * Sigma2;
    }else{
      svp = 1;
      Q = arma::zeros<arma::mat>(m,m);
    }

    //update S;
    S = nonnegASC_c(Q-Y/mu);

    //update w
    arma::mat Dmat = arma::eye<arma::mat>(n,n);
    arma::vec dvec(n);
    for (int i=0; i<n; ++i){
      dvec(i) = norm(S-P.slice(i), "fro");
      dvec(i) = dvec(i) * dvec(i);
    }
    arma::mat Amat = arma::mat(n, n+1);
    Amat.col(0) = arma::ones<arma::vec>(n);
    Amat.submat(0,1,n-1,n) = arma::eye<arma::mat>(n,n);
    arma::vec bvec = arma::zeros<arma::vec>(n+1);
    bvec(0) = 1;
    Rcpp::List qp = solve_QP(Dmat, -dvec, Amat, bvec, 1);
    w = as<arma::vec>(qp["solution"]);

    //update Y
    Y = Y + mu * (S - Q);

    //update mu
    // mu = 1;
    mu = min(rho*mu, 1e+10);
  }

  arma::mat U; arma::vec s; arma::mat V;
  svd(U,s,V,S);
  arma::vec pi = V.col(0);
  arma::vec Dist = pi / accu(pi) + 1e-10;
  arma::vec invDist = 1/Dist;
  arma::mat Distmat = diagmat(Dist);
  arma::mat invDistmat = diagmat(invDist);
  arma::mat sspi = sqrt(invDistmat);
  arma::mat spi = sqrt(Distmat);
  S = (spi * S * sspi + sspi * S.t() * spi)/2;

  return Rcpp::List::create(
    Rcpp::Named("S") = S,
    Rcpp::Named("f") = funV,
    Rcpp::Named("w") = w
  );
}

