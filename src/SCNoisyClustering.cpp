#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::sp_mat sparse_sum(Rcpp::List P, int length){
  arma::sp_mat out = as<arma::sp_mat>(P[0]);
  for(int i=1; i < length; ++i){
    out = out + as<arma::sp_mat>(P[i]);
  }
  return out;
}

// // [[Rcpp::export]]
// arma::mat dist_c(arma::mat X){
//   int n = X.n_rows;
//   arma::mat out = arma::zeros<arma::mat>(n,n);
//   for (int i=0; i<(n-1); ++i){
//     for (int j=(i+1); j<n; ++j){
//       arma::rowvec diff = X.row(i) - X.row(j);
//       out(i,j) = sqrt(accu(diff%diff));
//       out(j,i) = out(i,j);
//     }
//   }
//   return out;
// }

// [[Rcpp::export]]
arma::vec QminusPslice(arma::sp_mat Q,
                       Rcpp::List P,
                       int length){
  arma::vec out = arma::zeros<arma::vec>(length);
  for (int i = 1; i < length; ++i){
    arma::sp_mat Pmat = as<arma::sp_mat>(P[i]);
    out(i) = accu((Q-Pmat)%(Q-Pmat));
  }
  return out;
}
// [[Rcpp::export]]
arma::mat proj_c(arma::mat X, arma::mat ref){
  return((X.t() * ref / X.n_rows).t());
}


// [[Rcpp::export]]
arma::mat armapmax(arma::mat A, double bound){
  int n = A.n_rows;
  int p = A.n_cols;
  // A.elem(find_nonfinite(A)).zeros();
  arma::mat out = arma::mat(n,p);
  for (int i=0; i<n; ++i){
    for (int j=0; j<p; ++j){
      out(i,j) = max(A(i,j), bound);
    }
  }
  return out;
}

// [[Rcpp::export]]
arma::mat armapmin(arma::mat A, double bound){
  int n = A.n_rows;
  int p = A.n_cols;
  // A.elem(find_nonfinite(A)).zeros();
  arma::mat out = arma::mat(n,p);
  for (int i=0; i<n; ++i){
    for (int j=0; j<p; ++j){
      out(i,j) = min(A(i,j), bound);
    }
  }
  return out;
}


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







// [[Rcpp::export]]
arma::sp_mat get_kernel_matrix(arma::mat X,
                               arma::mat Diff,
                               int k,
                               double sigma){
  int N = X.n_cols;

  //compute and sort Diff
  arma::mat Diff_sort(N,N);
  for (int j = 0; j < N; ++j){
    Diff_sort.row(j) = sort(Diff.row(j));
  }
  arma::mat TT(N, 1);
  TT.col(0) = mean(Diff_sort.cols(1, k), 1);
  arma::mat Sig = arma::zeros<arma::mat>(N,N);
  Sig.each_col() = TT.col(0);
  Sig = (Sig + Sig.t())/2;
  // Sig.elem(find(Sig < arma::datum::eps)).zeros();
  // Sig = Sig + arma::datum::eps;
  // arma::mat Sig_valid = arma::zeros<arma::mat>(N,N);
  // Sig_valid.elem(find(Sig > arma::datum::eps)).ones();
  // Sig = Sig % Sig_valid + arma::datum::eps;
  arma::mat W = (normpdf(Diff, arma::zeros<arma::mat>(N,N), sigma*Sig));
  arma::mat K = (W + W.t())/2;
  arma::vec dinv = 1/sqrt(K.diag()+1);
  arma::mat G = K % (dinv * dinv.t());
  arma::mat G1 = arma::mat(N,N);
  G1.each_col() = G.diag();
  arma::mat D_Kernels = (G1 + G1.t() - 2*G)/2;
  arma::sp_mat out = arma::zeros<arma::sp_mat>(N,N);
  arma::mat distX = D_Kernels;
  arma::mat distX1(N,N);
  arma::umat idx(N,N);
  for (int j=0; j < N; ++j){
    distX1.row(j) = sort(distX.row(j));
    idx.row(j)    = sort_index(distX.row(j)).t();
  }
  arma::sp_mat A = arma::zeros<arma::sp_mat>(N,N);
  arma::mat di = distX1.cols(1,k+1);
  arma::umat id = idx.cols(1,k+1);
  arma::mat numerator(N, k+1);
  numerator.each_col() = di.col(k);
  numerator = numerator - di;
  arma::vec temp = k * di.col(k) - sum(di.cols(0,k-1), 1);
  arma::mat denominator(N, k+1);
  denominator.each_col() = temp;
  arma::mat div = numerator / denominator;
  arma::mat a(N, k+1);
  a.each_col() = arma::linspace(0, N-1, N);
  //cout << "check 7" << "\n";
  for (int row = 0; row < N; ++row){
    for (int col = 0; col < k+1; ++col){
      if(div(row,col) != arma::datum::inf){
        A(a(row,col), id(row,col)) = div(row,col);
      }
    }
  }
  //cout << "check 8" << "\n";
  out = (A + A.t())/2;
  return out;
}


// [[Rcpp::export]]
Rcpp::List sparse_scaledlasso_list_c(Rcpp::List P,
                                     int n,
                                     const double tau,
                                     const double gamma,
                                     double mu=1e-3,
                                     const double rho = 1.9,
                                     const int max_iter=100,
                                     const double eps = 1e-9,
                                     const bool verbose = false){

  arma::sp_mat samplemat = as<arma::sp_mat>(P[0]);
  int m = samplemat.n_rows;
  int p = samplemat.n_cols;

  arma::vec funV        = arma::zeros<arma::vec>(max_iter);
  arma::vec sigma       = arma::ones<arma::vec>(n);
  arma::mat Y           = arma::zeros<arma::mat>(m,p);
  arma::sp_mat Q        = arma::zeros<arma::sp_mat>(m,p);
  arma::sp_mat S        = sparse_sum(P,n);
  arma::sp_mat S_old    = arma::zeros<arma::sp_mat>(m,p);
  int step = 0;

  while(1){
    if(step > max_iter-1){
      cout << "reached max iteration \n";
      break;
    }
    step += 1;
    double max_inf_norm = norm(S-Q, "inf");
    double relChg = norm(S-S_old, "fro")/max(1.0, norm(S_old, "fro"));
    S_old = S;
    if(verbose & (step % 5==0)){
      cout << "iter" << step << ":"
      "\n relative change = " << relChg << "\n";
    }
    if(step > 1 & max_inf_norm < eps & relChg < eps){
      break;
    }

    //Update Qi;
    arma::mat M = mu*S + Y;
    double phi = sum(1/(m*m*sigma)) + mu + gamma;
    for (int i=0; i<n; ++i){
      M += as<arma::sp_mat>(P[i])/(m*m*sigma(i));
    }
    M = M/phi;
    double C = tau/phi;
    Q = armapmax(M-C, 0) + armapmin(M-C,0);

    //update S;
    S = arma::conv_to<arma::sp_mat>::from(nonnegASC_c(Q-Y/mu));

    //update sigma
    for (int i=0; i < n; ++i){
      sigma(i) = norm(as<arma::sp_mat>(P[i]) - Q, "fro");
      sigma(i) /= (m*m);
      sigma(i) = 1/sigma(i);
    }

    //update Y
    Y = Y + mu * (S - Q);

    //update mu
    mu = min(rho*mu, 1e+10);
  }

  arma::mat U; arma::vec s; arma::mat V;
  svd(U,s,V, arma::conv_to<arma::mat>::from(S));
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
    // Rcpp::Named("f") = funV,
    Rcpp::Named("sigma") = sigma
    // Rcpp::Named("L1_norm") = L1_norm,
    // Rcpp::Named("F_norm") = F_norm,
    // Rcpp::Named("fro_norm") = fro_norm,
    // Rcpp::Named("noise") = noise
  );
}

// [[Rcpp::export]]
arma::mat tsne_c(arma::mat X,
                 arma::mat initial_config,
                 int k,
                 int max_iter=1000,
                 double min_cost = 1e-3,
                 int epoch = 100){

  double momentum = 0.8;  //double final_momentum = 0.8;
  double epsilon = 500;   //double mom_switch_iter = 250;
  double min_gain = 0.01; double initial_P_gain = 4;
  double eps = arma::datum::eps;
  arma::vec gaincheck = arma::zeros<arma::vec>(1000);
  arma::mat ydata = initial_config;
  initial_P_gain = 1;

  arma::mat P = (X + X.t())/2;
  P.elem(find(P<eps)).fill(eps);
  P = P / accu(P);
  P = P*initial_P_gain;
  arma::mat grads = arma::zeros<arma::mat>(ydata.n_rows, ydata.n_cols);
  arma::mat incs = arma::zeros<arma::mat>(ydata.n_rows, ydata.n_cols);
  arma::mat gains = arma::ones<arma::mat>(ydata.n_rows, ydata.n_cols);
  arma::mat Q; arma::mat stiffness;
  arma::vec costs = arma::zeros<arma::vec>(max_iter);
  arma::mat ydata_old = ydata;
  for (int niter = 0; niter < max_iter; niter++){

    arma::vec sum_ydata = sum(ydata%ydata, 1);
    arma::mat tmp1 = -2 * ydata * ydata.t();
    arma::mat tmp2 = tmp1;
    for (int i=0; i < tmp2.n_cols; ++i){
      tmp2.row(i) = tmp2.row(i) + sum_ydata.t();
    }
    arma::mat tmp3 = tmp2;
    for (int i=0; i < tmp3.n_cols; ++i){
      tmp3.col(i) = tmp3.col(i) + sum_ydata + arma::ones<arma::vec>(sum_ydata.size());
    }
    arma::mat num = 1/(tmp3);
    num.diag().fill(0);
    Q = num/accu(num);
    Q.elem(find(Q<eps)).fill(eps);
    arma::mat stiffness = (P-Q) % num;
    grads = 4 * (diagmat(sum(stiffness, 0)) - stiffness) * ydata;
    gains = (gains + 0.2) % abs(sign(grads) != sign(incs)) +
      gains*0.8 % abs(sign(grads) == sign(incs));
    gains.elem(find(gains < min_gain)).fill(min_gain);

    incs = momentum*incs - epsilon*(gains%grads);
    ydata = ydata + incs;
    for (int i=0; i < ydata.n_cols; ++i){
      arma::vec tmpy = ydata.col(i);
      tmpy = tmpy - mean(tmpy);
      ydata.col(i) = tmpy;
    }
    ydata.elem(find(ydata<-100)).fill(-100);
    ydata.elem(find(ydata>100)).fill(100);

    //convergence criterion
    if(niter % epoch == 99){
      double cost = norm(ydata-ydata_old, "fro") / (ydata.n_rows * ydata.n_cols);
      if(cost < min_cost){
        break;
      }
    }
    ydata_old = ydata;
  }
  return ydata;
}

// // [[Rcpp::export]]
// arma::vec tsne_spectral_c(arma::mat A, int numClust, int numEigen = 0){
//   if(numEigen==0){numEigen = numClust;}
//   arma::vec rs = sum(A, 1);
//   arma::mat L  = diagmat(1/sqrt(rs)) * (diagmat(rs)-A) * diagmat(1/sqrt(rs));
//   arma::vec D; arma::mat U;
//   eig_sym(D, U, L);
//   arma::vec U_index = arma::linspace(U.n_cols, U.n_cols-numEigen+1);
//   arma::mat tmpU(A.n_cols, numEigen);
//   for (int i=0; i < numEigen-1; ++i){
//     tmpU.col(i) = U.col(U_index(i)-1);
//   }
//   arma::mat F_last = tsne_c(A, tmpU, numEigen);
//   arma::mat means;
//   kmeans(means, F_last, numClust, arma::random_subset, 10, false);
//   return means;
//
// }

// // [[Rcpp::depends("RcppArmadillo")]]
// // [[Rcpp::export]]
// Rcpp::List RMSC_c(arma::cube T,
//                   double lambda,
//                   double mu=1e-3,
//                   double rho = 1.9,
//                   int max_iter=100,
//                   double eps = 1e-9,
//                   bool verbose = false){
//   int m = T.n_rows; int p = T.n_cols; int n = T.n_slices;
//   if(m!=p){
//     Rcerr << "input matrix T must be square transition matrix\n";
//   }
//   arma::mat Z  = arma::zeros<arma::mat>(m,p);
//   arma::cube E = arma::zeros<arma::cube>(m,p,n);
//   E.randn();
//   arma::cube Y = arma::zeros<arma::cube>(m,p,n);
//   arma::cube B = arma::zeros<arma::cube>(m,p,n);
//   arma::mat Q  = arma::zeros<arma::mat>(m,p);
//   arma::mat P  = arma::zeros<arma::mat>(m,p);
//   arma::mat P_old = arma::zeros<arma::mat>(m,p);
//   P_old.randn();
//   int step = 0;
//   while(1){
//     step += 1;
//     double max_inf_norm = -1;
//     for (int i=0; i < n; ++i){
//       arma::mat Ti    = T.slice(i); arma::mat Ei = E.slice(i);
//       arma::mat diff  = Ti-Ei-P;
//       double inf_norm = norm(diff, "inf");
//       max_inf_norm    = max(max_inf_norm, inf_norm);
//     }
//     arma::vec s;
//     svd(s, P);
//     double funV = sum(s) + lambda * accu(abs(E));
//     double relChg = norm(P-P_old, "fro")/max(1.0, norm(P_old, "fro"));
//     P_old = P;
//     arma::mat tmp = P-Q;
//     double max_inf_norm2 = norm(tmp, "inf");
//     if(verbose && step%10==0){
//       cout << "iter" << step <<
//         ": \n max_inf_norm = " << max_inf_norm <<
//           "\n max_inf_norm2 = " << max_inf_norm2 <<
//             "\n relChg = " << relChg <<
//               "\n mu = " << mu <<
//                 "\n funV = " << funV << "\n";
//     }
//     if(step > 1 & max_inf_norm < eps){
//       break;
//     }
//     if(step > max_iter){
//       cout << "reached max iteration \n";
//       break;
//     }
//
//     //Update P
//     arma::mat temp = sum(T-E-Y/mu, 2);
//     arma::mat B = (Q- Z/mu + temp)/(n+1);
//     P = nonnegASC_c(B);
//     for (int i=0; i < m; ++i){
//       if(abs(sum(P.row(i)) - 1) > 1e-10){
//         Rcerr << "rowsum 1 error";
//       }
//     }
//
//     //Update Q;
//     arma::mat M = P+Z/mu;
//     double C = 1/mu;
//     arma::mat U; arma::vec Sigma; arma::mat V;
//     svd(U,Sigma,V,M);
//     arma::uvec svp2 = find(Sigma > C);
//     int svp = svp2.size();
//     if(svp >=2){
//       arma::vec Sigma2 = Sigma.subvec(0,svp-1) - C;
//       Q = U.cols(0,svp-1) * diagmat(Sigma2) * V.cols(0,svp-1).t();
//     }else if(svp==1){
//       double Sigma2 = Sigma(0) - C;
//       Q = U.col(0) * V.col(0).t() * Sigma2;
//     }else{
//       svp = 1;
//       Q = arma::zeros<arma::mat>(m,m);
//     }
//
//     //Update Ei and Yi
//     for (int i=0; i < n; ++i){
//       arma::mat C = T.slice(i) - P - Y.slice(i) / mu;
//       E.slice(i) = armapmax(C-lambda/mu, 0) + armapmin(C+lambda/mu,0);
//       Y.slice(i) = Y.slice(i) + mu * (P+E.slice(i)-T.slice(i));
//     }
//     //Update Z
//     Z = Z + mu*(P-Q);
//     mu = min(rho*mu, 1e+10);
//   }
//   arma::mat U; arma::vec s; arma::mat V;
//   svd(U,s,V,P.t());
//   arma::vec pi = V.col(0);
//   arma::vec Dist = pi / accu(pi);
//   arma::vec invDist = 1/Dist;
//   arma::mat Distmat = diagmat(Dist);
//   arma::mat invDistmat = diagmat(invDist);
//   arma::mat sspi = sqrt(invDistmat);
//   arma::mat spi = sqrt(Distmat);
//   P = (spi * P * sspi + sspi * P.t() * spi)/2;
//   return Rcpp::List::create(
//     Rcpp::Named("S") = P,
//     Rcpp::Named("E") = E
//   );
// }
//
// // [[Rcpp::depends("RcppArmadillo")]]
// // [[Rcpp::export]]
// Rcpp::List objective_1_c(const arma::cube P,
//                          const int numClust,
//                          double mu=1e-3,
//                          const double rho = 1.9,
//                          const int max_iter=100,
//                          const double eps = 1e-9,
//                          const bool verbose = false){
//   int m = P.n_rows; int p = P.n_cols; int n = P.n_slices;
//   arma::vec  funV  = arma::zeros<arma::vec>(max_iter);
//   arma::vec  fnorm = arma::zeros<arma::vec>(max_iter);
//   arma::vec  w     = arma::ones<arma::vec>(n)* 1/n;
//   arma::cube Q     = arma::zeros<arma::cube>(m,p,n);
//   arma::cube Y     = arma::zeros<arma::cube>(m,p,n);
//   arma::mat  S     = arma::zeros<arma::mat>(m,p);
//   arma::mat  S_old = arma::zeros<arma::mat>(m,p);
//   S_old.randn();
//   int step = 0;
//   while(1){
//     step += 1;
//     if(step > max_iter-1){
//       cout << "reached max iteration \n";
//       break;
//     }
//     double max_inf_norm = -1;
//     for (int i=0; i < n; ++i){
//       arma::mat Qi = Q.slice(i);
//       arma::mat diff = S-Qi;
//       double inf_norm = norm(diff, "inf");
//       max_inf_norm = max(max_inf_norm, inf_norm);
//     }
//     funV(step-1) = sum(w%w);
//     fnorm(step-1) = 0;
//     for (int i=0; i < n; ++i){
//       funV(step-1) += norm(Q.slice(i)-P.slice(i), "fro") * norm(Q.slice(i)-P.slice(i), "fro") * w(i);
//       fnorm(step-1) += norm(Q.slice(i)-P.slice(i), "fro") * norm(Q.slice(i)-P.slice(i), "fro") * w(i);
//     }
//     double relChg = norm(S-S_old, "fro")/max(1.0, norm(S_old, "fro"));
//     S_old = S;
//     if(verbose){
//       cout << "iter" << step <<
//         ": \n max_inf_norm = " << max_inf_norm <<
//           "\n relChg = " << relChg <<
//             "\n mu = " << mu <<
//               "\n funV = " << funV(step-1) << "\n";
//     }
//     if(step > 1 & max_inf_norm < eps & relChg < eps){
//       break;
//     }
//
//     //Update Qi;
//     for (int i = 0; i < n; ++i){
//       arma::mat tmp = mu*S + P.slice(i) * w(i) + Y.slice(i);
//       tmp = tmp / (w(i) + mu);
//       arma::mat U; arma::vec d; arma::mat V;
//       svd(U,d,V,tmp);
//       d.subvec(numClust, d.size()-1) = arma::zeros<arma::vec>(d.size()-numClust);
//       Q.slice(i) = U * diagmat(d) * V.t();
//     }
//
//     //update S;
//     arma::mat temp = sum(Q - Y/mu, 2) / n;
//     S = nonnegASC_c(temp);
//
//     //update w
//     arma::mat Dmat = arma::eye<arma::mat>(n,n);
//     arma::vec dvec(n);
//     for (int i=0; i<n; ++i){
//       dvec(i) = norm(Q.slice(i)-P.slice(i), "fro");
//       dvec(i) = -dvec(i) * dvec(i)/2;
//     }
//     arma::mat Amat = arma::mat(n, n+1);
//     Amat.col(0) = arma::ones<arma::vec>(n);
//     Amat.submat(0,1,n-1,n) = arma::eye<arma::mat>(n,n);
//     arma::vec bvec = arma::zeros<arma::vec>(n+1);
//     bvec(0) = 1;
//     Rcpp::List qp = solve_QP(Dmat, dvec, Amat, bvec, 1);
//     w = as<arma::vec>(qp["solution"]);
//
//     //update Y
//     for (int i = 0; i < n; ++i){
//       Y.slice(i) = Y.slice(i) + mu * (S - Q.slice(i));
//     }
//
//     //update mu
//     mu = min(rho*mu, 1e+10);
//   }
//   arma::mat U; arma::vec s; arma::mat V;
//   svd(U,s,V,S);
//   arma::vec pi = V.col(0);
//   arma::vec Dist = pi / accu(pi) + 1e-10;
//   arma::vec invDist = 1/Dist;
//   arma::mat Distmat = diagmat(Dist);
//   arma::mat invDistmat = diagmat(invDist);
//   arma::mat sspi = sqrt(invDistmat);
//   arma::mat spi = sqrt(Distmat);
//   S = (spi * S * sspi + sspi * S.t() * spi)/2;
//   return Rcpp::List::create(
//     Rcpp::Named("S") = S,
//     Rcpp::Named("f") = funV,
//     Rcpp::Named("fnorm") = fnorm,
//     Rcpp::Named("w") = w
//   );
// }
//
// // [[Rcpp::depends("RcppArmadillo")]]
// // [[Rcpp::export]]
// Rcpp::List objective_2_c(const arma::cube P,
//                          const int numClust,
//                          double mu=1e-3,
//                          const double rho = 1.9,
//                          const int max_iter=100,
//                          const double eps = 1e-9,
//                          const bool verbose = false){
//
//   int m = P.n_rows; int p = P.n_cols; int n = P.n_slices;
//   arma::vec funV = arma::zeros<arma::vec>(max_iter);
//   arma::vec sigma = arma::ones<arma::vec>(n);
//   arma::cube Q = arma::zeros<arma::cube>(m,p,n);
//   arma::cube Y = arma::zeros<arma::cube>(m,p,n);
//   arma::mat S = arma::zeros<arma::mat>(m,p);
//   arma::mat S_old = arma::zeros<arma::mat>(m,p);
//   S_old.randn();
//
//   int step = 0;
//   while(1){
//     arma::vec mmsigma = m*m*sigma;
//     step += 1;
//     if(step > max_iter-1){
//       cout << "reached max iteration \n";
//       break;
//     }
//     double max_inf_norm = -1;
//     for (int i=0; i < n; ++i){
//       arma::mat Qi = Q.slice(i);
//       arma::mat diff = S-Qi;
//       double inf_norm = norm(diff, "inf");
//       max_inf_norm = max(max_inf_norm, inf_norm);
//     }
//     funV(step-1) = accu(sigma/2);
//     for (int i=0; i < n; ++i){
//       funV(step-1) += norm(Q.slice(i)-P.slice(i), "fro")/(mmsigma(i));
//     }
//     double relChg = norm(S-S_old, "fro")/max(1.0, norm(S_old, "fro"));
//     S_old = S;
//     if(verbose){
//       cout << "iter" << step <<
//         ": \n max_inf_norm = " << max_inf_norm <<
//           "\n relChg = " << relChg <<
//             "\n mu = " << mu <<
//               "\n funV = " << funV(step-1) << "\n";
//     }
//     if(step > 1 & max_inf_norm < eps & relChg < eps){
//       break;
//     }
//
//     //Update Qi;
//     for (int i = 0; i < n; ++i){
//       arma::mat tmp = (mu*S + P.slice(i)/mmsigma(i) + Y.slice(i));
//       tmp = tmp *(mmsigma(i) / (mu*mmsigma(i) + 1));
//       arma::mat U; arma::vec d; arma::mat V;
//       svd(U,d,V,tmp);
//       d.subvec(numClust, d.size()-1) = arma::zeros<arma::vec>(d.size()-numClust);
//       Q.slice(i) = U * diagmat(d) * V.t();
//       //Q.slice(i) = Q.slice(i).each_col() / sum(Q.slice(i), 1);
//     }
//
//     //update S;
//     arma::mat tmp = sum(Q - Y/mu, 2) / n;
//     S = nonnegASC_c(tmp);
//
//     //update sigma
//     for (int i=0; i<n; ++i){
//       double tmpnum = norm(Q.slice(i)-P.slice(i), "fro");
//       sigma(i) = tmpnum*tmpnum/(m*m);
//     }
//
//     //update Y
//     for (int i = 0; i < n; ++i){
//       Y.slice(i) = Y.slice(i) + mu * (S - Q.slice(i));
//     }
//     //update mu
//     mu = min(rho*mu, 1e+10);
//   }
//   arma::mat U; arma::vec s; arma::mat V;
//   svd(U,s,V,S);
//   arma::vec pi = V.col(0);
//   arma::vec Dist = pi / accu(pi) + 1e-10;
//   arma::vec invDist = 1/Dist;
//   arma::mat Distmat = diagmat(Dist);
//   arma::mat invDistmat = diagmat(invDist);
//   arma::mat sspi = sqrt(invDistmat);
//   arma::mat spi = sqrt(Distmat);
//   S = (spi * S * sspi + sspi * S.t() * spi)/2;
//   return Rcpp::List::create(
//     Rcpp::Named("S") = S,
//     Rcpp::Named("f") = funV,
//     Rcpp::Named("sigma") = sigma
//   );
// }
//
//
//
//
// // [[Rcpp::depends("RcppArmadillo")]]
// // [[Rcpp::export]]
// Rcpp::List objective_1_c_one_aux(const arma::cube P,
//                                  const int numClust,
//                                  double mu=1e-3,
//                                  const double rho = 1.9,
//                                  const int max_iter=100,
//                                  const double eps = 1e-9,
//                                  const bool verbose = false){
//   int m = P.n_rows; int p = P.n_cols; int n = P.n_slices;
//   arma::vec funV   = arma::zeros<arma::vec>(max_iter);
//   arma::vec w      = arma::ones<arma::vec>(n)* 1/n;
//   arma::mat Q      = arma::zeros<arma::mat>(m,p);
//   arma::mat Y      = arma::zeros<arma::mat>(m,p);
//   arma::mat S      = arma::zeros<arma::mat>(m,p);
//   arma::mat S_old  = arma::zeros<arma::mat>(m,p);
//   S_old.randn();
//   S.randn();
//   S = nonnegASC_c(S);
//   int step = 0;
//   while(1){
//     step += 1;
//     if(step > max_iter-1){
//       cout << "reached max iteration \n";
//       break;
//     }
//     double max_inf_norm = norm(S-Q, "inf");
//     funV(step-1) = sum(w%w);
//     for (int i=0; i < n; ++i){
//       funV(step-1) += norm(Q-P.slice(i), "fro") * w(i);
//     }
//     double relChg = norm(S-S_old, "fro")/max(1.0, norm(S_old, "fro"));
//     S_old = S;
//     if(verbose){
//       cout << "iter" << step <<
//         ": \n max_inf_norm = " << max_inf_norm <<
//           "\n relChg = " << relChg <<
//             "\n mu = " << mu <<
//               "\n funV = " << funV(step-1) << "\n";
//     }
//     if(step > 1 & max_inf_norm < eps & relChg < eps){
//       break;
//     }
//
//     //Update Qi;
//     arma::mat tmp = S*mu + Y;
//     for (int i=0; i < n; ++i){
//       tmp = tmp + P.slice(i) * w(i);
//     }
//     tmp = tmp / (mu + 1);
//     arma::mat U; arma::vec d; arma::mat V;
//     svd(U,d,V,tmp);
//     d.subvec(numClust, d.size()-1) = arma::zeros<arma::vec>(d.size()-numClust);
//     Q = U * diagmat(d) * V.t();
//
//     //update S;
//     S = nonnegASC_c(Y/mu + Q);
//
//     //update w
//     arma::mat Dmat = arma::eye<arma::mat>(n,n);
//     arma::vec dvec(n);
//     for (int i=0; i<n; ++i){
//       dvec(i) = norm(Q-P.slice(i), "fro");
//       dvec(i) = dvec(i) * dvec(i) / 2;
//     }
//     arma::mat Amat = arma::mat(n, n+1);
//     Amat.col(0) = arma::ones<arma::vec>(n);
//     Amat.submat(0,1,n-1,n) = arma::eye<arma::mat>(n,n);
//     arma::vec bvec = arma::zeros<arma::vec>(n+1);
//     bvec(0) = 1;
//     Rcpp::List qp = solve_QP(Dmat, -dvec, Amat, bvec, 1);
//     w = as<arma::vec>(qp["solution"]);
//
//     //update Y
//     Y = Y + mu * (S - Q);
//
//     //update mu
//     mu = 1;
//     //mu = min(rho*mu, 1e+10);
//     cout << w.t() << "\n";
//   }
//   arma::mat U; arma::vec s; arma::mat V;
//   svd(U,s,V,S);
//   arma::vec pi = V.col(0);
//   arma::vec Dist = pi / accu(pi) + 1e-10;
//   arma::vec invDist = 1/Dist;
//   arma::mat Distmat = diagmat(Dist);
//   arma::mat invDistmat = diagmat(invDist);
//   arma::mat sspi = sqrt(invDistmat);
//   arma::mat spi = sqrt(Distmat);
//   S = (spi * S * sspi + sspi * S.t() * spi)/2;
//   return Rcpp::List::create(
//     Rcpp::Named("Q") = Q,
//     Rcpp::Named("S") = S,
//     Rcpp::Named("f") = funV,
//     Rcpp::Named("w") = w
//   );
// }
//
//
// // [[Rcpp::depends("RcppArmadillo")]]
// // [[Rcpp::export]]
// Rcpp::List objective_3_c(const arma::cube P,
//                          const double tau,
//                          double mu=1e-3,
//                          const double rho = 1.9,
//                          const int max_iter=100,
//                          const double eps = 1e-9,
//                          const bool verbose = false){
//   int m = P.n_rows; int p = P.n_cols; int n = P.n_slices;
//   arma::vec funV   = arma::zeros<arma::vec>(max_iter);
//   arma::vec w      = arma::ones<arma::vec>(n)* 1/n;
//   arma::mat Q      = arma::zeros<arma::mat>(m,p);
//   arma::mat Y      = arma::zeros<arma::mat>(m,p);
//   arma::mat S      = sum(P, 2)/n;
//   arma::mat S_old  = arma::zeros<arma::mat>(m,p);
//   int step = 0;
//   while(1){
//     step += 1;
//     if(step > max_iter-1){
//       cout << "reached max iteration \n";
//       break;
//     }
//     double max_inf_norm = norm(S-Q, "inf");
//     arma::mat U; arma::vec d; arma::mat V;
//     svd(U,d,V,S);
//     funV(step-1) = sum(w%w);
//     for (int i=0; i < n; ++i){
//       funV(step-1) += norm(Q-P.slice(i), "fro") * w(i);
//     }
//     funV(step-1) += sum(d);
//     double relChg = norm(S-S_old, "fro")/max(1.0, norm(S_old, "fro"));
//     S_old = S;
//     if(verbose){
//       cout << "iter" << step <<
//         ": \n max_inf_norm = " << max_inf_norm <<
//           "\n relChg = " << relChg <<
//             "\n mu = " << mu <<
//               "\n funV = " << funV(step-1) << "\n";
//     }
//     if(step > 1 & max_inf_norm < eps & relChg < eps){
//       break;
//     }
//
//     //Update Qi;
//     arma::mat M = mu*S + Y;
//     for (int i=0; i<n; ++i){
//       M += w(i) * P.slice(i);
//     }
//     M = M/(mu+2);
//     double C = tau/(mu+2);
//     svd(U,d,V,M);
//     arma::uvec svp2 = find(d > C);
//     int svp = svp2.size();
//     if(svp >=2){
//       arma::vec Sigma2 = d.subvec(0,svp-1) - C;
//       Q = U.cols(0,svp-1) * diagmat(Sigma2) * V.cols(0,svp-1).t();
//     }else if(svp==1){
//       double Sigma2 = d(0) - C;
//       Q = U.col(0) * V.col(0).t() * Sigma2;
//     }else{
//       svp = 1;
//       Q = arma::zeros<arma::mat>(m,m);
//     }
//
//     //update S;
//     S = nonnegASC_c(Q-Y/mu);
//
//     //update w
//     arma::mat Dmat = arma::eye<arma::mat>(n,n);
//     arma::vec dvec(n);
//     for (int i=0; i<n; ++i){
//       dvec(i) = norm(S-P.slice(i), "fro");
//       dvec(i) = dvec(i) * dvec(i);
//     }
//     arma::mat Amat = arma::mat(n, n+1);
//     Amat.col(0) = arma::ones<arma::vec>(n);
//     Amat.submat(0,1,n-1,n) = arma::eye<arma::mat>(n,n);
//     arma::vec bvec = arma::zeros<arma::vec>(n+1);
//     bvec(0) = 1;
//     Rcpp::List qp = solve_QP(Dmat, -dvec, Amat, bvec, 1);
//     w = as<arma::vec>(qp["solution"]);
//
//     //update Y
//     Y = Y + mu * (S - Q);
//
//     //update mu
//     // mu = 1;
//     mu = min(rho*mu, 1e+10);
//   }
//
//   arma::mat U; arma::vec s; arma::mat V;
//   svd(U,s,V,S);
//   arma::vec pi = V.col(0);
//   arma::vec Dist = pi / accu(pi) + 1e-10;
//   arma::vec invDist = 1/Dist;
//   arma::mat Distmat = diagmat(Dist);
//   arma::mat invDistmat = diagmat(invDist);
//   arma::mat sspi = sqrt(invDistmat);
//   arma::mat spi = sqrt(Distmat);
//   S = (spi * S * sspi + sspi * S.t() * spi)/2;
//
//   return Rcpp::List::create(
//     Rcpp::Named("S") = S,
//     Rcpp::Named("f") = funV,
//     Rcpp::Named("w") = w
//   );
// }
//
//
//
//
// // [[Rcpp::depends("RcppArmadillo")]]
// // [[Rcpp::export]]
// Rcpp::List nuclear_scaledlasso_c(const arma::cube P,
//                                  const double tau,
//                                  double mu=1e-3,
//                                  const double rho = 1.9,
//                                  const int max_iter=100,
//                                  const double eps = 1e-9,
//                                  const bool verbose = false){
//   int m = P.n_rows; int p = P.n_cols; int n = P.n_slices;
//   arma::vec funV     = arma::zeros<arma::vec>(max_iter);
//   arma::vec sigma    = arma::ones<arma::vec>(n);
//   arma::mat Y        = arma::zeros<arma::mat>(m,p);
//   arma::mat Q        = arma::zeros<arma::mat>(m,p);
//   arma::mat S        = sum(P, 2)/n;
//   arma::mat S_old    = arma::zeros<arma::mat>(m,p);
//   arma::vec nuc_norm = arma::zeros<arma::vec>(max_iter);
//   arma::vec fro_norm = arma::zeros<arma::vec>(max_iter);
//   arma::vec noise    = arma::zeros<arma::vec>(max_iter);
//   int step = 0;
//   while(1){
//     if(step > max_iter-1){
//       cout << "reached max iteration \n";
//       break;
//     }
//     step += 1;
//     double max_inf_norm = norm(S-Q, "inf");
//     arma::mat U; arma::vec d; arma::mat V;
//     svd(U,d,V,S);
//     noise(step-1) = sum(sigma)/2;
//     for (int i=0; i < n; ++i){
//       fro_norm(step-1) += norm(Q-P.slice(i), "fro")/(sigma(i) * m * m);
//     }
//     nuc_norm(step-1) += sum(d);
//     funV(step-1) = nuc_norm(step-1) + fro_norm(step-1) + noise(step-1);
//     double relChg = norm(S-S_old, "fro")/max(1.0, norm(S_old, "fro"));
//     S_old = S;
//     if(verbose){
//       cout << "iter" << step <<
//         ": \n max_inf_norm = " << max_inf_norm <<
//           "\n relChg = " << relChg <<
//             "\n mu = " << mu <<
//               "\n funV = " << funV(step-1) << "\n";
//     }
//     if(step > 1 & max_inf_norm < eps & relChg < eps){
//       break;
//     }
//
//     //Update Qi;
//     arma::mat M = mu*S + Y;
//     double phi = sum(1/(m*m*sigma)) + mu;
//     for (int i=0; i<n; ++i){
//       M += P.slice(i)/(m*m*sigma(i));
//     }
//     M = M/phi;
//     double C = tau/phi;
//     svd(U,d,V,M);
//     arma::uvec svp2 = find(d > C);
//     int svp = svp2.size();
//     if(svp >=2){
//       arma::vec Sigma2 = d.subvec(0,svp-1) - C;
//       Q = U.cols(0,svp-1) * diagmat(Sigma2) * V.cols(0,svp-1).t();
//     }else if(svp==1){
//       double Sigma2 = d(0) - C;
//       Q = U.col(0) * V.col(0).t() * Sigma2;
//     }else{
//       svp = 1;
//       Q = arma::zeros<arma::mat>(m,m);
//     }
//
//     //update S;
//     S = nonnegASC_c(Q-Y/mu);
//
//     //update sigma
//     for (int i=0; i < n; ++i){
//       sigma(i) = norm(P.slice(i) - Q, "fro");
//       sigma(i) /= (m*m);
//     }
//
//     //update Y
//     Y = Y + mu * (S - Q);
//
//     //update mu
//     mu = min(rho*mu, 1e+10);
//   }
//
//   arma::mat U; arma::vec s; arma::mat V;
//   svd(U,s,V,S);
//   arma::vec pi = V.col(0);
//   arma::vec Dist = pi / accu(pi) + 1e-10;
//   arma::vec invDist = 1/Dist;
//   arma::mat Distmat = diagmat(Dist);
//   arma::mat invDistmat = diagmat(invDist);
//   arma::mat sspi = sqrt(invDistmat);
//   arma::mat spi = sqrt(Distmat);
//   S = (spi * S * sspi + sspi * S.t() * spi)/2;
//
//   return Rcpp::List::create(
//     Rcpp::Named("S") = S,
//     Rcpp::Named("f") = funV,
//     Rcpp::Named("sigma") = sigma,
//     Rcpp::Named("nuc_norm") = nuc_norm,
//     Rcpp::Named("fro_norm") = fro_norm,
//     Rcpp::Named("noise") = noise
//   );
// }
//

// // [[Rcpp::export]]
// arma::cube corr_kernel_c(arma::mat X,
//                          arma::mat Diff,
//                          arma::vec allk_input,
//                          arma::vec sigma_input,
//                          int k = 0){
//   int N = X.n_rows;
//   int KK = 0;
//   if(k==0){
//     k = round(N/20);
//   }
//   arma::vec sigma = sigma_input; int slen = sigma.size();
//   arma::vec allk = allk_input;   int klen = allk.size();
//   int kerlen = slen*klen;
//
//   //compute and sort Diff
//   arma::mat Diff_sort(N,N);
//   for (int j = 0; j < N; ++j){
//     Diff_sort.row(j) = sort(Diff.row(j));
//   }
//
//   //compute combined kernels
//   arma::cube D_Kernels = arma::zeros<arma::cube>(N,N,kerlen);
//   for (int l = 0; l < klen; ++l){
//     arma::mat TT(N, 1);
//     TT.col(0) = mean(Diff_sort.cols(1, allk(l)), 1);
//     arma::mat Sig = arma::zeros<arma::mat>(N,N);
//     Sig.each_col() = TT.col(0);
//     Sig = (Sig + Sig.t())/2;
//     for (int j = 0; j < slen; ++j){
//       arma::mat W = (normpdf(Diff, arma::zeros<arma::mat>(N,N), sigma(j)*Sig));
//       D_Kernels.slice(KK) = (W + W.t())/2;
//       KK = KK + 1;
//     }
//   }
//   return D_Kernels;
//   for (int i=0; i < kerlen; ++i){
//     arma::mat K = D_Kernels.slice(i);
//     arma::vec dinv = 1/sqrt(K.diag()+1);
//     arma::mat G = K % (dinv * dinv.t());
//     arma::mat G1 = arma::mat(N,N);
//     G1.each_col() = G.diag();
//     arma::mat D_Kernels_tmp = (G1 + G1.t() - 2*G)/2;
//     D_Kernels.slice(i) = D_Kernels_tmp;
//   }
//   arma::cube out = arma::zeros<arma::cube>(N,N,kerlen);
//   for (int i=0; i < kerlen; ++i){
//     arma::mat distX = D_Kernels.slice(i);
//     arma::mat distX1(N,N);
//     arma::umat idx(N,N);
//     for (int j=0; j < N; ++j){
//       distX1.row(j) = sort(distX.row(j));
//       idx.row(j)    = sort_index(distX.row(j)).t();
//     }
//
//     //knn
//     arma::mat A = arma::zeros<arma::mat>(N,N);
//     arma::mat di = distX1.cols(1,k+1);
//     arma::umat id = idx.cols(1,k+1);
//     arma::mat numerator(N, k+1);
//     numerator.each_col() = di.col(k);
//     numerator = numerator - di;
//     arma::vec temp = k * di.col(k) - sum(di.cols(0,k-1), 1);
//     arma::mat denominator(N, k+1);
//     denominator.each_col() = temp;
//     arma::mat div = numerator / denominator;
//     arma::mat a(N, k+1);
//     a.each_col() = arma::linspace(0, N-1, N);
//     for (int row = 0; row < N; ++row){
//       for (int col = 0; col < k+1; ++col){
//         A(a(row,col), id(row,col)) = div(row,col);
//       }
//     }
//     A.elem(find_nonfinite(A)).zeros();
//     out.slice(i) = (A + A.t())/2;
//   }
//   return out;
// }
//
// // [[Rcpp::export]]
// arma::cube dist_kernel_c(arma::mat X,
//                          arma::vec allk_input,
//                          arma::vec sigma_input,
//                          int k = 0){
//   int N = X.n_rows;
//   int KK = 0;
//   if(k==0){
//     k = round(N/20);
//   }
//
//   arma::vec sigma = sigma_input; int slen = sigma.size();
//   arma::vec allk = allk_input;   int klen = allk.size();
//   int kerlen = slen*klen;
//
//   //compute and sort Diff
//   arma::mat Diff = dist_c(X);
//   arma::mat Diff_sort(N,N);
//   for (int j = 0; j < N; ++j){
//     Diff_sort.row(j) = sort(Diff.row(j));
//   }
//   //compute combined kernels
//   arma::cube D_Kernels = arma::zeros<arma::cube>(N,N,kerlen);
//   for (int l = 0; l < klen; ++l){
//     arma::mat TT(N, 1);
//     TT.col(0) = mean(Diff_sort.cols(1, allk(l)), 1);
//     arma::mat Sig = arma::zeros<arma::mat>(N,N);
//     Sig.each_col() = TT.col(0);
//     Sig = (Sig + Sig.t())/2;
//     for (int j = 0; j < slen; ++j){
//       arma::mat W = (normpdf(Diff, arma::zeros<arma::mat>(N,N), sigma(j)*Sig));
//       D_Kernels.slice(KK) = (W + W.t())/2;
//       KK = KK + 1;
//     }
//   }
//   for (int i=0; i < kerlen; ++i){
//     arma::mat K = D_Kernels.slice(i);
//     arma::vec dinv = 1/sqrt(K.diag()+1);
//     arma::mat G = K % (dinv * dinv.t());
//     arma::mat G1 = arma::mat(N,N);
//     G1.each_col() = G.diag();
//     arma::mat D_Kernels_tmp = (G1 + G1.t() - 2*G)/2;
//     D_Kernels.slice(i) = D_Kernels_tmp;
//   }
//   arma::cube out = arma::zeros<arma::cube>(N,N,kerlen);
//   for (int i=0; i < kerlen; ++i){
//     arma::mat distX = D_Kernels.slice(i);
//     arma::mat distX1(N,N);
//     arma::umat idx(N,N);
//     for (int j=0; j < N; ++j){
//       distX1.row(j) = sort(distX.row(j));
//       idx.row(j)    = sort_index(distX.row(j)).t();
//     }
//
//     //knn
//     arma::mat A = arma::zeros<arma::mat>(N,N);
//     arma::mat di = distX1.cols(1,k+1);
//     arma::umat id = idx.cols(1,k+1);
//     arma::mat numerator(N, k+1);
//     numerator.each_col() = di.col(k);
//     numerator = numerator - di;
//     arma::vec temp = k * di.col(k) - sum(di.cols(0,k-1), 1);
//     arma::mat denominator(N, k+1);
//     denominator.each_col() = temp;
//     arma::mat div = numerator / denominator;
//     arma::mat a(N, k+1);
//     a.each_col() = arma::linspace(0, N-1, N);
//     for (int row = 0; row < N; ++row){
//       for (int col = 0; col < k+1; ++col){
//         A(a(row,col), id(row,col)) = div(row,col);
//       }
//     }
//     A.elem(find_nonfinite(A)).zeros();
//     out.slice(i) = (A + A.t())/2;
//   }
//   return out;
// }
//
// // [[Rcpp::export]]
// arma::cube rank_kernel_c(arma::mat X,
//                          arma::mat Diff,
//                          arma::vec allk_input,
//                          arma::vec sigma_input,
//                          int k = 0){
//   int N = X.n_rows;   int KK = 0;
//   if(k==0){
//     k = round(N/20);
//   }
//
//   arma::vec sigma = sigma_input; int slen = sigma.size();
//   arma::vec allk = allk_input;   int klen = allk.size();
//   int kerlen = slen*klen;
//
//   //compute and sort Diff
//   arma::mat Diff_sort(N,N);
//   for (int j = 0; j < N; ++j){
//     Diff_sort.row(j) = sort(Diff.row(j));
//   }
//
//   //compute combined kernels
//   arma::cube D_Kernels = arma::zeros<arma::cube>(N,N,kerlen);
//   for (int l = 0; l < klen; ++l){
//     arma::mat TT(N, 1);
//     TT.col(0) = mean(Diff_sort.cols(1, allk(l)), 1);
//     arma::mat Sig = arma::zeros<arma::mat>(N,N);
//     Sig.each_col() = TT.col(0);
//     Sig = (Sig + Sig.t())/2;
//     for (int j = 0; j < slen; ++j){
//       arma::mat W = (normpdf(Diff, arma::zeros<arma::mat>(N,N), sigma(j)*Sig));
//       D_Kernels.slice(KK) = (W + W.t())/2;
//       KK = KK + 1;
//     }
//   }
//   for (int i=0; i < kerlen; ++i){
//     arma::mat K = D_Kernels.slice(i);
//     arma::vec dinv = 1/sqrt(K.diag()+1);
//     arma::mat G = K % (dinv * dinv.t());
//     arma::mat G1 = arma::mat(N,N);
//     G1.each_col() = G.diag();
//     arma::mat D_Kernels_tmp = (G1 + G1.t() - 2*G)/2;
//     D_Kernels.slice(i) = D_Kernels_tmp;
//   }
//   arma::cube out = arma::zeros<arma::cube>(N,N,kerlen);
//   for (int i=0; i < kerlen; ++i){
//     arma::mat distX = D_Kernels.slice(i);
//     arma::mat distX1(N,N);
//     arma::umat idx(N,N);
//     for (int j=0; j < N; ++j){
//       distX1.row(j) = sort(distX.row(j));
//       idx.row(j)    = sort_index(distX.row(j)).t();
//     }
//
//     //knn
//     arma::mat A = arma::zeros<arma::mat>(N,N);
//     arma::mat di = distX1.cols(1,k+1);
//     arma::umat id = idx.cols(1,k+1);
//     arma::mat numerator(N, k+1);
//     numerator.each_col() = di.col(k);
//     numerator = numerator - di;
//     arma::vec temp = k * di.col(k) - sum(di.cols(0,k-1), 1);
//     arma::mat denominator(N, k+1);
//     denominator.each_col() = temp;
//     arma::mat div = numerator / denominator;
//     arma::mat a(N, k+1);
//     a.each_col() = arma::linspace(0, N-1, N);
//     for (int row = 0; row < N; ++row){
//       for (int col = 0; col < k+1; ++col){
//         A(a(row,col), id(row,col)) = div(row,col);
//       }
//     }
//     A.elem(find_nonfinite(A)).zeros();
//     out.slice(i) = (A + A.t())/2;
//   }
//   return out;
// }

// // [[Rcpp::export]]
// Rcpp::List sparse_scaledlasso_c(const arma::cube P,
//                                 const double tau,
//                                 const double gamma,
//                                 double mu=1e-3,
//                                 const double rho = 1.9,
//                                 const int max_iter=100,
//                                 const double eps = 1e-7,
//                                 const bool verbose = false){
//   int m = P.n_rows; int p = P.n_cols; int n = P.n_slices;
//   arma::vec funV     = arma::zeros<arma::vec>(max_iter);
//   arma::vec sigma    = arma::ones<arma::vec>(n);
//   arma::mat Y        = arma::zeros<arma::mat>(m,p);
//   arma::mat Q        = arma::zeros<arma::mat>(m,p);
//   arma::mat S        = sum(P, 2)/n;
//   arma::mat S_old    = arma::zeros<arma::mat>(m,p);
//   arma::vec L1_norm = arma::zeros<arma::vec>(max_iter);
//   arma::vec F_norm = arma::zeros<arma::vec>(max_iter);
//   arma::vec fro_norm = arma::zeros<arma::vec>(max_iter);
//   arma::vec noise    = arma::zeros<arma::vec>(max_iter);
//   int step = 0;
//   while(1){
//     if(step > max_iter-1){
//       cout << "reached max iteration \n";
//       break;
//     }
//     step += 1;
//     double max_inf_norm = norm(S-Q, "inf");
//     L1_norm(step-1) = norm(S, 1);
//     F_norm(step-1) = norm(S, "fro");
//     noise(step-1) = sum(sigma)/2;
//     for (int i=0; i < n; ++i){
//       fro_norm(step-1) += norm(Q-P.slice(i), "fro")/(sigma(i) * m * m);
//     }
//     funV(step-1) = fro_norm(step-1) + noise(step-1) +
//       L1_norm(step-1) + F_norm(step-1);
//     double relChg = norm(S-S_old, "fro")/max(1.0, norm(S_old, "fro"));
//     S_old = S;
//     if(verbose & (step % 5==0)){
//       cout << "iter" << step << ":"
//        // ": \n max_inf_norm = " << max_inf_norm <<
//           "\n relative change = " << relChg << "\n";
//        //     "\n mu = " << mu <<
//       //        "\n funV = " << funV(step-1) << "\n";
//     }
//     if(step > 1 & max_inf_norm < eps & relChg < eps){
//       break;
//     }
//
//     //Update Qi;
//     arma::mat M = mu*S + Y;
//     double phi = sum(1/(m*m*sigma)) + mu + gamma;
//     for (int i=0; i<n; ++i){
//       M += P.slice(i)/(m*m*sigma(i));
//     }
//     M = M/phi;
//     double C = tau/phi;
//     Q = armapmax(M-C, 0) + armapmin(M-C,0);
//
//     //update S;
//     S = nonnegASC_c(Q-Y/mu);
//
//     //update sigma
//     for (int i=0; i < n; ++i){
//       sigma(i) = norm(P.slice(i) - Q, "fro");
//       sigma(i) /= (m*m);
//       sigma(i) = 1/sigma(i);
//     }
//
//     //update Y
//     Y = Y + mu * (S - Q);
//
//     //update mu
//     mu = min(rho*mu, 1e+10);
//   }
//
//   arma::mat U; arma::vec s; arma::mat V;
//   svd(U,s,V,S);
//   arma::vec pi = V.col(0);
//   arma::vec Dist = pi / accu(pi) + 1e-10;
//   arma::vec invDist = 1/Dist;
//   arma::mat Distmat = diagmat(Dist);
//   arma::mat invDistmat = diagmat(invDist);
//   arma::mat sspi = sqrt(invDistmat);
//   arma::mat spi = sqrt(Distmat);
//   S = (spi * S * sspi + sspi * S.t() * spi)/2;
//
//   return Rcpp::List::create(
//     Rcpp::Named("S") = S,
//     Rcpp::Named("f") = funV,
//     Rcpp::Named("sigma") = sigma,
//     Rcpp::Named("L1_norm") = L1_norm,
//     Rcpp::Named("F_norm") = F_norm,
//     Rcpp::Named("fro_norm") = fro_norm,
//     Rcpp::Named("noise") = noise
//   );
// }
// // [[Rcpp::export]]
// arma::vec get_rank(const arma::vec X){
//   arma::vec ind = arma::conv_to<arma::vec>::from(sort_index(X));
//   int n = X.size();
//   arma::vec out = arma::zeros<arma::vec>(n);
//   for (int i=0; i<n; ++i){
//     arma::uvec tmp = find(ind==i);
//     out(i) = tmp(0);
//   }
//   return out;
// }






// // [[Rcpp::export]]
// arma::mat spearman_c(const arma::mat X){
//   int n = X.n_rows;
//   arma::mat out = arma::ones<arma::mat>(n,n);
//   for (int i=0; i<(n-1); ++i){
//     for (int j=(i+1); j<n; ++j){
//       arma::rowvec x1 = X.row(i);
//       arma::rowvec x2 = X.row(j);
//       arma::vec rx1 = get_rank(x1.t())+1;
//       arma::vec rx2 = get_rank(x2.t())+1;
//       out(i,j) = arma::conv_to<double>::from(cov(rx1, rx2))/(stddev(rx1) * stddev(rx2));
//       out(j,i) = out(i,j);
//     }
//   }
//   return out;
// }
