#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
NumericMatrix prMpost(NumericMatrix phi_index, NumericMatrix phi, 
                          NumericMatrix omega, NumericVector rep_G,
                          int FF, int SS) {
  
  int p = phi_index.ncol();
  int N =  phi_index.nrow();
  NumericMatrix prMpost(N, SS);
  
  for (int i = 0; i < N; i++) {
    int Gi = rep_G[i] - 1;
    NumericVector prMposti(SS);
    double updatesum = 0;
    
    for (int m = 0; m < SS; m++) {
      double phiprod = 1;
      for (int k = 0; k < p; k++) {
        int phiindexjk = phi_index(i, k);
        phiprod *= phi(phiindexjk-1, Gi+(m*FF));
      }
      prMposti[m] = omega(Gi, m)*phiprod;
      updatesum += prMposti[m];
    }
    
    for (int m = 0; m < SS; m++) {
      prMpost(i, m) = prMposti[m]/updatesum;
    }
    
  }
  return prMpost;
}


