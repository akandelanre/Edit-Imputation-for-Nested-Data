#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
NumericMatrix prGpost(NumericMatrix phi_index, NumericMatrix lambda_index, 
                          NumericMatrix phi, NumericMatrix lambda, NumericMatrix omega,
                          NumericVector pii, int FF, int SS, NumericVector cn_i) {
  
  int p = phi_index.ncol();
  int n = lambda_index.nrow(), q = lambda_index.ncol();
  int houseindex = 0;
  NumericMatrix prGpost(n, FF);
  
  for (int i = 0; i < n; i++) {
    NumericVector prGposti(FF);
    double updatesum = 0;
    int houseindexend = houseindex + cn_i[i];
    
    for (int g = 0; g < FF; g++) {
      double prodhouse = 1;
      for (int j = houseindex; j < houseindexend; j++) {
        double sumomegaprodphi = 0;
        for (int m = 0; m < SS; m++) {
          double phiprod = 1;
          for (int k = 0; k < p; k++){
            int phiindexjk = phi_index(j, k);
            phiprod *= phi(phiindexjk-1, g+(m*FF));
          }
          sumomegaprodphi += phiprod*omega(g, m);
        }
        prodhouse *= sumomegaprodphi;
      }
      
      double prodlambda = 1;
      for (int kq = 0; kq < q; kq++) {
        int lambdaindexjk = lambda_index(i, kq);
        prodlambda *= lambda(lambdaindexjk-1, g);
      }
      
      prGposti[g] = pii[g]*prodlambda*prodhouse;
      updatesum += prGposti[g];
    }
    
    for (int g = 0; g < FF; g++) {
      prGpost(i, g) = prGposti[g]/updatesum;
    }
    
    houseindex += cn_i[i];
  }
  return prGpost;
}


