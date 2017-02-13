#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
NumericMatrix prZpost(NumericMatrix epsilon_indiv_index, NumericMatrix epsilon_house_index, 
                          NumericMatrix epsilon_indiv, NumericMatrix epsilon_house,
                          NumericVector zeta, int TT, NumericVector n_i) {
  
  int p = epsilon_indiv_index.ncol();
  int n = epsilon_house_index.nrow(), q = epsilon_house_index.ncol();
  int houseindex = 0;
  NumericMatrix prZpost(n, TT);
  
  for (int i = 0; i < n; i++) {
    NumericVector prZposti(TT);
    double updatesum = 0;
    int houseindexend = houseindex + n_i[i];
    
    for (int t = 0; t < TT; t++) {
      double prodhouse = 1;
      for (int j = houseindex; j < houseindexend; j++) {
          double epsilon_indivprod = 1;
          for (int k = 0; k < p; k++){
            int epsilon_indivindexjk = epsilon_indiv_index(j, k);
            epsilon_indivprod *= epsilon_indiv(epsilon_indivindexjk-1, t);
          }
        prodhouse *= epsilon_indivprod;
      }
      
      double epsilon_houseprod = 1;
      for (int kq = 0; kq < q; kq++) {
        int epsilon_houseindexjk = epsilon_house_index(i, kq);
        epsilon_houseprod *= epsilon_house(epsilon_houseindexjk-1, t);
      }
      
      prZposti[t] = zeta[t]*epsilon_houseprod*prodhouse;
      updatesum += prZposti[t];
    }
    
    for (int t = 0; t < TT; t++) {
      prZpost(i, t) = prZposti[t]/updatesum;
    }
    
    houseindex += n_i[i];
  }
  return prZpost;
}


