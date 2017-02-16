#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
NumericVector prHH(NumericMatrix phi_index,NumericMatrix phi,NumericVector G,
                          NumericVector M,int FF,int h) {
    
    int p = phi_index.ncol();
    int n = phi_index.nrow()/h;
    int houseindex = 0;
    NumericVector prEachComb(n);
    
    for (int i = 0; i < n; i++) {
        int houseindexend = houseindex + h;
        int g = G[0] - 1;
        double prodhouse = 1;
        int withinhouseindex = 0;
        for (int j = houseindex; j < houseindexend; j++) {
            int m = M[withinhouseindex] - 1;
            double phiprod = 1;
            for (int k = 0; k < p; k++){
                int phiindexjk = phi_index(j, k);
                phiprod *= phi(phiindexjk-1, g+(m*FF));
            }
            withinhouseindex += 1;
            prodhouse *= phiprod;
        }
        prEachComb[i] = prodhouse;
        houseindex += h;
    }
    return prEachComb;
}


