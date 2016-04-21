#include <Rcpp.h>
#include <algorithm>    // std::sort
#include <vector>       // std::vector
using namespace Rcpp;

// if we want to iterate over all possible filter thresholds (m), then the
// genefilter implementation of the optimal thresholding method is way too slow
// here a naive O(m^2) implementation is provided, also not terribly fast
// but at least sufficient for simulations with up to m \approx 40k

// [[Rcpp::export]]
Rcpp::IntegerVector ddhf_order(Rcpp::NumericVector sorted_p,
                Rcpp::IntegerVector filterorder,
                double alpha) {

    int m = sorted_p.size();
    Rcpp::IntegerVector rjs(m);

    std::vector<double> sortedp = Rcpp::as<std::vector<double> >(sorted_p);
    //std::sort(sortedp.begin(), sortedp.end());

    for (int i=0; i <= m-1; i++){

      if (i>0) {
        int deletepos = filterorder[i] - 1;
        sortedp.erase(sortedp.begin() + deletepos);

        for (int k=0; k < m; k ++){
          if (filterorder[k] > deletepos){
            --filterorder[k];
          }
        }
      }

      int m_filt = m - i;
      int j = m_filt-1;

      while ( (j>=0) && (sortedp[j] > alpha*(j+1)/m_filt) ) {
           j--;
      }

      rjs[i] = j+1;
    }

   //Rcpp::NumericVector tmp_vec = Rcpp::as<Rcpp::NumericVector>(pv_list[k])
   return rjs;
}
