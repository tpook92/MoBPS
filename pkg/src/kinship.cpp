
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix kinship_cpp(IntegerMatrix tmp1,
                          IntegerMatrix tmp2) {

  int m = tmp1.nrow();
  int n = tmp1.ncol();

  NumericMatrix K(n, n);

  for (int r = 0; r < m; r++) {
    for (int i = 0; i < n; i++) {
      int t1_i = tmp1(r, i);
      int t2_i = tmp2(r, i);

      for (int j = i; j < n; j++) {

        int val = 0;

        if (tmp1(r, j) == t1_i) val++;
        if (tmp2(r, j) == t1_i) val++;
        if (tmp1(r, j) == t2_i) val++;
        if (tmp2(r, j) == t2_i) val++;

        K(i, j) += val;
        if (i != j) K(j, i) += val;
      }
    }
  }

  return K;
}

