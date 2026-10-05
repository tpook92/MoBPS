#include <Rcpp.h>
using namespace Rcpp;

//' Compare an integer matrix to a value and return a double matrix
//'
//' For each element of \code{x}, returns 1.0 if it equals \code{value}
//' and 0.0 otherwise.
//'
//' @param x Integer matrix.
//' @param value Integer value to compare against.
//' @return A numeric (double) matrix containing only 0.0 and 1.0.
//' @export
// [[Rcpp::export]]
NumericMatrix equal_double(IntegerMatrix x, int value) {
    int nr = x.nrow();
    int nc = x.ncol();

    NumericMatrix out(nr, nc);

    R_xlen_t n = x.size();
    const int* px = x.begin();
    double* po = out.begin();

    for (R_xlen_t i = 0; i < n; ++i) {
        po[i] = (px[i] == value);
    }

    return out;
}
