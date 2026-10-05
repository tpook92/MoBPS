#include <Rcpp.h>

using namespace Rcpp;

// Process one sorted parent breakpoint vector against sorted porc.
// For each porc[i]:
//   start[i] = findInterval(porc[i], parent_rec)
//            = number of parent_rec values <= porc[i]
//   end[i]   = findInterval(porc[i], parent_rec, left.open = TRUE)
//            = number of parent_rec values < porc[i]
//
// Both input vectors must be sorted non-decreasingly.
template <typename T>
inline void interval_sweep(
    const T* porc,
    const int n_porc,
    const T* parent_rec,
    const int n_parent,
    IntegerVector& start,
    IntegerVector& end
) {
    int p_lt = 0;
    int p_le = 0;

    for (int i = 0; i < n_porc; ++i) {
        const T x = porc[i];

        while (p_lt < n_parent && parent_rec[p_lt] < x) {
            ++p_lt;
        }

        if (p_le < p_lt) {
            p_le = p_lt;
        }

        while (p_le < n_parent && parent_rec[p_le] <= x) {
            ++p_le;
        }

        end[i] = p_lt;
        start[i] = p_le;
    }
}

// [[Rcpp::export]]
List find_intervals_joint_cpp(
    SEXP porcSEXP,
    SEXP parent_rec1SEXP,
    SEXP parent_rec2SEXP
) {
    const int n_porc = Rf_length(porcSEXP);

    IntegerVector start1(n_porc);
    IntegerVector start2(n_porc);
    IntegerVector end1(n_porc);
    IntegerVector end2(n_porc);

    if (TYPEOF(porcSEXP) == INTSXP) {

        IntegerVector porc(porcSEXP);
        IntegerVector parent_rec1 = as<IntegerVector>(parent_rec1SEXP);
        IntegerVector parent_rec2 = as<IntegerVector>(parent_rec2SEXP);

        interval_sweep<int>(
            INTEGER(porc), porc.size(),
            INTEGER(parent_rec1), parent_rec1.size(),
            start1, end1
        );

        interval_sweep<int>(
            INTEGER(porc), porc.size(),
            INTEGER(parent_rec2), parent_rec2.size(),
            start2, end2
        );

    } else {

        NumericVector porc = as<NumericVector>(porcSEXP);
        NumericVector parent_rec1 = as<NumericVector>(parent_rec1SEXP);
        NumericVector parent_rec2 = as<NumericVector>(parent_rec2SEXP);

        interval_sweep<double>(
            REAL(porc), porc.size(),
            REAL(parent_rec1), parent_rec1.size(),
            start1, end1
        );

        interval_sweep<double>(
            REAL(porc), porc.size(),
            REAL(parent_rec2), parent_rec2.size(),
            start2, end2
        );
    }

    return List::create(
        _["start1"] = start1,
        _["start2"] = start2,
        _["end1"]   = end1,
        _["end2"]   = end2
    );
}
