#include <Rcpp.h>
#include <vector>
#include <algorithm>

using namespace Rcpp;

// [[Rcpp::export]]
List breeding_loop_cpp(
    SEXP porcSEXP,
    IntegerVector start1,
    IntegerVector start2,
    IntegerVector end1,
    IntegerVector end2,
    SEXP parent_rec1SEXP,
    SEXP parent_rec2SEXP,
    IntegerVector parent_ori1,
    IntegerVector parent_ori2,
    LogicalVector mut_bool,
    List store_mut,
    IntegerVector parent_mut1,
    IntegerVector parent_mut2,
    double recombination_precision,
    double tmp1 = -2
) {

    std::vector<int> new_mut;
    std::vector<int> new_origin;

    int activ = 1;

    if (recombination_precision == 1.0) {

        // Double path
        NumericVector porc = as<NumericVector>(porcSEXP);
        NumericVector parent_rec1 = as<NumericVector>(parent_rec1SEXP);
        NumericVector parent_rec2 = as<NumericVector>(parent_rec2SEXP);

        std::vector<double> new_poc;
        const int n_intervals = porc.size() - 1;

        for (int index = 0; index < n_intervals; ++index) {

            int tmps;
            int tmpe;

            const NumericVector* rec_activ;
            const IntegerVector* ori_activ;

            if (activ == 1) {
                tmps = start1[index];
                tmpe = end1[index + 1];
                rec_activ = &parent_rec1;
                ori_activ = &parent_ori1;
            } else {
                tmps = start2[index];
                tmpe = end2[index + 1];
                rec_activ = &parent_rec2;
                ori_activ = &parent_ori2;
            }

            const double tmp8 = porc[index + 1];
            const bool save1 = (tmp1 != tmp8);

            if (save1) {
                tmp1 = tmp8;

                if (tmps != tmpe) {
                    for (int j = tmps; j < tmpe; ++j) {
                        new_poc.push_back((*rec_activ)[j]);
                    }
                }

                new_poc.push_back(tmp1);
            }

            if (!mut_bool[activ - 1]) {

                NumericVector sm = store_mut[activ - 1];

                const IntegerVector& parent_mut =
                    (activ == 1) ? parent_mut1 : parent_mut2;

                for (int j = 0; j < sm.size(); ++j) {
                    if (sm[j] < porc[index + 1] &&
                        sm[j] > porc[index]) {
                        new_mut.push_back(parent_mut[j]);
                    }
                }
            }

            int start = tmps;
            int tmp3 = std::max(start, tmpe);

            if (start == 0) {
                start = 1;
            }

            tmp3 = std::min(
                static_cast<int>(ori_activ->size()),
                tmp3
            );

            if (porc[index + 1] > 0 && save1) {
                if (start <= tmp3) {
                    for (int j = start - 1; j < tmp3; ++j) {
                        new_origin.push_back((*ori_activ)[j]);
                    }
                }
            }

            activ = 3 - activ;
        }

        // Equivalent to:
        // new.poc <- new.poc[-length(new.poc)]
        if (!new_poc.empty()) {
            new_poc.pop_back();
        }

        return List::create(
            _["new.poc"] = wrap(new_poc),
            _["new.mut"] = wrap(new_mut),
            _["new.origin"] = wrap(new_origin),
            _["tmp1"] = tmp1
        );

    } else {

        // Integer path
        IntegerVector porc = as<IntegerVector>(porcSEXP);
        IntegerVector parent_rec1 = as<IntegerVector>(parent_rec1SEXP);
        IntegerVector parent_rec2 = as<IntegerVector>(parent_rec2SEXP);

        std::vector<int> new_poc;
        const int n_intervals = porc.size() - 1;

        // tmp1 originates as -2 and then is always assigned integer breakpoints here.
        int tmp1_int = static_cast<int>(tmp1);

        for (int index = 0; index < n_intervals; ++index) {

            int tmps;
            int tmpe;

            const IntegerVector* rec_activ;
            const IntegerVector* ori_activ;

            if (activ == 1) {
                tmps = start1[index];
                tmpe = end1[index + 1];
                rec_activ = &parent_rec1;
                ori_activ = &parent_ori1;
            } else {
                tmps = start2[index];
                tmpe = end2[index + 1];
                rec_activ = &parent_rec2;
                ori_activ = &parent_ori2;
            }

            const int tmp8 = porc[index + 1];
            const bool save1 = (tmp1_int != tmp8);

            if (save1) {
                tmp1_int = tmp8;

                if (tmps != tmpe) {
                    for (int j = tmps; j < tmpe; ++j) {
                        new_poc.push_back((*rec_activ)[j]);
                    }
                }

                new_poc.push_back(tmp1_int);
            }

            if (!mut_bool[activ - 1]) {

                NumericVector sm = store_mut[activ - 1];

                const IntegerVector& parent_mut =
                    (activ == 1) ? parent_mut1 : parent_mut2;

                // store_mut contains SNP positions from population$info$snp.position,
                // which remain numeric positions in the current R implementation.
                for (int j = 0; j < sm.size(); ++j) {
                    if (sm[j] < porc[index + 1] &&
                        sm[j] > porc[index]) {
                        new_mut.push_back(parent_mut[j]);
                    }
                }
            }

            int start = tmps;
            int tmp3 = std::max(start, tmpe);

            if (start == 0) {
                start = 1;
            }

            tmp3 = std::min(
                static_cast<int>(ori_activ->size()),
                tmp3
            );

            if (porc[index + 1] > 0 && save1) {
                if (start <= tmp3) {
                    for (int j = start - 1; j < tmp3; ++j) {
                        new_origin.push_back((*ori_activ)[j]);
                    }
                }
            }

            activ = 3 - activ;
        }

        // Equivalent to:
        // new.poc <- new.poc[-length(new.poc)]
        if (!new_poc.empty()) {
            new_poc.pop_back();
        }

        return List::create(
            _["new.poc"] = wrap(new_poc),
            _["new.mut"] = wrap(new_mut),
            _["new.origin"] = wrap(new_origin),
            _["tmp1"] = tmp1_int
        );
    }
}
