#include <Rcpp.h>
#include <vector>

using namespace Rcpp;

// [[Rcpp::export]]
SEXP merge_porc_cpp(
    SEXP porcSEXP,
    SEXP boundariesSEXP,
    SEXP finalBoundarySEXP
) {
    // Integer path
    if (TYPEOF(porcSEXP) == INTSXP) {

        IntegerVector porc(porcSEXP);
        IntegerVector boundaries = as<IntegerVector>(boundariesSEXP);

        const int final_boundary = as<int>(finalBoundarySEXP);

        std::vector<int> out;
        out.reserve(porc.size() + boundaries.size() + 2);

        // First sentinel
        out.push_back(-1);

        int i = 0;
        int j = 0;

        bool have_last = false;
        int last = 0;

        // Merge two sorted vectors while removing duplicates
        while (i < porc.size() && j < boundaries.size()) {

            int value;

            if (porc[i] < boundaries[j]) {
                value = porc[i++];
            } else if (boundaries[j] < porc[i]) {
                value = boundaries[j++];
            } else {
                value = porc[i];
                ++i;
                ++j;
            }

            if (!have_last || value != last) {
                out.push_back(value);
                last = value;
                have_last = true;
            }
        }

        while (i < porc.size()) {
            const int value = porc[i++];
            if (!have_last || value != last) {
                out.push_back(value);
                last = value;
                have_last = true;
            }
        }

        while (j < boundaries.size()) {
            const int value = boundaries[j++];
            if (!have_last || value != last) {
                out.push_back(value);
                last = value;
                have_last = true;
            }
        }

        // Final sentinel
        out.push_back(final_boundary + 1);

        return wrap(out);
    }

    // Double path
    NumericVector porc = as<NumericVector>(porcSEXP);
    NumericVector boundaries = as<NumericVector>(boundariesSEXP);

    const double final_boundary = as<double>(finalBoundarySEXP);

    std::vector<double> out;
    out.reserve(porc.size() + boundaries.size() + 2);

    // First sentinel
    out.push_back(-1.0);

    int i = 0;
    int j = 0;

    bool have_last = false;
    double last = 0.0;

    // Merge two sorted vectors while removing duplicates
    while (i < porc.size() && j < boundaries.size()) {

        double value;

        if (porc[i] < boundaries[j]) {
            value = porc[i++];
        } else if (boundaries[j] < porc[i]) {
            value = boundaries[j++];
        } else {
            value = porc[i];
            ++i;
            ++j;
        }

        if (!have_last || value != last) {
            out.push_back(value);
            last = value;
            have_last = true;
        }
    }

    while (i < porc.size()) {
        const double value = porc[i++];
        if (!have_last || value != last) {
            out.push_back(value);
            last = value;
            have_last = true;
        }
    }

    while (j < boundaries.size()) {
        const double value = boundaries[j++];
        if (!have_last || value != last) {
            out.push_back(value);
            last = value;
            have_last = true;
        }
    }

    // Final sentinel
    out.push_back(final_boundary + 1.0);

    return wrap(out);
}
