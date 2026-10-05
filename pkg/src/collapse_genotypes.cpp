#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
CharacterVector collapse_genotypes_cpp(const IntegerMatrix& geno) {
    const R_xlen_t n_markers = geno.nrow();
    const R_xlen_t n_individuals = geno.ncol();

    CharacterVector result(n_individuals);

    for (R_xlen_t individual = 0; individual < n_individuals; ++individual) {
        std::string genotype;
        genotype.resize(n_markers);

        for (R_xlen_t marker = 0; marker < n_markers; ++marker) {
            genotype[marker] =
                static_cast<char>('0' + geno(marker, individual));
        }

        result[individual] = genotype;
    }

    return result;
}