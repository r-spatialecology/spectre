// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// calculate_solution_commonness_site_rcpp
IntegerMatrix calculate_solution_commonness_site_rcpp(IntegerMatrix solution_matrix, IntegerMatrix solution_commonness, int site);
RcppExport SEXP _OurFOAM_calculate_solution_commonness_site_rcpp(SEXP solution_matrixSEXP, SEXP solution_commonnessSEXP, SEXP siteSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type solution_matrix(solution_matrixSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type solution_commonness(solution_commonnessSEXP);
    Rcpp::traits::input_parameter< int >::type site(siteSEXP);
    rcpp_result_gen = Rcpp::wrap(calculate_solution_commonness_site_rcpp(solution_matrix, solution_commonness, site));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_OurFOAM_calculate_solution_commonness_site_rcpp", (DL_FUNC) &_OurFOAM_calculate_solution_commonness_site_rcpp, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_OurFOAM(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
