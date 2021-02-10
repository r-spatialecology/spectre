// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// optimizer_min_conf
List optimizer_min_conf(const IntegerVector alpha_list, const unsigned total_gamma, const IntegerMatrix target, const unsigned max_iterations, const IntegerMatrix partial_solution, const IntegerMatrix fixed_species, const unsigned long seed, const bool verbose, const bool interruptible);
RcppExport SEXP _spectre_optimizer_min_conf(SEXP alpha_listSEXP, SEXP total_gammaSEXP, SEXP targetSEXP, SEXP max_iterationsSEXP, SEXP partial_solutionSEXP, SEXP fixed_speciesSEXP, SEXP seedSEXP, SEXP verboseSEXP, SEXP interruptibleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerVector >::type alpha_list(alpha_listSEXP);
    Rcpp::traits::input_parameter< const unsigned >::type total_gamma(total_gammaSEXP);
    Rcpp::traits::input_parameter< const IntegerMatrix >::type target(targetSEXP);
    Rcpp::traits::input_parameter< const unsigned >::type max_iterations(max_iterationsSEXP);
    Rcpp::traits::input_parameter< const IntegerMatrix >::type partial_solution(partial_solutionSEXP);
    Rcpp::traits::input_parameter< const IntegerMatrix >::type fixed_species(fixed_speciesSEXP);
    Rcpp::traits::input_parameter< const unsigned long >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< const bool >::type interruptible(interruptibleSEXP);
    rcpp_result_gen = Rcpp::wrap(optimizer_min_conf(alpha_list, total_gamma, target, max_iterations, partial_solution, fixed_species, seed, verbose, interruptible));
    return rcpp_result_gen;
END_RCPP
}
// calculate_solution_commonness_rcpp
IntegerMatrix calculate_solution_commonness_rcpp(const IntegerMatrix solution_matrix);
RcppExport SEXP _spectre_calculate_solution_commonness_rcpp(SEXP solution_matrixSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerMatrix >::type solution_matrix(solution_matrixSEXP);
    rcpp_result_gen = Rcpp::wrap(calculate_solution_commonness_rcpp(solution_matrix));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP run_testthat_tests();

static const R_CallMethodDef CallEntries[] = {
    {"_spectre_optimizer_min_conf", (DL_FUNC) &_spectre_optimizer_min_conf, 9},
    {"_spectre_calculate_solution_commonness_rcpp", (DL_FUNC) &_spectre_calculate_solution_commonness_rcpp, 1},
    {"run_testthat_tests", (DL_FUNC) &run_testthat_tests, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_spectre(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
