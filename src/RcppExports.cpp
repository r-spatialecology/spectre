// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// optimizer_min_conf
List optimizer_min_conf(IntegerVector alpha_list, const unsigned total_gamma, IntegerMatrix target, IntegerMatrix fixed_species, IntegerMatrix partial_solution, const unsigned max_iterations, unsigned long seed, bool verbose, bool interruptible);
RcppExport SEXP _spectre_optimizer_min_conf(SEXP alpha_listSEXP, SEXP total_gammaSEXP, SEXP targetSEXP, SEXP fixed_speciesSEXP, SEXP partial_solutionSEXP, SEXP max_iterationsSEXP, SEXP seedSEXP, SEXP verboseSEXP, SEXP interruptibleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type alpha_list(alpha_listSEXP);
    Rcpp::traits::input_parameter< const unsigned >::type total_gamma(total_gammaSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type target(targetSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type fixed_species(fixed_speciesSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type partial_solution(partial_solutionSEXP);
    Rcpp::traits::input_parameter< const unsigned >::type max_iterations(max_iterationsSEXP);
    Rcpp::traits::input_parameter< unsigned long >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< bool >::type interruptible(interruptibleSEXP);
    rcpp_result_gen = Rcpp::wrap(optimizer_min_conf(alpha_list, total_gamma, target, fixed_species, partial_solution, max_iterations, seed, verbose, interruptible));
    return rcpp_result_gen;
END_RCPP
}
// calculate_solution_commonness_rcpp
IntegerMatrix calculate_solution_commonness_rcpp(IntegerMatrix solution_matrix);
RcppExport SEXP _spectre_calculate_solution_commonness_rcpp(SEXP solution_matrixSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type solution_matrix(solution_matrixSEXP);
    rcpp_result_gen = Rcpp::wrap(calculate_solution_commonness_rcpp(solution_matrix));
    return rcpp_result_gen;
END_RCPP
}
// calc_error_random_solution
unsigned calc_error_random_solution(const unsigned n, IntegerVector alpha_list, const unsigned total_gamma, IntegerMatrix target, unsigned long seed);
RcppExport SEXP _spectre_calc_error_random_solution(SEXP nSEXP, SEXP alpha_listSEXP, SEXP total_gammaSEXP, SEXP targetSEXP, SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const unsigned >::type n(nSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type alpha_list(alpha_listSEXP);
    Rcpp::traits::input_parameter< const unsigned >::type total_gamma(total_gammaSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type target(targetSEXP);
    Rcpp::traits::input_parameter< unsigned long >::type seed(seedSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_error_random_solution(n, alpha_list, total_gamma, target, seed));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP run_testthat_tests();

static const R_CallMethodDef CallEntries[] = {
    {"_spectre_optimizer_min_conf", (DL_FUNC) &_spectre_optimizer_min_conf, 9},
    {"_spectre_calculate_solution_commonness_rcpp", (DL_FUNC) &_spectre_calculate_solution_commonness_rcpp, 1},
    {"_spectre_calc_error_random_solution", (DL_FUNC) &_spectre_calc_error_random_solution, 5},
    {"run_testthat_tests", (DL_FUNC) &run_testthat_tests, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_spectre(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
