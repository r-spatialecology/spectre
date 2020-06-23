// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

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
// calculate_solution_commonness_site_rcpp
IntegerMatrix calculate_solution_commonness_site_rcpp(const IntegerMatrix solution_matrix, const IntegerMatrix solution_commonness, const int site);
RcppExport SEXP _spectre_calculate_solution_commonness_site_rcpp(SEXP solution_matrixSEXP, SEXP solution_commonnessSEXP, SEXP siteSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerMatrix >::type solution_matrix(solution_matrixSEXP);
    Rcpp::traits::input_parameter< const IntegerMatrix >::type solution_commonness(solution_commonnessSEXP);
    Rcpp::traits::input_parameter< const int >::type site(siteSEXP);
    rcpp_result_gen = Rcpp::wrap(calculate_solution_commonness_site_rcpp(solution_matrix, solution_commonness, site));
    return rcpp_result_gen;
END_RCPP
}
// calculate_solution_commonness_species_site_rcpp
IntegerMatrix calculate_solution_commonness_species_site_rcpp(const IntegerMatrix solution_matrix, const IntegerMatrix solution_commonness, const int site, const int species);
RcppExport SEXP _spectre_calculate_solution_commonness_species_site_rcpp(SEXP solution_matrixSEXP, SEXP solution_commonnessSEXP, SEXP siteSEXP, SEXP speciesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerMatrix >::type solution_matrix(solution_matrixSEXP);
    Rcpp::traits::input_parameter< const IntegerMatrix >::type solution_commonness(solution_commonnessSEXP);
    Rcpp::traits::input_parameter< const int >::type site(siteSEXP);
    Rcpp::traits::input_parameter< const int >::type species(speciesSEXP);
    rcpp_result_gen = Rcpp::wrap(calculate_solution_commonness_species_site_rcpp(solution_matrix, solution_commonness, site, species));
    return rcpp_result_gen;
END_RCPP
}
// update_solution_commonness_site_rcpp
void update_solution_commonness_site_rcpp(const IntegerMatrix solution_matrix, IntegerMatrix& solution_commonness, const unsigned site);
RcppExport SEXP _spectre_update_solution_commonness_site_rcpp(SEXP solution_matrixSEXP, SEXP solution_commonnessSEXP, SEXP siteSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerMatrix >::type solution_matrix(solution_matrixSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix& >::type solution_commonness(solution_commonnessSEXP);
    Rcpp::traits::input_parameter< const unsigned >::type site(siteSEXP);
    update_solution_commonness_site_rcpp(solution_matrix, solution_commonness, site);
    return R_NilValue;
END_RCPP
}
// optimizer_min_conf
List optimizer_min_conf(IntegerVector alpha_list, const unsigned total_gamma, IntegerMatrix target, const unsigned max_iterations, const double energy_threshold, unsigned long seed, bool verbose);
RcppExport SEXP _spectre_optimizer_min_conf(SEXP alpha_listSEXP, SEXP total_gammaSEXP, SEXP targetSEXP, SEXP max_iterationsSEXP, SEXP energy_thresholdSEXP, SEXP seedSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type alpha_list(alpha_listSEXP);
    Rcpp::traits::input_parameter< const unsigned >::type total_gamma(total_gammaSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type target(targetSEXP);
    Rcpp::traits::input_parameter< const unsigned >::type max_iterations(max_iterationsSEXP);
    Rcpp::traits::input_parameter< const double >::type energy_threshold(energy_thresholdSEXP);
    Rcpp::traits::input_parameter< unsigned long >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(optimizer_min_conf(alpha_list, total_gamma, target, max_iterations, energy_threshold, seed, verbose));
    return rcpp_result_gen;
END_RCPP
}
// optimizer_min_conf0
List optimizer_min_conf0(IntegerVector alpha_list, const unsigned total_gamma, IntegerMatrix target, IntegerMatrix fixed_species, IntegerMatrix partial_solution, const unsigned max_iterations, const double energy_threshold, unsigned long seed, bool verbose);
RcppExport SEXP _spectre_optimizer_min_conf0(SEXP alpha_listSEXP, SEXP total_gammaSEXP, SEXP targetSEXP, SEXP fixed_speciesSEXP, SEXP partial_solutionSEXP, SEXP max_iterationsSEXP, SEXP energy_thresholdSEXP, SEXP seedSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type alpha_list(alpha_listSEXP);
    Rcpp::traits::input_parameter< const unsigned >::type total_gamma(total_gammaSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type target(targetSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type fixed_species(fixed_speciesSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type partial_solution(partial_solutionSEXP);
    Rcpp::traits::input_parameter< const unsigned >::type max_iterations(max_iterationsSEXP);
    Rcpp::traits::input_parameter< const double >::type energy_threshold(energy_thresholdSEXP);
    Rcpp::traits::input_parameter< unsigned long >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(optimizer_min_conf0(alpha_list, total_gamma, target, fixed_species, partial_solution, max_iterations, energy_threshold, seed, verbose));
    return rcpp_result_gen;
END_RCPP
}
// optimizer_min_conf1
List optimizer_min_conf1(IntegerVector alpha_list, const unsigned total_gamma, IntegerMatrix target, IntegerMatrix fixed_species, IntegerMatrix partial_solution, const unsigned max_iterations, const double energy_threshold, unsigned long seed, bool verbose);
RcppExport SEXP _spectre_optimizer_min_conf1(SEXP alpha_listSEXP, SEXP total_gammaSEXP, SEXP targetSEXP, SEXP fixed_speciesSEXP, SEXP partial_solutionSEXP, SEXP max_iterationsSEXP, SEXP energy_thresholdSEXP, SEXP seedSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type alpha_list(alpha_listSEXP);
    Rcpp::traits::input_parameter< const unsigned >::type total_gamma(total_gammaSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type target(targetSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type fixed_species(fixed_speciesSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type partial_solution(partial_solutionSEXP);
    Rcpp::traits::input_parameter< const unsigned >::type max_iterations(max_iterationsSEXP);
    Rcpp::traits::input_parameter< const double >::type energy_threshold(energy_thresholdSEXP);
    Rcpp::traits::input_parameter< unsigned long >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(optimizer_min_conf1(alpha_list, total_gamma, target, fixed_species, partial_solution, max_iterations, energy_threshold, seed, verbose));
    return rcpp_result_gen;
END_RCPP
}
// optimizer_min_conf2
List optimizer_min_conf2(IntegerVector alpha_list, const unsigned total_gamma, IntegerMatrix target, IntegerMatrix fixed_species, IntegerMatrix partial_solution, const unsigned max_iterations, const double energy_threshold, unsigned long seed, bool verbose);
RcppExport SEXP _spectre_optimizer_min_conf2(SEXP alpha_listSEXP, SEXP total_gammaSEXP, SEXP targetSEXP, SEXP fixed_speciesSEXP, SEXP partial_solutionSEXP, SEXP max_iterationsSEXP, SEXP energy_thresholdSEXP, SEXP seedSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type alpha_list(alpha_listSEXP);
    Rcpp::traits::input_parameter< const unsigned >::type total_gamma(total_gammaSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type target(targetSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type fixed_species(fixed_speciesSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type partial_solution(partial_solutionSEXP);
    Rcpp::traits::input_parameter< const unsigned >::type max_iterations(max_iterationsSEXP);
    Rcpp::traits::input_parameter< const double >::type energy_threshold(energy_thresholdSEXP);
    Rcpp::traits::input_parameter< unsigned long >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(optimizer_min_conf2(alpha_list, total_gamma, target, fixed_species, partial_solution, max_iterations, energy_threshold, seed, verbose));
    return rcpp_result_gen;
END_RCPP
}
// optimizer_backtracking
List optimizer_backtracking(IntegerVector alpha_list, const unsigned total_gamma, IntegerMatrix target, const unsigned max_iterations, bool verbose);
RcppExport SEXP _spectre_optimizer_backtracking(SEXP alpha_listSEXP, SEXP total_gammaSEXP, SEXP targetSEXP, SEXP max_iterationsSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type alpha_list(alpha_listSEXP);
    Rcpp::traits::input_parameter< const unsigned >::type total_gamma(total_gammaSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type target(targetSEXP);
    Rcpp::traits::input_parameter< const unsigned >::type max_iterations(max_iterationsSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(optimizer_backtracking(alpha_list, total_gamma, target, max_iterations, verbose));
    return rcpp_result_gen;
END_RCPP
}
// calc_energy
double calc_energy(const IntegerMatrix solution_commonness, const IntegerMatrix solution_commonness_target);
RcppExport SEXP _spectre_calc_energy(SEXP solution_commonnessSEXP, SEXP solution_commonness_targetSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerMatrix >::type solution_commonness(solution_commonnessSEXP);
    Rcpp::traits::input_parameter< const IntegerMatrix >::type solution_commonness_target(solution_commonness_targetSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_energy(solution_commonness, solution_commonness_target));
    return rcpp_result_gen;
END_RCPP
}
// gen_init_solution
IntegerMatrix gen_init_solution(const IntegerVector alpha_list, const unsigned gamma_diversity);
RcppExport SEXP _spectre_gen_init_solution(SEXP alpha_listSEXP, SEXP gamma_diversitySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerVector >::type alpha_list(alpha_listSEXP);
    Rcpp::traits::input_parameter< const unsigned >::type gamma_diversity(gamma_diversitySEXP);
    rcpp_result_gen = Rcpp::wrap(gen_init_solution(alpha_list, gamma_diversity));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_sample
Rcpp::IntegerVector rcpp_sample(Rcpp::IntegerVector x, int size, bool replace);
RcppExport SEXP _spectre_rcpp_sample(SEXP xSEXP, SEXP sizeSEXP, SEXP replaceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type size(sizeSEXP);
    Rcpp::traits::input_parameter< bool >::type replace(replaceSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_sample(x, size, replace));
    return rcpp_result_gen;
END_RCPP
}
// get_swap_rows_rcpp
IntegerVector get_swap_rows_rcpp(const IntegerVector& v, bool cpp_index);
RcppExport SEXP _spectre_get_swap_rows_rcpp(SEXP vSEXP, SEXP cpp_indexSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerVector& >::type v(vSEXP);
    Rcpp::traits::input_parameter< bool >::type cpp_index(cpp_indexSEXP);
    rcpp_result_gen = Rcpp::wrap(get_swap_rows_rcpp(v, cpp_index));
    return rcpp_result_gen;
END_RCPP
}
// get_swap_rows_rcpp_bruteforce
IntegerVector get_swap_rows_rcpp_bruteforce(const IntegerVector& v, bool cpp_index);
RcppExport SEXP _spectre_get_swap_rows_rcpp_bruteforce(SEXP vSEXP, SEXP cpp_indexSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerVector& >::type v(vSEXP);
    Rcpp::traits::input_parameter< bool >::type cpp_index(cpp_indexSEXP);
    rcpp_result_gen = Rcpp::wrap(get_swap_rows_rcpp_bruteforce(v, cpp_index));
    return rcpp_result_gen;
END_RCPP
}
// which_not_vec
IntegerVector which_not_vec(const IntegerVector x, const IntegerVector y);
RcppExport SEXP _spectre_which_not_vec(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(which_not_vec(x, y));
    return rcpp_result_gen;
END_RCPP
}
// which_not
IntegerVector which_not(const IntegerVector x, const int y);
RcppExport SEXP _spectre_which_not(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< const int >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(which_not(x, y));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_spectre_calculate_solution_commonness_rcpp", (DL_FUNC) &_spectre_calculate_solution_commonness_rcpp, 1},
    {"_spectre_calculate_solution_commonness_site_rcpp", (DL_FUNC) &_spectre_calculate_solution_commonness_site_rcpp, 3},
    {"_spectre_calculate_solution_commonness_species_site_rcpp", (DL_FUNC) &_spectre_calculate_solution_commonness_species_site_rcpp, 4},
    {"_spectre_update_solution_commonness_site_rcpp", (DL_FUNC) &_spectre_update_solution_commonness_site_rcpp, 3},
    {"_spectre_optimizer_min_conf", (DL_FUNC) &_spectre_optimizer_min_conf, 7},
    {"_spectre_optimizer_min_conf0", (DL_FUNC) &_spectre_optimizer_min_conf0, 9},
    {"_spectre_optimizer_min_conf1", (DL_FUNC) &_spectre_optimizer_min_conf1, 9},
    {"_spectre_optimizer_min_conf2", (DL_FUNC) &_spectre_optimizer_min_conf2, 9},
    {"_spectre_optimizer_backtracking", (DL_FUNC) &_spectre_optimizer_backtracking, 5},
    {"_spectre_calc_energy", (DL_FUNC) &_spectre_calc_energy, 2},
    {"_spectre_gen_init_solution", (DL_FUNC) &_spectre_gen_init_solution, 2},
    {"_spectre_rcpp_sample", (DL_FUNC) &_spectre_rcpp_sample, 3},
    {"_spectre_get_swap_rows_rcpp", (DL_FUNC) &_spectre_get_swap_rows_rcpp, 2},
    {"_spectre_get_swap_rows_rcpp_bruteforce", (DL_FUNC) &_spectre_get_swap_rows_rcpp_bruteforce, 2},
    {"_spectre_which_not_vec", (DL_FUNC) &_spectre_which_not_vec, 2},
    {"_spectre_which_not", (DL_FUNC) &_spectre_which_not, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_spectre(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
