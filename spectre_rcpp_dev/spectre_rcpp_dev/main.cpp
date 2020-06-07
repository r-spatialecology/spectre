#include <iostream>
#include <RInside.h>
#include <Rcpp.h>
#include "../../src/optimizer.h"
using namespace Rcpp;


int main()
{
    RInside R;
    R.parseEvalQ("alpha15 <- readRDS(\"~/sebastian@hanss.info/2_ECOMOD/1_Projects/spectre2/spectre/data/alpha15.rds\")");
    IntegerVector alpha_list = R["alpha15"];
    R.parseEvalQ("load(\"/home/bacchus/sebastian@hanss.info/2_ECOMOD/1_Projects/spectre2/spectre/data/estimated_gamma.rda\")");
    const unsigned total_gamma = R["estimated_gamma"];
    R.parseEvalQ("target15 <- readRDS(\"~/sebastian@hanss.info/2_ECOMOD/1_Projects/spectre2/spectre/data/target15.rds\")");
    IntegerMatrix target = R["target15"];

    R.parseEvalQ("load(\"~/sebastian@hanss.info/2_ECOMOD/1_Projects/spectre2/spectre/data/random_solution.rda\")");
    auto res = optimizer_min_conf2(alpha_list, 15, target, 15000);
//    IntegerMatrix sol = R.parseEval("current_solution");
//    IntegerVector alpha_list = R.parseEval("alpha_list");
//    int total_gamma = R.parseEval("estimated_gamma");
//    IntegerMatrix target_matrix = R.parseEval("target_matrix");

//    auto sol_com = calculate_solution_commonness_rcpp(sol);
//    //auto energy = calc_energy(target_matrix, sol_com);
//    //R["result"]
//    auto result = mh_optimizer(alpha_list, total_gamma, target_matrix);
//    R.parseEvalQ("saveRDS(result, \"result.rds\"");
////    IntegerMatrix solution_matrix = R.parseEval("matrix(sample(0:1, 48, replace = TRUE), 6, 8)");
////    IntegerMatrix solution_commonness = R.parseEval("matrix(sample(0:8, 64, replace = TRUE), 8, 8)");
////    int site = 6;
//    //IntegerMatrix solution = calculate_solution_commonness_site_rcpp(
//    //            solution_matrix, solution_commonness, site);


    return 0;
}
