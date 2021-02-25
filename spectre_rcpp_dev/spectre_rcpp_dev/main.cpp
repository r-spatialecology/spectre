#include <iostream>
#include <RInside.h>
#include <Rcpp.h>
#include <random>
#include "../../src/minconf.h"
//#include "../../src/optimizer.h"
//using namespace Rcpp;


int main()
{
    std::random_device rd;
    long long seed = rd();
//    std::vector<unsigned> alpha_list = {2, 1, 2};
//    unsigned gamma = 3;
//    std::vector<int> target = {  -10, 0, 2 ,
//                                 0, -10, 0 ,
//                                 2, 0, -10  };

    std::vector<unsigned> alpha_list = {2, 1, 2};
    unsigned gamma = 3;
    std::vector<int> target = {-10, 0, 0, 0, -10, 0, 0, 0, -10};
    MinConf mc(alpha_list, gamma, target, std::vector<int>(), std::vector<int>(), 0, -10);
    mc.target[0][1] = 0.667;
    mc.target[0][2] = 0;
    mc.target[1][2] = 0;
    mc.solution = {{0, 1, 1}, {0, 1, 0}, {0, 1, 1}};
    mc.calc_min_conflict_species(1);
    // now let's benchmark:
    //    BENCHMARK("MinConf") {
   // mc.optimize(5000, false, false); // will be finished after <<5k steps
    //    };


   /* RInside R;
    R.parseEvalQ("alpha15 <- readRDS(\"~/sebastian@hanss.info/2_ECOMOD/1_Projects/spectre2/spectre/data/alpha15.rds\")");
    IntegerVector alpha_list = R["alpha15"];
    R.parseEvalQ("load(\"/home/bacchus/sebastian@hanss.info/2_ECOMOD/1_Projects/spectre2/spectre/data/estimated_gamma.rda\")");
    const unsigned total_gamma = R["estimated_gamma"];
    R.parseEvalQ("target15 <- readRDS(\"~/sebastian@hanss.info/2_ECOMOD/1_Projects/spectre2/spectre/data/target15.rds\")");
    IntegerMatrix target = R["target15"];

    R.parseEvalQ("load(\"~/sebastian@hanss.info/2_ECOMOD/1_Projects/spectre2/spectre/data/random_solution.rda\")");
    auto res = optimizer_min_conf(alpha_list, total_gamma, target, IntegerMatrix(), IntegerMatrix(),
                                  50, 0, false, false);
    /*/
    //    IntegerMatrix sol = R.parseEval("current_solution");
    //    IntegerVector alpha_list = R.parseEval("alpha_list");
    //    int total_gamma = R.parseEval("estimated_gamma");
    //    IntegerMatrix target_matrix = R.parseEval("target_matrix");

    //    auto sol_com = calculate_solution_commonness_rcpp(sol);
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
