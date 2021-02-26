#include "../../src/minconf.h"
#include "../../src/optimizer.h"
#include <RInside.h>
#include <Rcpp.h>
#include <iostream>
#include <random>
using namespace Rcpp;

int main() {
  std::random_device rd;
  long long seed = rd();
  std::vector<unsigned> alpha_list = {2, 1, 2};
  unsigned gamma = 3;
  std::vector<int> target = {-10, 0, 2, 0, -10, 0, 2, 0, -10};

  // MinConf mc(alpha_list, gamma, target, std::vector<int>(),
  // std::vector<int>(), seed);

  // now let's benchmark:
  //    BENCHMARK("MinConf") {
  // mc.optimize(5000, false, false); // will be finished after <<5k steps
  //    };

  RInside R;
  R.parseEvalQ(
      "dev <- spectre::run_optimization_min_conf(alpha_list = "
      "spectre::alpha_list, total_gamma = spectre::estimated_gamma, target = "
      "spectre::target_matrix, max_iterations = 8, seed = 1)");
  R.parseEvalQ(
      "x <- spectre:::calculate_solution_commonness_rcpp(dev$optimized_grid)");
  //  IntegerVector alpha_list = R["spectre::alpha_list"];
  //  const unsigned gamma = R["spectre::estimated_gamma"];
  //  IntegerMatrix target = R["spectre::target_matrix"];

  // auto res = calculate_solution_commonness_rcpp(mc.solution);

  //    IntegerMatrix sol = R.parseEval("current_solution");
  //    IntegerVector alpha_list = R.parseEval("alpha_list");
  //    int total_gamma = R.parseEval("estimated_gamma");
  //    IntegerMatrix target_matrix = R.parseEval("target_matrix");

  //    auto sol_com = calculate_solution_commonness_rcpp(sol);
  //    //R["result"]
  //    auto result = mh_optimizer(alpha_list, total_gamma, target_matrix);
  //    R.parseEvalQ("saveRDS(result, \"result.rds\"");
  ////    IntegerMatrix solution_matrix = R.parseEval("matrix(sample(0:1,
  /// 48,
  /// replace = TRUE), 6, 8)"); /    IntegerMatrix solution_commonness =
  /// R.parseEval("matrix(sample(0:8, 64, replace = TRUE), 8, 8)"); /    int
  /// site = 6;
  //    //IntegerMatrix solution = calculate_solution_commonness_site_rcpp(
  //    //            solution_matrix, solution_commonness, site);

  return 0;
}
