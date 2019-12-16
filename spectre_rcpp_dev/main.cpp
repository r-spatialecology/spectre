#include <iostream>
#include <RInside.h>
#include <Rcpp.h>
using namespace Rcpp;

#include "../src/optimizer.h"
#include "../src/calculate_solution_commonness.h"


int main()
{
    RInside R;
    R.parseEvalQ("library(dplyr)");
    R.parseEvalQ("testdata <- matrix(c(0,0,0,1,0,0,0,1,0,0,1,0,0,0,1,0,0,1,0,0,1,1,1,1,1), nrow = 5, ncol = 5)");
   // R.parseEvalQ("saveRDS(testdata, \"result.rds\")");
    IntegerMatrix testdata = R.parseEval("testdata");
    NumericVector species_prop = R.parseEval("c(2/5, 1/5, 3/5, 2/5, 2/5)");
    IntegerVector alpha_list = R.parseEval("testdata %>% as_tibble() %>% summarise_all(funs(sum)) %>% slice(1) %>% unlist(., use.names = FALSE)");
    unsigned total_gamma = R.parseEval("testdata %>% as_tibble() %>% filter_all(any_vars(sum(.) != 0)) %>% nrow()");
    IntegerMatrix target_matrix = calculate_solution_commonness_rcpp(testdata);
    List test = mh_optimizer(alpha_list, total_gamma, target_matrix, species_prop, 20000, 0.0, .01);

    //    auto sol_com = calculate_solution_commonness_rcpp(sol);
    //    //auto energy = calc_energy(target_matrix, sol_com);
    R["result"] = test;
    //    auto result = mh_optimizer(alpha_list, total_gamma, target_matrix);
    R.parseEvalQ("saveRDS(result, \"./result.rds\")");
    ////    IntegerMatrix solution_matrix = R.parseEval("matrix(sample(0:1, 48, replace = TRUE), 6, 8)");
    ////    IntegerMatrix solution_commonness = R.parseEval("matrix(sample(0:8, 64, replace = TRUE), 8, 8)");
    ////    int site = 6;
    //    //IntegerMatrix solution = calculate_solution_commonness_site_rcpp(
    //    //            solution_matrix, solution_commonness, site);


    return 0;
}
