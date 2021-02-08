#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#define CATCH_CONFIG_ENABLE_BENCHMARKING // enable benchmarking
#include "catch.hpp"
#include "../../src/minconf.h"

TEST_CASE("MinConf") {

    std::random_device rd;
    long long seed = rd();
    std::vector<unsigned> alpha_list = {2, 1, 2};
    unsigned gamma = 3;
    std::vector<int> target = {  -10, 0, 2 ,
                                 0, -10, 0 ,
                                 2, 0, -10  };
    MinConf mc(alpha_list, gamma, target, seed);

    // now let's benchmark (will stop automatically when it found a solution after a few iterations)
    BENCHMARK("MinConf test1234567890test1234567890test1234567890test1234567890") { return mc.optimize(5000, false, false); };
    //    BENCHMARK("MinConf optimize 3 sites, 3 species, <<5k iterations") {
    //        return mc.optimize(5000, false, false); // will be finished after <<5k steps
    //    };


    //    const std::vector<unsigned> alpha_list_large = {15, 13, 13, 12, 21, 12, 15, 14, 14, 14, 13, 15, 12, 11, 13, 10, 14, 14, 18, 15,
    //                                                    13, 10, 12, 13, 12, 10, 12, 14, 15, 17, 11, 13, 10, 12, 11, 12, 9, 16, 14, 12, 10,
    //                                                    11, 12, 12, 14, 9, 13, 11, 14, 11, 15, 11, 14, 10, 12, 15, 14, 12, 19, 13, 11,
    //                                                    15, 11, 11, 12, 11, 12, 19, 12, 11, 13, 14, 11, 10, 11, 10, 13, 11, 11, 11, 10,
    //                                                    11, 13, 12, 13, 12, 12, 10, 17, 12, 11, 12, 13, 15, 13, 15, 14, 17, 13, 14};
    //    const unsigned gamma_large = 139;
    //    const std::vector<int> target_large = {-1,2,2,2,2,2,2,2,2,2,1,2,2,1,2,2,2,1,3,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,1,2,1,2,2,2,2,2,2,1,2,1,1,1,2,2,2,2,1,2,2,2,2,1,2,1,1,2,2,1,1,1,2,2,1,2,3,2,2,3,1,2,3,1,1,1,3,1,2,1,2,2,2,3,4,4,2,3,2,4,4,4,2,1,2,2,
    //                                           -1,-1,3,3,2,2,2,3,3,3,2,2,2,2,2,2,3,1,3,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,3,2,2,4,2,2,4,2,2,4,2,1,2,3,2,3,2,2,2,2,3,5,5,2,4,2,5,6,6,2,2,2,3,
    //                                           -1,-1,-1,2,1,1,1,2,2,2,1,1,1,1,1,1,1,1,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,1,1,2,1,1,2,1,1,2,1,1,1,2,1,1,1,1,1,1,2,2,3,1,2,1,3,3,3,1,1,1,2,
    //                                           -1,-1,-1,-1,2,2,2,2,2,2,1,1,2,1,2,2,2,1,3,2,2,2,1,1,2,2,2,2,2,2,2,2,2,2,1,2,1,2,2,2,2,1,2,1,2,1,1,1,2,2,2,2,1,2,2,2,2,1,2,1,1,2,2,1,1,1,2,2,1,1,3,2,1,3,1,1,3,1,1,1,2,1,2,1,2,2,2,2,3,4,2,3,2,4,4,4,2,1,2,2,
    //                                           -1,-1,-1,-1,-1,2,1,2,2,2,1,1,1,1,2,1,2,1,2,1,2,1,1,1,1,1,1,1,1,2,2,2,1,2,1,1,1,2,2,1,2,1,1,1,1,1,1,1,1,1,2,2,1,2,2,2,1,1,2,1,1,1,2,1,1,1,1,2,1,1,3,2,1,3,1,1,3,1,1,1,2,1,2,1,1,2,1,2,3,3,1,2,2,3,4,4,2,1,2,2,
    //                                           -1,-1,-1,-1,-1,-1,1,2,2,2,1,1,1,1,2,1,2,1,2,1,2,1,1,1,1,1,1,1,1,2,1,2,1,2,1,1,1,2,2,1,2,1,1,1,1,1,1,1,1,1,2,2,1,1,2,1,1,1,2,1,1,1,1,1,1,1,1,2,1,1,3,1,1,3,1,1,3,1,1,1,2,1,2,1,1,1,1,2,3,3,1,2,2,3,4,3,2,1,2,2,
    //                                           -1,-1,-1,-1,-1,-1,-1,2,2,2,1,1,1,1,1,1,2,1,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,1,1,1,1,1,1,1,1,1,1,1,1,1,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,1,1,2,1,1,2,1,1,3,1,1,1,2,1,2,1,1,1,1,2,3,3,1,2,1,3,3,3,1,1,1,2,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,3,3,2,2,2,2,2,2,3,1,3,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,3,2,2,4,2,2,4,2,2,4,2,1,2,3,2,3,2,2,2,2,3,5,5,2,3,2,5,6,5,2,1,2,3,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,3,2,2,2,2,2,2,3,2,4,2,2,2,2,2,2,2,2,2,2,3,2,2,2,2,2,2,2,3,2,2,2,2,2,2,2,2,2,2,2,2,2,3,2,2,3,2,2,2,2,2,1,2,2,2,2,2,2,3,2,2,4,2,2,4,2,2,5,2,1,2,3,2,3,2,2,2,2,3,5,5,2,4,2,5,6,6,2,2,2,3,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,2,2,2,2,2,2,3,2,4,2,2,2,2,2,2,2,2,2,2,3,2,2,2,2,2,2,2,3,2,2,2,2,2,2,2,2,2,2,2,2,2,3,2,2,3,2,2,2,2,2,1,2,2,2,2,2,2,3,2,2,4,2,2,4,2,2,5,2,1,2,3,2,3,2,2,2,2,3,5,5,2,4,2,5,6,6,2,2,2,3,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,2,2,2,2,2,2,1,3,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,1,2,2,2,2,1,2,2,1,2,4,2,2,4,2,2,4,2,1,2,3,1,2,2,2,2,2,3,4,5,2,3,2,5,5,5,2,1,2,3,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,2,2,3,2,3,2,4,2,3,2,2,2,2,2,2,2,2,3,2,3,2,3,2,2,2,3,3,2,3,2,2,2,2,2,2,2,2,2,3,3,2,2,3,2,2,2,3,2,2,2,2,2,2,2,2,3,2,2,5,2,2,5,2,2,5,2,2,2,4,2,3,2,2,2,2,4,5,6,2,4,3,6,7,6,3,2,3,3,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,2,2,2,1,3,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,1,2,1,2,2,2,2,2,2,1,2,1,1,1,2,2,2,2,1,2,2,2,2,1,2,1,1,2,2,1,1,1,2,2,1,2,3,2,2,3,1,1,3,1,1,1,3,1,2,1,2,2,2,3,4,4,2,3,2,4,4,4,2,1,2,2,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,2,2,2,1,3,2,2,2,1,1,2,2,2,2,2,2,2,2,2,2,1,2,1,2,2,2,2,1,2,1,2,1,1,1,2,2,2,2,1,2,2,2,2,1,2,1,1,2,2,1,1,1,2,2,1,2,3,2,1,3,1,1,3,1,1,1,2,1,2,1,2,2,2,2,3,4,2,3,2,4,4,4,2,1,2,2,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,2,2,1,3,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,1,1,1,2,2,2,2,1,2,2,2,2,2,2,1,1,2,2,2,2,1,2,2,1,2,3,2,2,3,2,2,4,1,1,1,3,1,2,1,2,2,2,3,4,4,2,3,2,4,5,5,2,1,2,2,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,2,1,3,2,2,2,1,1,2,2,2,2,2,2,2,2,2,2,1,2,1,2,2,1,2,1,2,1,2,1,1,1,2,2,2,2,1,2,2,2,2,1,2,1,1,2,2,1,1,1,2,2,1,1,3,2,1,3,1,1,3,1,1,1,2,1,2,1,2,2,2,2,3,4,1,3,2,4,4,4,2,1,2,2,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,3,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,1,1,1,2,2,2,2,1,2,2,2,2,2,2,1,1,2,2,2,2,1,2,2,1,2,3,2,2,3,2,2,4,1,1,2,3,1,2,1,2,2,2,3,4,4,2,3,2,4,5,5,2,1,2,2,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,4,2,2,2,2,2,2,2,2,2,2,3,2,2,2,2,2,2,2,3,2,2,2,2,2,2,2,2,2,2,2,2,2,3,2,2,3,2,2,2,2,2,1,2,2,2,2,2,2,3,2,2,4,2,2,4,2,2,5,2,1,2,3,2,3,2,2,2,2,3,5,5,2,4,2,5,6,6,3,2,2,3,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,2,2,2,1,1,2,2,2,2,2,2,2,2,2,2,1,2,1,2,2,2,2,1,2,1,2,1,1,1,2,2,2,2,1,2,2,2,2,1,2,1,1,2,2,1,1,1,2,2,1,2,3,2,1,3,1,1,3,1,1,1,3,1,2,1,2,2,2,3,3,4,2,3,2,4,4,4,2,1,2,2,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,1,2,2,2,2,1,2,2,1,2,4,2,2,4,2,2,4,2,1,2,3,1,2,2,2,2,2,3,4,4,2,3,2,4,5,5,2,1,2,2,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,1,2,1,1,1,2,2,2,2,1,2,2,2,2,1,2,1,1,2,2,2,2,1,2,2,1,2,3,2,2,3,2,2,4,1,1,1,3,1,2,1,2,2,2,3,4,4,2,3,2,4,5,4,2,1,2,2,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,1,2,2,2,2,1,2,3,1,2,4,2,2,4,2,2,4,2,1,2,3,1,2,2,2,2,2,3,4,5,2,3,2,5,5,5,2,1,2,3,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,2,2,2,2,2,2,3,2,2,2,2,2,2,2,3,2,2,2,2,2,2,2,2,2,2,2,2,2,3,2,2,3,2,2,2,2,2,1,2,2,2,2,2,2,3,2,2,4,2,2,4,2,2,5,2,1,2,3,2,3,2,2,2,2,3,5,5,2,4,2,5,6,6,3,2,2,3,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,2,2,2,2,2,3,2,2,2,2,2,2,2,3,3,2,2,2,2,2,2,2,2,2,2,2,3,3,2,2,3,2,2,2,2,2,1,2,2,2,2,2,2,3,2,2,5,2,2,4,2,2,5,2,2,2,4,2,3,2,2,2,2,4,5,6,2,4,2,5,6,6,3,2,3,3,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,2,1,1,1,1,1,1,1,2,1,1,1,1,1,1,1,1,1,1,1,1,1,2,1,1,2,1,1,1,1,1,1,1,1,1,1,1,1,2,1,1,2,1,1,2,1,1,3,1,1,1,2,1,2,1,1,1,1,2,3,3,1,2,1,3,3,3,2,1,1,2,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,2,2,2,2,2,2,2,2,2,2,2,3,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,1,2,2,2,2,1,2,3,1,2,4,2,2,4,2,2,4,2,1,2,3,1,2,2,2,2,2,3,4,5,2,3,2,5,5,5,2,1,2,3,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,1,1,1,2,2,2,2,1,2,2,2,2,2,2,1,1,2,2,2,2,1,2,2,1,2,3,2,2,3,2,2,4,1,1,1,3,1,2,1,2,2,2,3,4,4,2,3,2,4,5,5,2,1,2,2,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,2,3,2,2,2,2,2,2,2,3,2,2,2,2,2,2,2,2,2,2,2,2,2,3,2,2,3,2,2,2,2,2,1,2,2,2,2,2,2,3,2,2,4,2,2,4,2,2,5,2,1,2,4,2,3,2,2,2,2,4,5,5,2,4,2,5,6,6,3,2,2,3,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,1,1,2,2,2,2,1,2,2,2,2,2,2,1,1,2,2,2,2,1,2,2,1,2,4,2,2,3,2,2,4,1,1,2,3,1,2,1,2,2,2,3,4,4,2,3,2,4,5,5,2,1,2,2,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,1,2,2,2,2,1,2,2,1,2,4,2,2,4,2,2,4,2,1,2,3,1,2,2,2,2,2,3,4,5,2,3,2,4,5,5,2,1,2,3,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,2,2,2,2,2,2,3,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,1,2,2,2,2,1,2,3,1,2,4,2,2,4,2,2,4,2,1,2,3,1,2,2,2,2,2,3,4,5,2,3,2,5,5,5,2,1,2,3,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,1,2,2,1,1,1,1,1,1,1,1,1,1,1,2,2,1,1,2,1,1,1,1,1,1,1,1,1,1,1,1,2,1,1,3,1,1,2,1,1,3,1,1,1,2,1,2,1,1,1,1,2,3,3,1,2,1,3,3,3,2,1,2,2,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,2,1,2,1,2,2,2,2,1,2,1,2,1,1,1,2,2,2,2,1,2,2,2,2,1,2,1,1,2,2,1,1,1,2,2,1,2,3,2,1,3,1,1,3,1,1,1,3,1,2,1,2,2,2,3,4,4,2,3,2,4,4,4,2,1,2,2,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,2,2,2,3,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,3,2,2,4,2,2,4,2,2,4,2,1,2,3,2,3,2,2,2,2,3,5,5,2,4,2,5,6,6,2,2,2,3,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,2,2,3,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,3,2,2,4,2,2,4,2,2,5,2,1,2,3,2,3,2,2,2,2,3,5,5,2,4,2,5,6,6,2,2,2,3,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,2,1,1,1,1,1,1,1,1,1,1,1,1,1,2,1,1,2,1,1,1,1,1,1,1,1,1,1,1,1,2,1,1,2,1,1,2,1,1,3,1,1,1,2,1,2,1,1,1,1,2,3,3,1,2,1,3,3,3,2,1,1,2,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,2,2,1,2,1,1,1,2,1,1,1,2,2,2,2,1,2,2,2,2,1,2,1,1,1,2,1,1,1,2,2,1,1,3,2,1,3,1,1,3,1,1,1,2,1,2,1,2,2,2,2,3,4,1,2,2,4,4,4,2,1,2,2,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,1,2,2,2,2,1,2,2,1,2,4,2,2,4,2,2,4,2,1,2,3,1,2,2,2,2,2,3,4,5,2,3,2,5,5,5,2,1,2,3,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,1,2,2,2,2,1,2,2,1,2,4,2,2,4,2,2,4,2,1,2,3,1,2,2,2,2,2,3,4,4,2,3,2,4,5,5,2,1,2,2,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,3,2,2,4,2,2,4,2,2,4,2,1,2,3,2,3,2,2,2,2,3,5,5,2,4,2,5,6,6,2,2,2,3,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,1,1,2,1,1,2,1,1,2,1,1,1,2,1,2,1,1,1,1,2,3,3,1,2,1,3,3,3,1,1,1,2,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,2,1,2,1,1,1,2,2,2,2,1,2,2,2,2,1,2,1,1,2,2,1,1,1,2,2,1,2,3,2,1,3,1,1,3,1,1,1,2,1,2,1,2,2,2,2,3,4,2,3,2,4,4,4,2,1,2,2,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,2,1,1,1,2,2,2,2,1,2,2,2,2,1,2,1,1,2,2,1,1,1,2,2,1,2,3,2,2,3,1,1,3,1,1,1,3,1,2,1,2,2,2,3,4,4,2,3,2,4,4,4,2,1,2,2,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,2,1,1,1,2,2,2,2,1,2,2,2,2,1,2,1,1,2,2,1,1,1,2,2,1,2,3,2,2,3,1,1,3,1,1,1,3,1,2,1,2,2,2,3,4,4,2,3,2,4,4,4,2,1,2,2,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,2,2,2,2,2,3,3,2,2,3,2,2,2,3,2,2,2,2,2,2,2,2,3,2,2,5,2,2,5,2,2,5,2,2,2,4,2,3,2,2,2,2,4,5,6,2,4,3,6,7,6,3,2,3,3,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,2,2,2,2,2,2,2,2,2,2,2,2,2,1,1,2,2,2,2,1,2,3,1,2,4,2,2,4,2,2,4,2,1,2,3,1,2,2,2,2,2,3,4,5,2,3,2,5,5,5,2,1,2,3,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,3,2,2,4,2,2,4,2,2,5,2,1,2,3,2,3,2,2,2,2,3,5,5,2,4,2,5,6,6,2,2,2,3,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,3,1,2,4,2,2,4,2,2,4,2,1,2,3,1,2,2,2,2,2,3,5,5,2,3,2,5,6,5,2,1,2,3,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,2,2,3,2,2,3,2,2,2,2,2,1,2,2,2,2,2,2,3,2,2,4,2,2,4,2,2,5,2,1,2,4,2,3,2,2,2,2,4,5,5,2,4,2,5,6,6,3,2,2,3,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,2,2,1,2,2,2,2,1,2,1,1,2,2,1,2,1,2,2,1,2,3,2,2,3,2,2,4,1,1,1,3,1,2,1,2,2,2,3,4,4,2,3,2,4,5,4,2,1,2,2,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,3,2,2,3,2,2,2,2,2,1,2,2,2,2,2,2,3,2,2,5,2,2,4,2,2,5,2,2,2,4,2,3,2,2,2,2,4,5,6,2,4,2,5,6,6,3,2,3,3,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,2,2,2,2,1,2,1,1,2,2,1,1,1,2,2,1,2,3,2,1,3,1,1,3,1,1,1,2,1,2,1,2,2,2,2,3,4,2,3,2,4,4,4,2,1,2,2,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,2,2,2,2,2,2,1,1,2,2,2,2,1,2,2,1,2,3,2,2,3,2,2,4,1,1,2,3,1,2,1,2,2,2,3,4,4,2,3,2,4,5,5,2,1,2,2,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,2,2,2,2,2,1,1,2,2,2,2,1,2,3,1,2,4,2,2,4,2,2,4,2,1,2,3,1,2,2,2,2,2,3,4,5,2,3,2,5,5,5,2,1,2,3,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,2,2,2,2,2,1,2,2,2,2,2,2,3,2,2,4,2,2,4,2,2,4,2,1,2,3,2,3,2,2,2,2,3,5,5,2,3,2,5,6,5,2,1,2,3,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,2,2,2,2,1,2,2,2,2,2,2,3,2,2,5,2,2,4,2,2,5,2,2,2,4,2,3,2,2,2,2,4,5,6,2,4,2,5,6,6,3,2,3,3,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,2,2,1,1,2,2,2,2,1,2,2,1,2,3,2,2,3,2,2,4,1,1,2,3,1,2,1,2,2,2,3,4,4,2,3,2,4,5,5,2,1,2,2,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,2,2,1,2,2,2,2,2,2,3,2,2,4,2,2,4,2,2,5,2,1,2,3,2,3,2,2,2,2,3,5,5,2,4,2,5,6,6,3,2,2,3,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,2,2,2,3,2,2,2,3,3,2,2,5,3,2,5,2,2,6,2,2,2,4,2,3,2,3,3,3,4,6,6,2,4,3,6,7,7,3,2,3,3,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,2,2,2,2,2,2,3,2,2,5,2,2,4,2,2,5,2,2,2,4,2,3,2,2,2,2,4,5,6,2,4,2,5,6,6,3,2,3,3,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,2,2,2,2,2,2,3,2,2,4,2,2,4,2,2,5,2,1,2,3,2,3,2,2,2,2,3,5,5,2,4,2,5,6,6,2,2,2,3,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,2,2,2,2,2,3,2,2,5,2,2,4,2,2,5,2,2,2,4,2,3,2,2,2,2,4,5,6,2,4,2,5,6,6,3,2,3,3,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,2,2,2,2,3,1,2,4,2,2,4,2,2,4,2,1,2,3,1,2,2,2,2,2,3,5,5,2,3,2,5,6,5,2,1,2,3,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,2,2,2,3,1,2,4,2,2,4,2,2,4,2,1,2,3,1,2,2,2,2,2,3,5,5,2,3,2,5,6,5,2,1,2,3,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,2,2,3,2,2,4,2,2,4,2,2,5,2,1,2,3,2,3,2,2,2,2,3,5,5,2,4,2,5,6,6,2,2,2,3,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,2,3,1,2,4,2,2,4,2,2,4,2,1,2,3,1,2,2,2,2,2,3,5,5,2,3,2,5,6,5,2,1,2,3,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,3,2,2,4,2,2,4,2,2,4,2,1,2,3,2,3,2,2,2,2,3,5,5,2,4,2,5,6,6,2,2,2,3,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,2,2,5,3,2,5,2,2,6,2,2,2,4,2,3,2,3,3,3,4,6,6,2,4,3,6,7,7,3,2,3,3,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,2,4,2,2,4,2,2,4,2,1,2,3,1,3,2,2,2,2,3,5,5,2,3,2,5,6,5,2,1,2,3,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,4,2,2,4,2,2,4,2,1,2,3,1,2,2,2,2,2,3,5,5,2,3,2,5,6,5,2,1,2,3,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,2,2,4,2,2,5,2,1,2,3,2,3,2,2,2,2,3,5,5,2,4,2,5,6,6,2,2,2,3,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,2,4,2,2,5,2,1,2,3,2,3,2,2,2,2,3,5,5,2,4,2,5,6,6,2,2,2,3,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,4,2,2,4,2,1,2,3,1,2,2,2,2,2,3,5,5,2,3,2,5,6,5,2,1,2,3,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,2,2,4,2,1,2,3,1,3,2,2,2,2,3,5,5,2,3,2,5,6,5,2,2,2,3,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,2,4,2,1,2,3,1,2,2,2,2,2,3,5,5,2,3,2,5,6,5,2,1,2,3,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,4,2,1,2,3,1,2,2,2,2,2,3,4,5,2,3,2,5,5,5,2,1,2,3,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,2,1,2,3,2,3,2,2,2,2,3,5,5,2,4,2,5,6,6,2,2,2,3,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,2,3,1,3,2,2,2,2,3,5,5,2,3,2,5,6,5,2,1,2,3,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,2,3,1,2,2,2,2,2,3,5,5,2,3,2,5,6,5,2,1,2,3,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,3,2,3,2,2,2,2,3,5,5,2,3,2,5,6,6,2,2,2,3,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,2,2,2,2,2,3,4,5,2,3,2,5,5,5,2,1,2,3,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,3,2,2,2,2,3,5,5,2,3,2,5,6,5,2,2,2,3,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,2,2,2,2,3,5,5,2,4,2,5,6,6,2,2,2,3,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,2,2,2,3,5,5,2,4,2,5,6,6,2,2,2,3,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,2,2,3,5,5,2,4,2,5,6,6,2,2,2,3,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,2,3,4,5,2,3,2,5,6,5,2,1,2,3,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,3,5,5,2,4,2,5,6,6,2,2,2,3,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,4,5,2,3,2,5,5,5,2,1,2,3,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,6,2,4,3,6,7,7,3,2,3,3,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,2,4,2,5,6,6,2,2,2,3,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,3,2,5,5,5,2,1,2,3,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,2,5,6,6,2,2,2,3,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,5,6,6,2,2,2,3,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,7,6,3,2,3,3,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,6,2,2,2,3,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,3,2,3,3,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,2,2,3,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,3,3,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,3,
    //                                           -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1
    //                                          };
    //    MinConf mc_large(alpha_list_large, gamma_large, target_large, seed, std::vector<int>(), std::vector<int>(), -1);


    //    BENCHMARK("MinConf optimize 100 sites, 139 species, 5 iterations") {
    //        return mc_large.optimize(5, false, false);
    //    };

    //    BENCHMARK("Calculate commonness for a 100x100 solution") {
    //        return mc_large.calculate_commonness(mc_large.solution, alpha_list_large.size());
    //    };
}
