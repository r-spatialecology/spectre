#include "catch.h"
#include "test-minconf.h"
#include <RInside.h>

void TestMinConf::test_add_species_min_conf(unsigned site, const std::vector<std::vector<int> > &target)
{
    add_species_min_conf(site, target);
}

void TestMinConf::test_remove_species_max_conf(unsigned site, const std::vector<std::vector<int> > &target)
{
    remove_species_max_conf(site, target);
}

std::vector<unsigned> TestMinConf::test_calc_min_conflict_species(const unsigned site, const std::vector<unsigned> free_species, const std::vector<std::vector<int> > &target)
{
    return calc_min_conflict_species(site, free_species, target);
}

std::vector<unsigned> TestMinConf::test_calc_max_conflict_species(const unsigned site,
                                                                  const std::vector<unsigned> pesent_species,
                                                                  const std::vector<std::vector<int> > &target)
{
    return calc_max_conflict_species(site, pesent_species, target);
}

std::vector<std::vector<int> > TestMinConf::test_calculate_commonness()
{
    return calculate_commonness();
}

double TestMinConf::test_calc_energy(const std::vector<std::vector<int> > &commonness, const std::vector<std::vector<int> > &target)
{
    return calc_energy(commonness, target);
}

void TestMinConf::test_gen_init_solution()
{
    void gen_init_solution();
}


TEST_CASE("Tests for the MinConf class") {
    std::vector<unsigned> alpha_list = {2, 1, 2};
    unsigned gamma = 3;
    std::vector<int> target = {  -10, 0, 2 ,
                                 0, -10, 0 ,
                                 2, 0, -10  };

    TestMinConf mc(alpha_list, gamma, target);
    mc.solution = {{0, 1, 1},
                   {1, 0, 0},
                   {0, 1, 1}};

    SECTION("calc_energy() ") {
        std::vector<std::vector<int> > commonness = {{-1, 2, 2},
                                                     {2, -1, 3},
                                                     {2, 3, -1}};
        std::vector<std::vector<int> > target = {{-1, 2, 2},
                                                 {2, -1, 3},
                                                 {2, 3, -1}};

        CHECK(mc.test_calc_energy(commonness, target) == Approx(0.0));
        target[2][1] = 2;
        CHECK(mc.test_calc_energy(commonness, target) == Approx(1.0));
        commonness[0][1] = 5;
        CHECK(mc.test_calc_energy(commonness, target) == Approx(4.0));
    }

    SECTION("calculate_commonness()") {
        std::vector<std::vector<int> > commonness = { { -1, 0, 2 },
                                                      { 0, -1, 0 },
                                                      { 2, 0, -1 } };
        CHECK(mc.test_calculate_commonness() == commonness);
    }

    SECTION("calc_max_conflict_species") {
        mc.solution = {{0, 1, 1},
                       {1, 0, 0},
                       {0, 1, 1}};
        unsigned site = 0;
        std::vector<unsigned> present_species = {1, 2};
        std::vector<std::vector<int> > target = { { -1, 0, 2 },
                                                  { 0, -1, 0 },
                                                  { 2, 0, -1 } };
        CHECK(mc.test_calc_max_conflict_species(site, present_species, target) == present_species);
        mc.solution = {{1, 0, 1},
                       {1, 0, 0},
                       {0, 1, 1}};
        present_species = {0, 2};
        CHECK(mc.test_calc_max_conflict_species(site, present_species, target) == std::vector<unsigned>{0});
        mc.solution = {{1, 1, 1},
                       {1, 0, 0},
                       {0, 1, 1}};
        present_species = {0, 1, 2};
        CHECK(mc.test_calc_max_conflict_species(site, present_species, target) == std::vector<unsigned>{0});
    }

    SECTION("calc_min_conflict_species") {
        mc.solution = {{0, 1, 1},
                       {1, 0, 0},
                       {0, 1, 1}};
        unsigned site = 1;
        std::vector<unsigned> absent_species = {1,2};
        std::vector<std::vector<int> > target = { { -1, 0, 2 },
                                                  { 0, -1, 0 },
                                                  { 2, 0, -1 } };
        CHECK(mc.test_calc_min_conflict_species(site, absent_species, target) == absent_species);
        mc.solution = {{0, 1, 1},
                       {0, 1, 0},
                       {0, 1, 1}};
        absent_species[0] = 0;
        CHECK(mc.test_calc_min_conflict_species(site, absent_species, target) == std::vector<unsigned>{0});
    }

    SECTION("remove_species_max_conf") {
        unsigned site = 0;
        mc.solution = {{1, 1, 1},
                       {1, 0, 0},
                       {0, 1, 1}};
        std::vector<std::vector<int> > expected_solution = {{0, 1, 1},
                                                            {1, 0, 0},
                                                            {0, 1, 1}};
        std::vector<std::vector<int> > target = { { -1, 0, 2 },
                                                  { 0, -1, 0 },
                                                  { 2, 0, -1 } };
        mc.test_remove_species_max_conf(site, target);
        CHECK(mc.solution == expected_solution);
        mc.solution = {{0, 0, 0},
                       {1, 0, 0},
                       {0, 1, 1}};
        expected_solution = {{0, 0, 0},
                             {1, 0, 0},
                             {0, 1, 1}};
        mc.test_remove_species_max_conf(site, target);
        CHECK(mc.solution == expected_solution);
    }

    SECTION("add_species_min_conf") {
        unsigned site = 0;
        mc.solution = {{0, 0, 1},
                       {1, 0, 0},
                       {0, 1, 1}};
        std::vector<std::vector<int> > expected_solution = {{0, 1, 1},
                                                            {1, 0, 0},
                                                            {0, 1, 1}};
        std::vector<std::vector<int> > target = { { -1, 0, 2 },
                                                  { 0, -1, 0 },
                                                  { 2, 0, -1 } };
        mc.test_add_species_min_conf(site, target);
        CHECK(mc.solution == expected_solution);
        mc.test_add_species_min_conf(site, target);
        expected_solution = {{1, 1, 1},
                             {1, 0, 0},
                             {0, 1, 1}};
        CHECK(mc.solution == expected_solution);
        mc.test_add_species_min_conf(site, target);
        CHECK(mc.solution == expected_solution);
    }
}

TEST_CASE("optimizer recognizes a fitting solution") {
    using namespace Rcpp;
    RInside R;
    R.parseEvalQ("total_gamma_sim <- 30; n_sites_sim <- 25; mean_alpha_sim <- 15; sd_sim <- 5");
    R.parseEvalQ("alpha_list_sim <- round(rnorm(n = n_sites_sim, mean = mean_alpha_sim, sd = sd_sim))");
    R.parseEvalQ("alpha_list_sim[alpha_list_sim < 1] <- 1");
    IntegerMatrix true_solution(30, 25);
    const auto alpha_list = as<std::vector<unsigned> >(R.parseEval("alpha_list_sim"));
    const unsigned n_sites = 25;
    const unsigned gamma = 30;
    for (unsigned site = 0; site < n_sites; site++) {
        const unsigned alpha = static_cast<unsigned>(alpha_list[site]);
        IntegerVector site_solution(gamma);
        for (unsigned species = 0; species < alpha; species++) {
            site_solution[species] = 1;
        }
        std::random_shuffle(site_solution.begin(), site_solution.end());
        for (unsigned species = 0; species < gamma; species++) {
            true_solution(species, site) = site_solution[species];
        }
    }
    R["true_solution"] = true_solution;
    IntegerMatrix target = R.parseEval("spectre:::calculate_solution_commonness_rcpp(true_solution)");
}
