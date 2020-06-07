#include "catch.h"
#include "test-minconf.h"

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
        CHECK(mc.test_calc_energy(commonness, target) == Approx(1.0/10));
        commonness[0][1] = 5;
        CHECK(mc.test_calc_energy(commonness, target) == Approx(4.0/10));
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

TEST_CASE("optimizer recognizes a solution") {
    std::vector<unsigned> alpha_list = {14, 12, 8, 12, 15, 10, 10, 15, 5, 13, 11, 11, 11, 8, 8, 13, 13, 11, 6, 12, 9, 11, 8, 12, 10};
    unsigned gamma = 30;
    std::vector<int> target = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                               5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                               3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                               7, 7, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                               6, 9, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                               6, 5, 4, 5, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                               7, 4, 4, 5, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                               8, 4, 6, 7, 6, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                               2, 2, 2, 2, 2, 1, 2, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                               4, 2, 4, 3, 5, 3, 3, 6, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                               5, 4, 3, 6, 6, 3, 3, 6, 2, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                               5, 3, 2, 6, 4, 5, 3, 7, 0, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                               4, 3, 3, 2, 6, 3, 5, 7, 3, 6, 2, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                               4, 4, 1, 3, 3, 5, 2, 3, 0, 4, 0, 3, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                               3, 2, 4, 3, 3, 2, 3, 3, 1, 5, 3, 2, 3, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                               6, 3, 2, 6, 7, 2, 4, 8, 2, 7, 5, 4, 5, 4, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                               5, 8, 5, 6, 6, 5, 5, 5, 3, 5, 3, 5, 5, 4, 5, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                               6, 5, 3, 6, 6, 3, 3, 3, 2, 5, 7, 4, 2, 2, 4, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1,
                               3, 3, 1, 3, 1, 3, 3, 4, 1, 3, 1, 3, 2, 3, 1, 1, 3, 1, -1, -1, -1, -1, -1, -1, -1,
                               5, 6, 4, 6, 7, 5, 6, 6, 1, 5, 4, 5, 4, 3, 5, 5, 7, 4, 2, -1, -1, -1, -1, -1, -1,
                               7, 2, 1, 4, 5, 4, 3, 6, 3, 3, 3, 4, 5, 3, 1, 4, 2, 3, 2, 1, -1, -1, -1, -1, -1,
                               7, 5, 1, 5, 7, 4, 2, 3, 1, 4, 4, 5, 3, 3, 3, 3, 5, 5, 2, 4, 5, -1, -1, -1, -1,
                               3, 3, 1, 3, 4, 2, 1, 1, 2, 4, 2, 3, 3, 2, 2, 3, 4, 3, 0, 5, 2, 4, -1, -1, -1,
                               7, 3, 2, 4, 7, 3, 5, 6, 2, 7, 5, 4, 4, 3, 3, 7, 4, 6, 1, 5, 5, 5, 4, -1, -1,
                               3, 4, 2, 3, 5, 2, 3, 5, 2, 6, 2, 4, 4, 3, 2, 5, 4, 1, 3, 5, 2, 4, 3, 5, -1};

    TestMinConf mc(alpha_list, gamma, target);
    //    mc.solution = { {1, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0},
    //                    {0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0},
    //                    {0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0},
    //                    {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0},
    //                    {0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 1},
    //                    {0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0},
    //                    {1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0},
    //                    {0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    //                    {0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1},
    //                    {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1},
    //                    {1, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0},
    //                    {1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0},
    //                    {1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0},
    //                    {1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0},
    //                    {0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1},
    //                    {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0},
    //                    {1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0},
    //                    {0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0},
    //                    {1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 1},
    //                    {0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0},
    //                    {0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1},
    //                    {1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1},
    //                    {0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1},
    //                    {1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1},
    //                    {1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0},
    //                    {0, 0, 0, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0},
    //                    {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1},
    //                    {1, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0},
    //                    {0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0},
    //                    {1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0} };


    //    SECTION("Min-conf-1") {
    //        int iter = mc.optimize1(5, 0.0, 0, false);
    //        auto energy = mc.energy_vector;
    //        CHECK(iter == 5);
    //    }

}
