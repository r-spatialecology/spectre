#include "test-minconf.h"

void TestMinConf::test_add_species_min_conf(unsigned site, const std::vector<std::vector<int> > &target)
{
  add_species_min_conf(site, target);
}

std::vector<unsigned> TestMinConf::test_calc_min_conflict_species(const unsigned site, const std::vector<unsigned> free_species, const std::vector<std::vector<int> > &target)
{
  return calc_min_conflict_species(site, free_species, target);
}


void TestMinConf::test_update_solution_commonness()
{
  update_solution_commonness();
}

double TestMinConf::test_calc_error(const std::vector<std::vector<int> > &commonness, const std::vector<std::vector<int> > &target)
{
  return calc_error(commonness, target);
}

void TestMinConf::test_gen_init_solution()
{
  void gen_init_solution();
}

std::vector<unsigned> TestMinConf::test_calc_missing_species()
{
  return calc_missing_species();
}

std::vector<std::vector<int> > TestMinConf::getTarget() const
{
  return target;
}


context("MinConf does not skrew up the target") {
  std::random_device rd;
  long long seed = rd();
  std::vector<unsigned> alpha_list = {14, 12, 8};
  unsigned gamma = 30;

  std::vector<int> target = {TestMinConf::NA, TestMinConf::NA, TestMinConf::NA,
                             5, TestMinConf::NA, TestMinConf::NA,
                             3, 4, TestMinConf::NA};

  std::vector<std::vector<int> > expected_target = { {TestMinConf::NA, 5, 3},
                                                     {TestMinConf::NA, TestMinConf::NA, 4},
                                                     {TestMinConf::NA, TestMinConf::NA, TestMinConf::NA} };

  TestMinConf min_conf(alpha_list, gamma, target, std::vector<int>(), std::vector<int>(), seed);
  auto calc_target = min_conf.getTarget();
  expect_true(calc_target == expected_target);
}


context("Tests for the MinConf class") {
  std::random_device rd;
  long long seed = rd();
  std::vector<unsigned> alpha_list = {2, 1, 2};
  unsigned gamma = 3;
  std::vector<int> target = {  -10, 0, 2 ,
                               0, -10, 0 ,
                               2, 0, -10  };

  test_that("calc_missing_species()") {
    TestMinConf mc(alpha_list, gamma, target, std::vector<int>(), std::vector<int>(), seed);
    mc.solution = {{1, 1, 1},
                   {0, 0, 0},
                   {1, 0, 0}};
    const auto missing_species = mc.test_calc_missing_species();

    expect_true(missing_species.size() == 3);
    expect_true(missing_species[0] == 0);
    expect_true(missing_species[1] == 1);
    expect_true(missing_species[2] == 1);
  }

  TestMinConf mc(alpha_list, gamma, target, std::vector<int>(), std::vector<int>(), seed);
  mc.solution = {{0, 1, 1},
                 {1, 0, 0},
                 {0, 1, 1}};

  test_that("calc_error() ") {
    std::vector<std::vector<int> > commonness = {{-1, 2, 2},
                                                 {2, -1, 3},
                                                 {2, 3, -1}};
    std::vector<std::vector<int> > target = {{-1, 2, 2},
                                             {2, -1, 3},
                                             {2, 3, -1}};

    expect_true(mc.test_calc_error(commonness, target) == Approx(0.0));
    target[2][1] = 2;
    expect_true(mc.test_calc_error(commonness, target) == Approx(1.0));
    commonness[0][1] = 5;
    expect_true(mc.test_calc_error(commonness, target) == Approx(4.0));
  }

  test_that("calculate_commonness()") {
    std::vector<std::vector<int> > commonness = { { -1, 0, 2 },
                                                  { -1, -1, 0 },
                                                  { -1, -1, -1 } };
    mc.test_update_solution_commonness();
    expect_true(mc.commonness == commonness);
  }

  test_that("calc_min_conflict_species") {
    mc.solution = {{0, 1, 1},
                   {1, 0, 0},
                   {0, 1, 1}};
    unsigned site = 1;
    std::vector<unsigned> absent_species = {1,2};
    std::vector<std::vector<int> > target = { { -1, 0, 2 },
                                              { 0, -1, 0 },
                                              { 2, 0, -1 } };
    expect_true(mc.test_calc_min_conflict_species(site, absent_species, target) == absent_species);
    mc.solution = {{0, 1, 1},
                   {0, 1, 0},
                   {0, 1, 1}};
    absent_species[0] = 0;
    expect_true(mc.test_calc_min_conflict_species(site, absent_species, target) == std::vector<unsigned>{0});
  }

  test_that("add_species_min_conf") {
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
    expect_true(mc.solution == expected_solution);
    mc.test_add_species_min_conf(site, target);
    expected_solution = {{1, 1, 1},
                         {1, 0, 0},
                         {0, 1, 1}};
    expect_true(mc.solution == expected_solution);
    mc.test_add_species_min_conf(site, target);
    expect_true(mc.solution == expected_solution);
  }
}
