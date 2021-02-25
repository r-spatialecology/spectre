#include "test-minconf.h"

void TestMinConf::test_add_species_min_conf(
    unsigned site, const std::vector<std::vector<float>> &target) {
  this->target = target;
  add_species_min_conf(site);
}

std::vector<unsigned> TestMinConf::test_calc_min_conflict_species(
    const unsigned site, const std::vector<std::vector<float>> &target) {
  this->target = target;
  return calc_min_conflict_species(site);
}

void TestMinConf::test_update_solution_commonness() {
  update_solution_commonness();
}

float TestMinConf::test_calc_error(
    const std::vector<std::vector<float>> &target) {
  this->target = target;
  return calc_error();
}

std::vector<unsigned> TestMinConf::test_calc_missing_species() {
    return calc_missing_species();
}

std::vector<std::vector<float>> TestMinConf::getTarget() const {
  return target;
}

context("MinConf does not skrew up the target") {
  std::random_device rd;
  long long seed = rd();
  std::vector<unsigned> alpha_list = {14, 12, 8};
  unsigned gamma = 30;

  std::vector<int> target = {TestMinConf::NA,
                             TestMinConf::NA,
                             TestMinConf::NA,
                             5,
                             TestMinConf::NA,
                             TestMinConf::NA,
                             3,
                             4,
                             TestMinConf::NA};

  std::vector<std::vector<float>> expected_target = {
      {TestMinConf::NA, 0.615385, 0.727273},
      {TestMinConf::NA, TestMinConf::NA, 0.6},
      {TestMinConf::NA, TestMinConf::NA, TestMinConf::NA}};

  TestMinConf min_conf(alpha_list, gamma, target, std::vector<int>(),
                       std::vector<int>(), seed);
  auto calc_target = min_conf.getTarget();
  for (unsigned site = 0; site < expected_target.size(); site++) {
    for (unsigned other_site = 0; other_site < expected_target.size();
         other_site++) {
      expect_true(calc_target[site][other_site] ==
                    Approx(expected_target[site][other_site]));
    }
  }
}

context("Tests for the MinConf class") {
  std::random_device rd;
  long long seed = rd();
  std::vector<unsigned> alpha_list = {2, 1, 2};
  unsigned gamma = 3;
  std::vector<int> target = {-10, 0, 2, 0, -10, 0, 2, 0, -10};

  test_that("calc_missing_species()") {
    TestMinConf mc(alpha_list, gamma, target, std::vector<int>(),
                   std::vector<int>(), seed);
    mc.solution = {{1, 1, 1}, {0, 0, 0}, {1, 0, 0}};
    const auto missing_species = mc.test_calc_missing_species();

    expect_true(missing_species.size() == 3);
    expect_true(missing_species[0] == 0);
    expect_true(missing_species[1] == 1);
    expect_true(missing_species[2] == 1);
  }

  TestMinConf mc(alpha_list, gamma, target, std::vector<int>(),
                 std::vector<int>(), seed);
  mc.solution = {{0, 1, 1}, {1, 0, 0}, {0, 1, 1}};

  test_that("calc_error() ") {
    std::vector<std::vector<float>> target = {{TestMinConf::NA, 1, 0},
                                              {TestMinConf::NA, TestMinConf::NA, 1},
                                              {TestMinConf::NA, TestMinConf::NA, TestMinConf::NA}};

    expect_true(mc.test_calc_error(target) == Approx(0.0));
    target[0][1] = 2;
    expect_true(mc.test_calc_error(target) == Approx(1.0));
  }

  test_that("calculate_commonness()") {
    std::vector<std::vector<int>> commonness = {
        {TestMinConf::NA, 0, 2},
        {TestMinConf::NA, TestMinConf::NA, 0},
        {TestMinConf::NA, TestMinConf::NA, TestMinConf::NA}};
    mc.test_update_solution_commonness();
    expect_true(mc.commonness == commonness);
  }

  test_that("calc_min_conflict_species") {
    mc.solution = {{0, 1, 1}, {0, 0, 0}, {0, 1, 1}};
    unsigned site = 1;
    std::vector<unsigned> expected_species = {1, 2};
    std::vector<std::vector<float>> target = {{TestMinConf::NA, 0.667, 0},
                                              {TestMinConf::NA, TestMinConf::NA, 0},
                                              {TestMinConf::NA, TestMinConf::NA, TestMinConf::NA}};
    expect_true(mc.test_calc_min_conflict_species(site, target) ==
                expected_species);
    mc.solution = {{0, 1, 1}, {0, 1, 0}, {0, 1, 0}};
    site = 2;
    expected_species = {2};
    expect_true(mc.test_calc_min_conflict_species(site, target) ==
                expected_species);
  }

  test_that("add_species_min_conf") {
    unsigned site = 0;
    mc.solution = {{0, 0, 1}, {0, 0, 0}, {0, 1, 1}};
    std::vector<std::vector<int>> expected_solution = {
        {0, 1, 1}, {0, 0, 0}, {0, 1, 1}};
    std::vector<std::vector<float>> target = {{TestMinConf::NA, 0, 0},
                                              {TestMinConf::NA, TestMinConf::NA, 0},
                                              {TestMinConf::NA, TestMinConf::NA, TestMinConf::NA}};
    mc.test_add_species_min_conf(site, target);
    expect_true(mc.solution == expected_solution);
    mc.test_add_species_min_conf(site, target);
    expected_solution = {{1, 1, 1}, {0, 0, 0}, {0, 1, 1}};
    expect_true(mc.solution == expected_solution);
    mc.test_add_species_min_conf(site, target);
    expect_true(mc.solution == expected_solution);
  }
}
