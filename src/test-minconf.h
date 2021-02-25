#ifndef TESTMINCONF_H
#define TESTMINCONF_H
#include <testthat.h>
#include "minconf.h"


class TestMinConf : public MinConf
{
public:
    using MinConf::MinConf;
    void test_add_species_min_conf(unsigned site, const std::vector<std::vector<float> > &target);
    std::vector<unsigned> test_calc_min_conflict_species(const unsigned site,
                                                    const std::vector<std::vector<float> > &target);
    void test_update_solution_commonness();
    float test_calc_error(const std::vector<std::vector<float> > &target);
    std::vector<std::vector<float> > getTarget() const;
    std::vector<unsigned> test_calc_missing_species();
    static const int NA = -2147483648; // Rcpp NA_INTEGER == R_NaInt
};

#endif // TESTMINCONF_H
