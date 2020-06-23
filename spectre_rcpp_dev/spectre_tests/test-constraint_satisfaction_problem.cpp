#include "catch.h"
#include "test-constraint_satisfaction_problem.h"

std::vector<std::vector<int> > Test_Constraint_satisfaction_problem::getTarget() const
{
    return target;
}


TEST_CASE("CSP does not skrew up the target") {
    std::vector<unsigned> alpha_list = {14, 12, 8};
    unsigned gamma = 30;
    std::vector<int> target = {-1, -1, -1,
                               5, -1, -1,
                               3, 4, -1};

    std::vector<std::vector<int> > expected_target = { {-1, 5, 3},
                                                       {5, -1, 4},
                                                       {3, 4, -1} };

    Test_Constraint_satisfaction_problem csp(alpha_list, gamma, target);
    auto calc_target = csp.getTarget();
    CHECK(calc_target == expected_target);
}
