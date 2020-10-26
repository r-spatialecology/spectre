#ifndef TEST_CONSTRAINT_SATISFACTION_PROBLEM_H
#define TEST_CONSTRAINT_SATISFACTION_PROBLEM_H
#include "../../src/constraint_satisfaction_problem.h"


class Test_Constraint_satisfaction_problem : public Constraint_satisfaction_problem
{
public:
    using Constraint_satisfaction_problem::Constraint_satisfaction_problem;

    std::vector<std::vector<int> > getTarget() const;
};

#endif // TEST_CONSTRAINT_SATISFACTION_PROBLEM_H
