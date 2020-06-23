#ifndef BACKTRACKING_H
#define BACKTRACKING_H
#include "constraint_satisfaction_problem.h"

class Backtracking : public Constraint_satisfaction_problem
{
public:
    using Constraint_satisfaction_problem::Constraint_satisfaction_problem;
    int optimize(long max_steps_ = 5000);
    std::vector<unsigned> solved_sites;

private:
    bool iterate_species(std::vector<int> &site);
    bool check_consistency(const int site);
    int backtrack(const int site);
    void reset_site(int site);
    template <class T>
    T count_before(const T x, const std::vector<int> &site, const unsigned species);
};

#endif // BACKTRACKING_H
