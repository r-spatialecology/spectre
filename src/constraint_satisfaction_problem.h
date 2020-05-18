#ifndef CONSTRAINT_SATISFACTION_PROBLEM_H
#define CONSTRAINT_SATISFACTION_PROBLEM_H
#include <vector>

class Constraint_satisfaction_problem
{
public:
    Constraint_satisfaction_problem(std::vector<unsigned> alpha_list,
                                    unsigned gamma_div, std::vector<int> target);

    int optimize(long max_steps_ = 5000);
    std::vector<std::vector<int> > solution;
    std::vector<unsigned> solved_sites;

private:
    long max_steps = 0;

    std::vector<std::vector<int> > target;
    const std::vector<unsigned> alpha_list;
    const unsigned gamma_div;
    const unsigned n_sites;
    bool iterate_species(std::vector<int> &site);
    bool check_consistency(const int site);
    int backtrack(const int site);
    void reset_site(int site);
    template <class T>
    T count_before(const T x, const std::vector<int> &site, const unsigned species);
};

#endif // CONSTRAINT_SATISFACTION_PROBLEM_H
