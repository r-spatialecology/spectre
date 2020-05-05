#include "constraint_satisfaction_problem.h"
#include "Rcpp.h"

using namespace Rcpp;

Constraint_satisfaction_problem::Constraint_satisfaction_problem(std::vector<unsigned> alpha_list,
                                                                 unsigned gamma_div,
                                                                 std::vector<int> target_)
    : alpha_list(alpha_list), gamma_div(gamma_div), n_sites(alpha_list.size())
{
    solution.resize(n_sites);
    target.resize(n_sites);

    for (unsigned site = 0; site < n_sites; site++) {
        solution[site].resize(gamma_div);
        for (unsigned species = 0; species < alpha_list[site]; species++) {
            solution[site][species] = 1;
        }

        // convert target matrix to a more convenient format
        target[site].resize(n_sites);
        for (unsigned other_site = 0; other_site < n_sites; other_site++) {
            if (site == other_site) {
                target[site][other_site] = -1; // i.e. NA
            } else if (target_[other_site * n_sites + site] < 0) { // i.e. NA
                target[site][other_site] = target_[site * n_sites + other_site];
            } else {
                target[site][other_site] = target_[other_site * n_sites + site];
            }
        }
    }
}

int Constraint_satisfaction_problem::optimize(long max_steps_)
{
    int site = 1;
    max_steps = max_steps_;

    if (max_steps_ > 0) { // if max steps are given
        while (max_steps-- > 0 && site < n_sites && site >= 0) {
            if (check_consistency(site)) {
                site++;
                Rcout << "\n optimized to site: " << site;
            } else if (!iterate_species(solution[site])) {
                site = backtrack(site);
                Rcout << "\n backtracked to site: " << site;
            }
        }
    } else { // no max_steps are given, just try your best forever
        while (site < n_sites && site >= 0) {
           // max_steps--;
            if (check_consistency(site)) {
                site++;
                Rcout << "\n optimized to site: " << site;
            } else if (!iterate_species(solution[site])) {
                site = backtrack(site);
                Rcout << "\n backtracked to site: " << site;
            }
        }
    }

    return max_steps;
}

bool Constraint_satisfaction_problem::iterate_species(std::vector<int> &site)
{
    for (int species = gamma_div - 1; species >= 0; species--) {
        if (site[species] == 1) {
            if (species < gamma_div - 1) { // there is a zero after the one
                site[species] = 0;
                site[species + 1] = 1;
                return true;
            } else {
                const unsigned species_before = count_before(1, site, species);
                const unsigned first_in_row = species - species_before;
                const unsigned zeros_before = count_before(0, site, first_in_row);
                if (zeros_before) { // last one, but zeros before and not only zeros before
                    if (first_in_row - zeros_before == 0)
                        return false; // only zeros before; end of iterations

                    for (unsigned s = 0; s <= species_before; s++) {
                        const unsigned spec = first_in_row - zeros_before + s;
                        site[first_in_row + s] = 0;
                        site[spec] = 1;
                    }

                    site[first_in_row - zeros_before - 1] = 0;
                    site[first_in_row - zeros_before + species_before + 1] = 1;

                    return true;
                }
            }
            return false; // only ones in this site (all species present)
        }
    }
    return false; // only zeros in this site (no species present)
}

bool Constraint_satisfaction_problem::check_consistency(const int site)
{

    for (unsigned other_site = 0; other_site < site; other_site++) {
        int common_species = 0;
        for (unsigned species = 0; species < gamma_div; species++) {
            common_species += solution[site][species] * solution[other_site][species];
        }
        if (common_species != target[site][other_site]) {
            return false;
        }
    }

    return true;
}

int Constraint_satisfaction_problem::backtrack(const int site)
{
    if (site == 0)
        return -1; // i.e. NA

    int backtracked_site = site;
    do {
        // reset current site
        reset_site(backtracked_site);
        backtracked_site--;
        if (!iterate_species(solution[backtracked_site])) { // is actually an additional step
            max_steps--; // could become < 0 here...
            backtrack(backtracked_site); // become recursive
        }
        max_steps--; // could become < 0 here...

    } while (!check_consistency(backtracked_site));

    return backtracked_site;
}

void Constraint_satisfaction_problem::reset_site(int site)
{
    for (unsigned species = 0; species < gamma_div; species++) {
        if (species < alpha_list[site]) {
            solution[site][species] = 1;
        } else {
            solution[site][species] = 0;
        }
    }
}

template <class T>
T Constraint_satisfaction_problem::count_before(const T x,
                                                const std::vector<int> &site,
                                                const unsigned species)
{
    unsigned counter = 0;
    for (int pos = species - 1; pos >= 0; pos--) {
        if (site[pos] != x) {
            return counter;
        }
        counter++;
    }
    return counter;
}
