#include "calculate_solution_commonness.h"
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]

IntegerMatrix calculate_solution_commonness_rcpp(const IntegerMatrix solution_matrix) {  
  
  // create results matrix, filled with NAs
  int nsites = solution_matrix.ncol();
  IntegerMatrix commonness_mat(nsites, nsites);
  std::fill(commonness_mat.begin(), commonness_mat.end(), NumericVector::get_na() ) ;
  
  // iterate over the lower triangular matrix without diagonal
  // and calculate # of common species between sites
#pragma omp parallel for schedule(dynamic)
  for (unsigned col = 0; col < nsites - 1; col++) {
    for (unsigned row = col + 1; row < nsites; row++) {
      commonness_mat(row, col) = sum(solution_matrix(_, row) *
        solution_matrix(_, col));
    }
  }
  
  return commonness_mat ;
  
}

// MSP 
NumericMatrix calculate_solution_sorensen_rcpp(const IntegerMatrix solution_matrix) { 
  
  // create result commonness matrix, filled with NAs
  int nsites = solution_matrix.ncol();
  IntegerMatrix commonness_mat(nsites, nsites);
  std::fill(commonness_mat.begin(), commonness_mat.end(), NumericVector::get_na() ) ;
  
  // iterate over the lower triangular matrix without diagonal
  // and calculate # of common species between sites
#pragma omp parallel for schedule(dynamic)
  for (unsigned col = 0; col < nsites - 1; col++) {
    for (unsigned row = col + 1; row < nsites; row++) {
      commonness_mat(row, col) = sum(solution_matrix(_, row) *
        solution_matrix(_, col));
    }
  }
  
  int nspecies = solution_matrix.nrow();
  
  // create results sorensen matrix, filled with NAs
  NumericMatrix sorensen_mat(nsites, nsites);
  std::fill(sorensen_mat.begin(), sorensen_mat.end(), NumericVector::get_na() ) ;
  
  for (unsigned other_site = 0; other_site < nsites; other_site++) {
    for (unsigned site = other_site + 1; site < nsites; site++) {
      if (site == other_site) {
        sorensen_mat(site, other_site) = NA_INTEGER;
        continue;
      } else {
        // sorensen richness variables for site a and site b
        int site_a = 0;
        int site_b = 0;
        for (unsigned species = 0; species < nspecies; species++) { // richness per site 
          if (solution_matrix(species, site)) { // present species at current site
            site_a ++;
          }
          if (solution_matrix(species, other_site)) { // present species at other site
            site_b ++;
          } 
        }
        // sorensen commonness variable 
        int c_temp = commonness_mat(site, other_site);
        sorensen_mat(site, other_site) = 1 - 2.0 * c_temp / (site_a + site_b);
      }
    }
  }
  return sorensen_mat ;
}
// MSP end

std::vector<int> calculate_solution_commonness(const std::vector<int> &solution_matrix, 
                                               const unsigned n_sites,
                                               const unsigned n_species) {
  
  // create results vector
  std::vector<int> result(n_sites * n_sites, 0);
  
  // iterate over the lower triangular matrix without diagonal
  // and calculate # of common species between sites
#pragma omp parallel for
  for (unsigned site = 0; site < n_sites; site++) {
    update_solution_commonness_site(solution_matrix, result, n_sites, n_species, site);
  }
  
  return result ;
}

// MSP start
std::vector<double> calculate_solution_sorensen(const std::vector<int> &solution_matrix,
                                                const unsigned n_sites,
                                                const unsigned n_species) {
  
  // create commonness vector
  std::vector<int> commonness_vector(n_sites * n_sites, 0);
  
  // iterate over the lower triangular matrix without diagonal
  // and calculate # of common species between sites
#pragma omp parallel for
  for (unsigned site = 0; site < n_sites; site++) {
    update_solution_commonness_site(solution_matrix, commonness_vector, n_sites, n_species, site); 
  }
  // calculate sorensen matrix from commonness_vector vector 
  // create sorensen vector
  std::vector<double> sorensen_vector(n_sites * n_sites, 0);
  
  for (unsigned site = 0; site < n_sites; site++) {
    update_solution_sorensen_site(solution_matrix, commonness_vector, sorensen_vector, n_sites, n_species, site); 
  }
  
  return sorensen_vector;
}
// MSP end 



std::vector<int> calculate_solution_commonness_site(const std::vector<int> &solution_matrix,
                                                    const std::vector<int> &solution_commonness,
                                                    const unsigned n_sites,
                                                    const unsigned n_species,
                                                    const unsigned site) {
  // create results vector
  std::vector<int> result = solution_commonness;
  update_solution_commonness_site(solution_matrix, result, n_sites, n_species, site);
  
  return result ;
}

// MSP
std::vector<double> calculate_solution_sorensen_site(const std::vector<int> &solution_matrix, // 
                                                     const std::vector<int> &solution_commonness,
                                                     const std::vector<double> &solution_sorensen,
                                                     const unsigned n_sites,
                                                     const unsigned n_species,
                                                     const unsigned site) {
  std::vector<int> result_commonness = solution_commonness;
  update_solution_commonness_site(solution_matrix, result_commonness, n_sites, n_species, site);
  
  // create results vector
  std::vector<double> result_sorensen = solution_sorensen;
  update_solution_sorensen_site(solution_matrix, result_commonness, result_sorensen, n_sites, n_species, site);
  
  return result_sorensen ;
}
// MSP end 


IntegerMatrix calculate_solution_commonness_site_rcpp(const IntegerMatrix solution_matrix, // after swich-try 
                                                      const IntegerMatrix solution_commonness, // commonness matrix before switch 
                                                      const int site) {
  
  const int nrows = solution_commonness.nrow();
  const int site_ = site - 1; // because C++ starts indexing at zero
  IntegerMatrix new_solution_commonness = solution_commonness;
#pragma omp parallel
{
#pragma omp for nowait schedule(static)
  for (int j = 0; j < site_; j++) {
    new_solution_commonness(site_, j) = sum(solution_matrix(_, j) *
      solution_matrix(_, site_));
  }
  
#pragma omp for schedule(static)
  for (int j = site_ + 1; j < nrows; j++) {
    new_solution_commonness(j, site_) = sum(solution_matrix(_, j) *
      solution_matrix(_, site_));
  }
}
return new_solution_commonness;
}

// MSP: BUG! ,but not used at the moment 
NumericMatrix calculate_solution_sorensen_site_rcpp(const IntegerMatrix solution_matrix, 
                                                    const NumericMatrix solution_sorensen,
                                                    const int site) {
  
  const int nrows = solution_sorensen.nrow();
  const int site_ = site - 1; // because C++ starts indexing at zero
  NumericMatrix new_solution_sorensen = solution_sorensen;
#pragma omp parallel
{
#pragma omp for nowait schedule(static)
  for (int j = 0; j < site_; j++) {
    new_solution_sorensen(site_, j) = sum(solution_matrix(_, j) *
      solution_matrix(_, site_));
  }
  
#pragma omp for schedule(static)
  for (int j = site_ + 1; j < nrows; j++) {
    new_solution_sorensen(j, site_) = sum(solution_matrix(_, j) *
      solution_matrix(_, site_));
  }
}
return new_solution_sorensen;
}

// MSP end


IntegerMatrix calculate_solution_commonness_species_site_rcpp(const IntegerMatrix solution_matrix,
                                                              const IntegerMatrix solution_commonness,
                                                              const int site, const int species)
{
  const int nrows = solution_commonness.nrow();
  const int site_ = site - 1; // because C++ starts indexing at zero
  IntegerMatrix new_solution_commonness = solution_commonness;
#pragma omp parallel
{
#pragma omp for nowait schedule(static)
  for (int j = 0; j < site_; j++) {
    new_solution_commonness(site_, j) = sum(solution_matrix(_, j) *
      solution_matrix(_, site_));
  }
  
#pragma omp for schedule(static)
  for (int j = site_ + 1; j < nrows; j++) {
    new_solution_commonness(j, site_) = sum(solution_matrix(_, j) *
      solution_matrix(_, site_));
  }
}
return new_solution_commonness;
}


void update_solution_commonness_site_rcpp(const IntegerMatrix solution_matrix,
                                          IntegerMatrix &solution_commonness,
                                          const unsigned site) {
  
  const int nrows = solution_commonness.nrow();
  const int site_ = site - 1; // because C++ starts indexing at zero
  
  for (int j = 0; j < site_; j++) {
    solution_commonness(site_, j) = sum((solution_matrix(_, j) +
      solution_matrix(_, site_)) > 1.0);
  }
  
  for (int j = site_ + 1; j < nrows; j++) {
    solution_commonness(j, site_) = sum((solution_matrix(_, j) +
      solution_matrix(_, site_)) > 1.0);
  }
  
}

void update_solution_commonness_site(const std::vector<int> &solution_matrix, // MSP comment: currently under use
                                     std::vector<int> &solution_commonness,
                                     const unsigned n_sites,
                                     const unsigned n_species,
                                     const unsigned site) {
  for (unsigned other_site = 0; other_site < n_sites; other_site++) {
    if (site == other_site) {
      solution_commonness[site * n_sites + other_site] = NA_INTEGER;
      continue;
    } else {
      for (unsigned species = 0; species < n_species; species++) {
        if (!solution_matrix[site * n_species + species]) { // no species at current site
          continue;
        } else if (!solution_matrix[other_site * n_species + species]) { // no species at other site
          continue;
        } else {
          solution_commonness[site * n_sites + other_site]++;
        }
      }
    }
  }
}

// MSP start 
void update_solution_sorensen_site(const std::vector<int> &solution_matrix, 
                                   std::vector<int> &commonness_vector,
                                   std::vector<double> &sorensen_vector,
                                   const unsigned n_sites,
                                   const unsigned n_species,
                                   const unsigned site) {
  for (unsigned other_site = 0; other_site < n_sites; other_site++) {
    if (site == other_site) {
      sorensen_vector[site * n_sites + other_site] = NA_REAL; // MSP latest change... 
      continue;
    } else {
      // sorensen richness variables for site a and site b
      int site_a = 0;
      int site_b = 0;
      for (unsigned species = 0; species < n_species; species++) { // richness per site 
        if (solution_matrix[site * n_species + species]) { // present species at current site
          site_a ++;
        }
        if (solution_matrix[other_site * n_species + species]) { // present species at other site
          site_b ++;
        } 
      }
      // sorensen commonness variable 
      int c_temp = commonness_vector[site * n_sites + other_site];
      sorensen_vector[site * n_sites + other_site] = 1 - 2.0 * c_temp / (site_a + site_b);
    }
  }
}

// MSP end

/*** R
# Big
new_solution <- matrix(sample(0:1, 480000, replace = TRUE), 6000, 8000)
solution_commonness <- matrix(sample(0:8, 640000, replace = TRUE), 8000, 8000)

# "Small"
new_solution <- matrix(sample(0:1, 200000, replace = TRUE), 2000, 3000)
solution_commonness <- matrix(sample(0:8, 90000, replace = TRUE), 3000, 3000)

random_col <- 68

calculate_solution_commonness_site_R <- function(solution_matrix, solution_commonness, site){
  
  for(j in 1:nrow(solution_commonness)){
    if (j < site){
      solution_commonness[site,j] <- sum((solution_matrix[,j] + solution_matrix[,site]) > 1)
    }
    if (j > site){
      solution_commonness[j,site] <- sum((solution_matrix[,j] + solution_matrix[,site]) > 1)
    }
    
  }
  return(solution_commonness)
}

calculate_solution_commonness_site_R2 <- function(solution_matrix, solution_commonness, site){
  
  for (j in 1:site) {
    solution_commonness[site, j] = sum((solution_matrix[, j] + solution_matrix[, site]) > 1)
  }
  
  for (j in site:nrow(solution_commonness)) {
    solution_commonness[j, site] = sum((solution_matrix[, j] + solution_matrix[, site]) > 1)
  }
  
  return(solution_commonness)
}

# Identical:
identical(calculate_solution_commonness_site_R(new_solution, solution_commonness, random_col), spectre:::calculate_solution_commonness_site_rcpp(new_solution, solution_commonness, random_col))
# Not identical:
identical(calculate_solution_commonness_site_R2(new_solution, solution_commonness, random_col), calculate_solution_commonness_site_R(new_solution, solution_commonness, random_col))

bm <- bench::mark(spectre:::calculate_solution_commonness_rcpp_p(new_solution),
                  spectre:::calculate_solution_commonness_rcpp(new_solution),
                  spectre:::calculate_solution_commonness_rcpp_old(new_solution),
                  check = FALSE) # check FALSE is needed because new variant calculates the lower triangular matrix, only

# 1 spectre:::calculate_solution_commonness_rcpp_p(new_solution)   3.07s  3.07s   0.326      34.3MB        0     1     0
# 2 spectre:::calculate_solution_commonness_rcpp(new_solution)     1.12m  1.12m   0.0148     34.3MB        0     1     0
# 3 spectre:::calculate_solution_commonness_rcpp_old(new_solution) 2.54m  2.54m   0.00657    34.3MB        0     1     0

bm2 <- bench::mark(
  spectre:::calculate_solution_commonness_site_rcpp_p(new_solution, solution_commonness, random_col),
  spectre:::calculate_solution_commonness_site_rcpp(new_solution, solution_commonness, random_col), iterations = 100)

# 1 spectre:::calculate_solution_commonness_site_rcpp_p(new_solution, solution_commonness, random_col)  2.05ms  3.2ms
# 2 spectre:::calculate_solution_commonness_site_rcpp(new_solution, solution_commonness, random_col)   46.82ms 46.9ms

*/
