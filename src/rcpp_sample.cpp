#include "rcpp_sample.h"

IntegerVector rcpp_sample(IntegerVector x, int size, bool replace) {

  return sample(x, size, replace);

}

IntegerVector which_not(const IntegerVector x, const int y) {
  IntegerVector v = seq(1, x.size());
  return v[x != y];
}

IntegerVector which_not_vec(const IntegerVector x, const IntegerVector y) {
  IntegerVector v = seq(1, x.size());
  return v[x != y];
}

IntegerVector get_swap_rows_rcpp(const IntegerVector &v, bool cpp_index) {
  RNGScope scope; // needed for debugging
  unsigned species1 = sample(IntegerVector(seq(1, v.size())), 1)[0];
  unsigned species2 = sample(which_not(v, v[species1 - 1]), 1)[0];  // R starts counting at 1
  if(cpp_index) { // convert to zero-based indexing
      species1--;
      species2--;
  }
  return(IntegerVector::create(species1, species2));
}

IntegerVector get_swap_rows_rcpp_bruteforce(const IntegerVector &v, bool cpp_index) {
  RNGScope scope; // needed for debugging
  IntegerVector rows = seq(1, v.size());
  unsigned species1 = sample(rows, 1)[0];
  unsigned species2 = sample(rows, 1)[0];

  ///TODO: This will run infinite if the vector contains only zeros or only ones (i.e. cannot swapped)
  while (v[species1 - 1] == v[species2 - 1]) {  // R starts counting at 1
    species2 = sample(rows, 1)[0];
  }
  if(cpp_index) { // convert to zero-based indexing
      species1--;
      species2--;
  }
  return(IntegerVector::create(species1, species2));
}

/*** R
rcpp_sample(x = 1:100, size = 5)
rcpp_sample(x = 1:5, size = 10, replace = TRUE)

bench::mark(rcpp_sample(x = 1:100, size = 5),
            sample(x = 1:100, size = 5),
            check = FALSE, relative = TRUE, iterations = 10000)
# 1 rcpp_sample(x = 1:100, size = 5)  1      2.01      1         1.15      NaN 10000     0      183ms <int … <df[,…
# 2 sample(x = 1:100, size = 5)       1.56   1         2.95      1         Inf  9999     1       62ms <int … <df[,…
                                                                                                                                                                                                                              # … with 2 more variables: time <list>, gc <list>
bench::mark(rcpp_sample(x = 1:10, size = 20, replace = TRUE),
            sample(x = 1:10, size = 20, replace = TRUE),
            check = FALSE, relative = TRUE, iterations = 10000)
# 1 rcpp_sample(x = 1:10, size = 20, replace = TRUE)  1      1         1            1      NaN 10000     0    152.2ms
# 2 sample(x = 1:10, size = 20, replace = TRUE)       1.69   1.02      1.54         1      Inf  9999     1     98.7ms

mat <- matrix(sample(0:1, 480000, replace = TRUE), 6000, 8000)
site <- 68

x <- spectre:::get_swap_rows(mat[, site])
mat[x[1], site] != mat[x[2], site]
x <- spectre:::get_swap_rows_rcpp(mat[, site])
mat[x[1], site] != mat[x[2], site]
x <- spectre:::get_swap_rows_rcpp_bruteforce(mat[, site])
mat[x[1], site] != mat[x[2], site]

bench::mark(spectre:::get_swap_rows(mat[, site]),
            spectre:::get_swap_rows_rcpp(mat[, site]),
            spectre:::get_swap_rows_rcpp_bruteforce(mat[, site]),
            check = FALSE, relative = TRUE, iterations = 10000)

# Large, "dense" matrix (i.e. 'mat' from above)
# 1 spectre:::get_swap_rows(mat[, site])                  2.96   3.45      1         1.84     1.13  9985    15
# 2 spectre:::get_swap_rows_rcpp(mat[, site])             3.45   3.95      1.11      1.81     1     9988    12
# 3 spectre:::get_swap_rows_rcpp_bruteforce(mat[, site])  1      1         4.03      1        1.81  9994     6

# small, "sparse" matrix (i.e. 'current_solution' from vignette)
# 1 spectre:::get_swap_rows(mat[, site])                  1.68   2.16      1         1.67     1.36  9997     3
# 2 spectre:::get_swap_rows_rcpp(mat[, site])             1.17   1.18      2.20      1.22     1     9999     1
# 3 spectre:::get_swap_rows_rcpp_bruteforce(mat[, site])  1      1         2.43      1        2.21  9998     2
*/
