% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_optimization_min_conf.R
\name{run_optimization_min_conf}
\alias{run_optimization_min_conf}
\title{run_optimization_min_conf}
\usage{
run_optimization_min_conf(
  alpha_list,
  total_gamma,
  target,
  max_iterations,
  partial_solution = NULL,
  fixed_species = NULL,
  autostop = 0,
  seed = NA,
  verbose = TRUE,
  interruptible = TRUE
)
}
\arguments{
\item{alpha_list}{\code{Matrix} of predicted alpha diversity (species richness) in
each cell.}

\item{total_gamma}{Total number of species present throughout the entire
landscape.}

\item{target}{Pairwise matrix of species in common between each site by site
pair. Only the upper triangle of the matrix is actually needed.}

\item{max_iterations}{The maximum number of iterations that the optimization
algorithm may run through before stopping.}

\item{partial_solution}{Can be either the result of a previous optimization
run (see \code{value}) or an (initial) \code{matrix} of species presences and
absences for each site in the landscape. The total number of presences must
match the estimated species richness of each site. If a result of a previous
optimization is used, its \code{optimized_grid} is used as initial matrix and
its \code{error} data frame will be extended with the new iterations.}

\item{fixed_species}{Fixed partial solution with species that are considered
as given. Those species are not going to be changed during optimization.}

\item{autostop}{The optimizer will stop after this number of iterations with no
improvement. Default: \code{0} means auto stop is disabled.}

\item{seed}{Seed for random number generator. Seed must be a positive integer value.
\code{seed = NA} means that a random integer is used as seed.}

\item{verbose}{If \code{TRUE} (default), a progress report is printed during
the optimization run.}

\item{interruptible}{Allow a run to be interrupted before completion.
\code{FALSE} increases the performance.#'}
}
\value{
A species presence/absence \code{matrix} of the study landscape.
}
\description{
Generate an optimized estimate of community composition
(species presences and absences) for every site in the study area.
}
\details{
\code{run_optimization_min_conf} is the core function of the
\code{spectre} package. The underlying algorithm of this function is
adapted from Mokany et al. (2011). A pairwise commonness matrix (having the
same structure as the \code{target} matrix) is calculated from the
\code{partial_solution} matrix and the value difference with the
\code{target} determined. If a difference is present and depending on the
set stopping criteria the algorithm continues. A random site in the
presence/absence matrix is selected, and a random presence record at this
site replaced with an absence. Every absence in the selected site is then
individually flipped to a presence and the value difference with the
objective recorded. The presence record which resulted in the lowest value
difference (minimum conflict) is retained. This cycle continues, with a
random site selected every iteration, until the pairwise commonness and
objective matrices match or the algorithm runs beyond the
\code{max_iterations}.
}
\references{
Mokany, K., Harwood, T.D., Overton, J.M., Barker, G.M., &
Ferrier, S. (2011). Combining \eqn{\alpha} and \eqn{\beta} diversity models to fill
gaps in our knowledge of biodiversity. Ecology Letters, 14(10), 1043-1051.
}
