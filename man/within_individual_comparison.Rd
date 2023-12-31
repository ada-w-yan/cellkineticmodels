% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stats_functions.R
\name{within_individual_comparison}
\alias{within_individual_comparison}
\title{compute 95\% CIs for ratios between lymphocyte parameters between cell populations}
\usage{
within_individual_comparison(fit)
}
\arguments{
\item{fit}{output of fit_all_stan}
}
\value{
named list with three elements: pb_w, dstar and delay.
Each element is a named list with two elements: quantiles and p_value.
quantiles is a matrix with 3 rows and p(p+1)/2 columns, where p is the number of cell
populations for the individual in question.  Eech column shows the 2.5th percentile,
50th percentile and 97.5th percentile for the ratio between parameter values for
pairs of cell populations.
p_value is a vector of length p(p+1)/2.  Each element is 2 * min(q, 1 - q), where q is the
proportion of draws from the posterior whose ratios are less than 1, for a particular
pair of cell populations.
}
\description{
compute 95\% CIs for ratios between lymphocyte parameters between cell populations
}
