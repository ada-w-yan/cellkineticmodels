% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stats_functions.R
\name{preferred_model_order}
\alias{preferred_model_order}
\title{short summary of preferred model order}
\usage{
preferred_model_order(compare_output)
}
\arguments{
\item{compare_output}{output of loo::loo_compare}
}
\value{
a string of form "model X > model Y > model Z", listing models
from most to least preferred
}
\description{
short summary of preferred model order
}
