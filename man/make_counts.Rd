% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/load_and_clean.R
\name{make_counts}
\alias{make_counts}
\title{Find daily butterfly counts for each regioin}
\usage{
make_counts(di)
}
\arguments{
\item{di}{A dataframe produced by read_ufb_datasheet}
}
\value{
A dataframe with daily observation counts for each region.
This is appropriate for generating label columns when joined with
daily abundance estimates from distance software.  Note that the
two subregions for Baldy Chato are reversed in the output of distance
so they need to be manually corrected.
}
\description{
Find daily butterfly counts for each regioin
}
\examples{
tt <- read_ufb_datasheet("Distance Input 2015.xlsx")
make_counts(tt)
}

