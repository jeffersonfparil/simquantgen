% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim.R
\name{fn_simulate_genotypes}
\alias{fn_simulate_genotypes}
\title{Simulate genotype data with some linkage disequilibrium between neighbouring loci}
\usage{
fn_simulate_genotypes(
  n = 100,
  l = 500,
  ploidy = 2,
  n_alleles = 2,
  min_allele_freq = 0.01,
  n_chr = 5,
  max_pos = 1.35e+08,
  dist_bp_at_50perc_r2 = 5e+06,
  n_threads = 2,
  verbose = FALSE,
  show_correlation_heatmap = FALSE
)
}
\arguments{
\item{n}{number of samples (individuals or pools)}

\item{l}{number of loci}

\item{ploidy}{ploidy level of individual samples or the number of individuals multiplied by their ploidy to simulate pools}

\item{n_alleles}{number of alleles per locus}

\item{min_allele_freq}{minimum minor allele frequency}

\item{n_chr}{number of chromosomes}

\item{max_pos}{total length of the genome}

\item{dist_bp_at_50perc_r2}{distance in bases at which the estimated linkage between loci at both ends is at 50\%}

\item{n_threads}{number of computing cores or threads to use in parallel simulation of genotypes}

\item{verbose}{show simulation messages?}

\item{show_correlation_heatmap}{show correlation heatmap?}
}
\value{
genotype matrix with $n$ rows, $l x (n_alleles-1)$ columns, and named rows and columns
}
\description{
Simulate genotype data with some linkage disequilibrium between neighbouring loci
}
\examples{
G_1 = fn_simulate_genotypes()
G_2 = fn_simulate_genotypes(n=150, l=10000, ploidy=4, n_alleles=3, verbose=TRUE)
G_3 = fn_simulate_genotypes(n=200, l=2e5, ploidy=4, n_alleles=3, n_chr=10, max_pos=2.2e9, dist_bp_at_50perc_r2=5e6, n_threads=32, verbose=TRUE)
}
