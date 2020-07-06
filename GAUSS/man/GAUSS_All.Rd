\name{GAUSS_All}
\alias{GAUSS_All}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Run GAUSS association test for given summary and GMT files (in parallel).
}
\description{
Run and output GAUSS p-values
}
\usage{
GAUSS_All(summary_file, gene_name, pv_name, output_file, gmt, verbose, ags, parallel, jobs, jobfile, start, stop)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{summary_file}{
Summary file containing the gene-names and corresponding p-values
}
  \item{gene_name}{
  Column number for gene names in summary file
}
  \item{pv_name}{
  Column number for p-value in summary file
}
  \item{output_file}{
  prefix for output file names
}
  \item{gmt}{
  GMT file containing the list of gene-set, one gene-set per line
}
  \item{verbose}{
  Print extra output; default = TRUE
}
  \item{ags}{
  control arguments for running GAUSS
}
  \item{path}{
  path to files for GAUSS
}
  \item{parallel}{
  logical indicating whether parallel jobs are to be created; default = FALSE
}
  \item{jobs}{
  number of parallel jobs to be created if parallel = TRUE
}
  \item{jobfile}{
  name of the output file with parallel jobs for GAUSS if parallel = TRUE
}
\item{start}{
  starting point of the GMT file if parallel = TRUE
}
\item{stop}{
  stopping point of the GMT file if parallel = TRUE
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The pvalues for the GAUSS association test is stored in $output_file.out and a log file is generated at $output_file.log
}
\references{

}
\author{
Diptavo Dutta
}
