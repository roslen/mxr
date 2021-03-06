% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mxr_create_annovar_inputs.R
\name{mxr_create_annovar_inputs}
\alias{mxr_create_annovar_inputs}
\title{Create Input File for Annovar.}
\usage{
mxr_create_annovar_inputs(target_list = "", alleles_from_reference = "",
  out_prefix = "", verbose = FALSE)
}
\arguments{
\item{target_list}{The "target_list" file created by a call to
mxr_create_region_file().}

\item{alleles_from_reference}{The "alleles_from_reference" file created by a
call to mxr_create_region_file().}

\item{out_prefix}{The output prefix.}

\item{verbose}{(Optional) Show verbose output. (DEFAULT=FALSE)}
}
\value{
TRUE if the PLINK run completed successfully. FALSE, otherwise.
}
\description{
\code{mxr_create_annovar_inputs} generates the input files needed by annovar.
}
\details{
This function needs to be run before annotating the significant SNPs. It
creates three files: 1. the target_list defined from the
alleles_from_reference and the, 2. the avinput, and 3.
avinput.corrected_chroms files.
}

