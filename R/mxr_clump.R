#' Post-GWAS Processing of Significant Peak.
#'
#' \code{mxr_clump} determines the clumps of SNPs based on LD along the region
#' that registers a significant association signal with the phenotype. This uses
#' the \code{--clump} command of plink.
#' \url{https://www.cog-genomics.org/plink2}
#'
#' Before running this function, make sure that the mxr_associate function has
#' already been called and that the results are successfully obtained without
#' warnings (that imply fatal errors) or errors.
#'
#' @param tped_prefix The path to the transposed '12' recoded plink file. The A1
#'   allele (the minor allele -- or effect allele) is recoded as "1", while A2
#'   allele is recoded as "2". The complete call to plink is \code{'plink2
#'   --bfile [bed_prefix] --recode12 --output-missing-genotype 0 --transpose
#'   --out [tped_prefix]'}.
#' @param pheno_file The path to the phenotype file. There are three required
#'   columns: FID (family ID), IID (individual ID), and phenotype column. It is
#'   recommended that the phenotype column be approximately normally
#'   distributed, in which case, the warpedLMM software would be very useful.
#'   IMPORTANT: The samples need to be sorted according to how they are arranged
#'   in the genotype file. Missing phenotype values must be represented as 'NA'.
#' @param kin_file The path to the kinship matrix file computed using emmax-kin.
#'   Use the Balding-Nichols algorithm. Generate the kinship matrix by a call to
#'   emmax-kin \code{'emmax-kin-intel64 -v -d 10 [tped_prefix]'}.
#' @param out_prefix Path and prefix of the output files.
#' @param cov_file (Optional) Path to the covariate file. The first two columns
#'   of the covariate file must be the FID and IID of the samples. The third
#'   column must be a vector of '1'. Quantitative covariates start from the
#'   fourth column onward.
#' @param verbose (Optional) Show verbose output. (DEFAULT=FALSE)
#' @return TRUE if the EMMAX run completed successfully. FALSE, otherwise.
#'
#' @export
mxr_clump <- function(tped_prefix, pheno_file, kin_file,
                          out_prefix, cov_file="", verbose=F) {
}
