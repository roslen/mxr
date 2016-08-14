#' Associate the genotype with the phenotype.
#'
#' \code{mxr_associate} outputs the results of the mixed linear model
#' association between the SNPs and the phenotype. It also outputs the REsidual
#' Maximun Likelihood (REML) estimates as well as the pseudoheritability
#' estimate. All calculations are done using EMMAX
#' \url{http://genome.sph.umich.edu/wiki/EMMAX}
#'
#' Before running this function, make sure to compute EMMAX to compute the
#' association of each SNP, coded as minor allele additive dosage (using
#' --recode A in PLINK), with the phenotype.
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
mxr_associate <- function(tped_prefix, pheno_file, kin_file,
                          out_prefix, cov_file="", verbose=F) {
   if (cov_file != "")
      not_found <- which(!file.exists(c(paste(tped_prefix,c("tped","tfam"),sep="."),
                                        pheno_file, kin_file, cov_file)))
   else
      not_found <- which(!file.exists(c(tped_prefix, pheno_file, kin_file)))

   # Check which parameter file is not found on disk
   if (length(missing)>0)
      stop(paste(paste(c(tped_prefix, pheno_file, kin_file)[not_found], collapse=","),
                 "not found."))

   OUTPUT_DIRECTORY <- dirname(out_prefix)
   OUTPUT_PREFIX <- basename(out_prefix)
   EMMAX2 <- findApplication("emmax2")

   # Run the association
   result = tryCatch({
      system(paste(EMMAX2, ifelse(verbose, "-v", ""), "-d 10",
                   "-t", tped_prefix,
                   "-p", pheno_file,
                   "-k", kin_file,
                   ifelse(cov_file != "", paste("-c", cov_file), ""),
                   "-o", out_prefix
      ))
   }, warning = function(w) {
      return ("")
   }, error = function(e) {
      stop("")
   }, finally = {
      #cleanup-code
   })

   # If execution managed to reach this line, then everything went well.
   return (TRUE)
}


