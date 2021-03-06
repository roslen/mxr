#' Post-GWAS Processing of Significant Association Peaks.
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
#' @param emmax2_results The fullpath and filename of the results of
#'   mxr_associate.
#' @param genotype_prefix The fullpath and prefix of the bim/bed/fam plink
#'   files.
#' @param snp_field The column header for the SNP column in the results file
#'   from mxr_associate. (Default: "snp_id")
#' @param clump_field The column header for the p-value column in the results
#'   file from mxr_associate. (Default: "p")
#' @param clump_p1 The maximum p-value of the index variants that define each
#'   clump. Index variants need to have values between 0 and this p-value.
#'   (Default: 0.00001)
#' @param clump_p2 The maximum association p-value of SNPs that belonging to the
#'   clump formed by an index variant. (Default: 0.01)
#' @param clump_kb The maximum distance of a clumped SNP from its index variant.
#'   (Default: 200)
#' @param clump_r2 The r-squared value of the clumped SNP with its index
#'   variant. (Default: 0.50)
#' @param clump_range_border_kb The number of kilobases past the region
#'   clump-range region. (Default: 20)
#' @param gene_list The list of genes must have the following columns: chromosome,
#'   start (bp), end (bp), locus_id. These fields must be tab separated and have an appropriately named column header. The chromosome must be coded to correspond to the plink data.
#'   @param out_prefix Path and prefix of the output files.
#' @param verbose (Optional) Show verbose output. (DEFAULT=FALSE)
#' @return TRUE if the PLINK run completed successfully. FALSE, otherwise.
#'
#' @export
mxr_clump <- function(emmax2_results,
                      genotype_prefix,
                      snp_field = "snp_id",
                      clump_field = "p",
                      clump_p1 = 0.00001,
                      clump_p2 = 0.01,
                      clump_kb = 200,
                      clump_r2 = 0.80,
                      clump_range_border_kb = 20,
                      gene_list = "",
                      out_prefix = "",
                      verbose = FALSE) {

   # Get the indices of input files that are not present
   not_found <- which(!file.exists(c(emmax2_results,
                                     paste(genotype_prefix,c("bim","bed","fam"),sep="."))))

   # If any of the input files are not present, then stop gracefully.
   if (length(not_found)>0) {
      stop(paste(paste(c(emmax2_results,
                         paste(genotype_prefix,c("bim","bed","fam"),sep="."))[not_found],
                       collapse=","),
                 "not found."))
   }

   OUTPUT_DIRECTORY <- dirname(out_prefix)
   OUTPUT_PREFIX <- basename(out_prefix)
   PLINK2 <- findApplication("plink2")

   # check if PLINK2 is present
   if (PLINK2=="") stop("plink2 is not found.")

   # Run the clump
   result = tryCatch({
      system(paste(PLINK2,
                   "--clump", emmax2_results,
                   "--bfile", genotype_prefix,
                   "--make-founders",
                   "--clump-snp-field", snp_field,
                   "--clump-field", clump_field,
                   "--clump-p1", clump_p1,
                   "--clump-p2", clump_p2,
                   "--clump-kb", clump_kb,
                   "--clump-r2", clump_r2,
                   "--r2 dprime",
                   ifelse(gene_list!="",
                          paste("--clump-range", gene_list,
                                "--clump-range-border", clump_range_border_kb),
                          ""),
                   "--out", out_prefix,
                   ifelse(!verbose,
                          paste("2>",paste0(out_prefix,".hv.errlog"),"1> /dev/null"),
                          "")
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
