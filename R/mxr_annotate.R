#' Annotate the significant SNPs.
#'
#' \code{mxr_annotate} annotates the significant SNPs using annovar.
#'
#' Before running this function, \code{mxr_create_annovar_inputs} needs to be
#' run first to generate the necessary input files. Note that it is important
#' that the chromosome designation matches exactly the way the chromosome number
#' is designated in the annotation file for the reference genome. For example,
#' "chr01" is not the same as "Chr1", or "Chr01", moreso "1".
#'
#' @param avinput Make sure that the chromosome designations match those in the
#'   reference genome annotation. (e.g. use the *.corrected_chroms file).
#' @param annovar_reference_db The complete path (with trailing "/") to the
#'   annovar local database containing the reference genome.
#' @param database_build_number The database build number of the annovar db.
#' @param aa_matrix_file The grantham.matrix file or any other amino acid
#'   substitution "relevance" file.
#' @param out_prefix The output prefix.
#' @param verbose (Optional) Show verbose output. (DEFAULT=FALSE)
#' @return TRUE if the PLINK run completed successfully. FALSE, otherwise.
#'
#' @export
mxr_annotate <- function(avinput = "",
                         annovar_reference_db = "~/bin/annovar/rice_msu7/",
                         database_build_number = "msu7b",
                         aa_matrix_file = "~/bin/annovar/rice_msu7/grantham.matrix",
                         out_prefix = "",
                         verbose = FALSE) {

   ANNOVAR <- findApplication("annotate_variation.pl")

   if (verbose) cat("Annotating...")
   system(paste(ANNOVAR, "-out", out_prefix,
                "-build", database_build_number,
                avinput, annovar_reference_db,
                "-aamatrixfile", aa_matrix_file,
                "&>", paste0(out_prefix,".log")
   ))
   if (verbose) cat("DONE.\n")

   return (TRUE)
}



