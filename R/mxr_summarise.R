#' Summarize the GWAS results into an Excel file.
#'
#' \code{mxr_summarise} creates an Excel workbook with worksheets for the SNPs
#' and linked SNPs, SNPs and associated genes, annotation of SNPs, annotation of
#' genes, and the identified tag SNPs. A worksheet is also prepared that
#' contains the parameters used in the GWAS run.
#'
#' \code{mxr_summarise} requires the results of \code{mxr::clump}.
#'
#' @param trait The directory containing the trait phenotype. All GWAS results
#'   are also stored into this directory that includes the outputs of this
#'   function.
#' @param clumped_genes The clumped.genes file outputted by the plink file.
#' @param reference_annotation The excel file containing the annotation file
#'   that will be used to annotate the genes.
#' @param out_prefix Path and prefix of the output files.
#' @param verbose (Optional) Show verbose output. (DEFAULT=FALSE)
#' @return TRUE if the PLINK run completed successfully. FALSE, otherwise.
#'
#' @export
mxr_summarise <- function(trait = "",
                          clump_file = "",
                          clumped_genes = "",
                          reference_annotation = "",
                          output_prefix ="",
                          variant_function = "",
                          exonic_variant_function = "",
                          verbose = FALSE) {
   if (!dir.exists(trait)) stop(paste(trait, "directory does not exist."))

   OUTPUT_DIRECTORY <- dirname(output_prefix)
   OUTPUT_PREFIX <- basename(output_prefix)

   TRAIT <- trait
   PATH <- TRAIT

   cat("Retrieving list of genes...")
   genes <- data.table::fread(clumped_genes, header=F, stringsAsFactors=F, data.table=F)
   cat("DONE.\n")

   cat("Retrieving functional annotation of genes...")
   annot <- readxl::read_excel(reference_annotation)
   cat("DONE.\n")


   # Get the annotation of the genes from the annotation file
   cat("Saving annotated genes to disk...")
   res <- dplyr::inner_join(genes, annot, by=c("V1"="GENE_ID"))
   colnames(res)[1:2] <- c("LOCUS","CHROM")
   write.table(res, file=paste0(clumped_genes,".annotation"),
               append=F, quote=F, sep="\t", col.names=T, row.names=F)
   cat("DONE.\n")



   cat("Loading functional annotations of significant SNPs...")
   var <- data.table::fread(variant_function,
                            data.table = F,
                            header = F, sep = "\t", stringsAsFactors = F)
   exn <- tryCatch(data.table::fread(exonic_variant_function,
                                     data.table = F,
                                     colClasses = c("character", "character", "character",
                                                    "character", "numeric", "numeric",
                                                    "character", "character", "character"),
                                     header = F, sep = "\t", stringsAsFactors = F),
                   error=function(e) NULL)
   cat("DONE.\n")

   if (!is.null(exn)) {
      cat("Consolidating annotations...")
      res <- dplyr::left_join(var,
                              exn,
                              by = c("V3"="V4", "V4"="V5", "V5"="V6", "V6"="V7", "V7"="V8"))
      write.table(res, file = paste0(variant_function,".consolidated"),
                  append = F, quote = F,
                  sep = "\t", row.names = F, col.names = F)
      cat("DONE.\n")
   }



   ## Portion to create the excel workbook
   Sys.setenv(R_ZIPCMD="/usr/bin/zip")   # without this openxlsx will cry




   return (TRUE)
}
