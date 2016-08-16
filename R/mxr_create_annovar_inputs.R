#' Create Input File for Annovar.
#'
#' \code{mxr_create_annovar_inputs} generates the input files needed by annovar.
#'
#' This function needs to be run before annotating the significant SNPs. It
#' creates three files: 1. the target_list defined from the
#' alleles_from_reference and the, 2. the avinput, and 3.
#' avinput.corrected_chroms files.
#'
#' @param target_list The "target_list" file created by a call to
#'   mxr_create_region_file().
#' @param alleles_from_reference The "alleles_from_reference" file created by a
#'   call to mxr_create_region_file().
#' @param out_prefix The output prefix.
#' @param verbose (Optional) Show verbose output. (DEFAULT=FALSE)
#' @return TRUE if the PLINK run completed successfully. FALSE, otherwise.
#'
#' @export
mxr_create_annovar_inputs <- function(target_list = "",
                                      alleles_from_reference = "",
                                      out_prefix = "",
                                      verbose = FALSE) {

   CAT <- findApplication("cat")
   AWK <- findApplication("awk")
   SED <- findApplication("sed")


   if (verbose) cat("Reading the alleles from reference...")
   list1 <- read.table(alleles_from_reference, header = FALSE,
                       quote = "", stringsAsFactors = FALSE)
   list1 <- unique(list1)
   if (verbose) cat("DONE.\n")

   if (verbose) cat("Reading the target list...")
   list2 <- read.table(target_list, header = FALSE,
                       quote = "", stringsAsFactors = FALSE)
   list2 <- unique(list2)
   if (verbose) cat("DONE.\n")

   if (verbose) cat("Matching the two lists...")
   list3 <- dplyr::inner_join(list1, list2, by=c("V1"="V1", "V2"="V2"))
   if (verbose) cat("DONE.\n")

   if (verbose) cat("Writing to disk...")
   write.table(list3,
               file=paste0(out_prefix, ".alleles_from_reference.target_list"),
               quote = FALSE, append = FALSE,
               col.names = FALSE, row.names = FALSE, sep = "\t")
   if (verbose) cat("DONE.\n")

   ###
   # Proceed with the rest of the script
   if (verbose) cat("Creating the annotation input file...")
   system(paste(CAT, paste0(out_prefix, ".alleles_from_reference.target_list"), "|",
                AWK, "'{split($3,a,\"=\"); if (a[2] != $6) print $1 \"\\t\" $2 \"\\t\" $2 \"\\t\" a[2] \"\\t\" $6 \"\\t\" \";A1 in .bim is reference allele.\"; else print $1 \"\\t\" $2 \"\\t\" $2 \"\\t\" a[2] \"\\t\" $5 \"\\t\" \"\"; }'", ">",
                paste0(out_prefix,".avinput")
   ))
   if (verbose) cat("DONE.\n")


   # IMPORTANT: Possibly reformat the "CHROM" field of the *.avinput file because the rice_msu7 format writes it as 'Chr#'
   if (verbose) cat("Reformating avinput file...")
   system(paste(CAT, paste0(out_prefix,".avinput"), "|",
                SED, "'s/chr0/Chr/g'", "|",
                SED, "'s/chr/Chr/g'", ">",
                paste0(out_prefix,".avinput",".corrected_chroms")
   ))
   if (verbose) cat("DONE.\n")

}



