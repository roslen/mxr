#' Create the Region File from the Clumps.
#'
#' \code{mxr_create_region_file} creates the region file that will be used to
#' extract the corresponding information from SNPs stored in the SNP reference
#' file (that contains the reference alleles collated from processing the VCF
#' files.)
#'
#' The output of this file will be used in annotating the GWAS results.
#'
#' @param genotype_prefix The fullpath and prefix of the bim/bed/fam plink
#'   files.
#' @param reference_alleles The genome-wide list of reference allele for each
#'   SNP. It has three mandatory columns: chr, bp, Ref_allele. Chr column may
#'   have values like "chr01"..."chr12", bp are just the base pair positions of
#'   the SNPs, and the Ref_allele column is formatted as "Ref=<nuc>".
#' @param out_prefix Path and prefix of the output files.
#' @param verbose (Optional) Show verbose output. (DEFAULT=FALSE)
#' @return TRUE if the PLINK run completed successfully. FALSE, otherwise.
#'
#' @export
mxr_create_region_file <- function(genotype_prefix, reference_alleles = "",
                                   out_prefix = "", verbose = FALSE) {

   if (!file.exists(paste0(out_prefix,".clumped"))) stop(paste0(out_prefix,
                                                                ".clumped",
                                                                " does not exist."))

   OUTPUT_DIRECTORY <- dirname(out_prefix)
   OUTPUT_PREFIX <- basename(out_prefix)
   SED <- findApplication("sed")
   TAIL <- findApplication("tail")
   AWK <- findApplication("awk")
   XARGS <- findApplication("xargs")
   SORT <- findApplication("sort")
   UNIQ <- findApplication("uniq")
   TABIX <- findApplication("tabix")


   # Create the REGION_FILE
   if (verbose) cat("Extracting target SNPs from plink.bim file...")
   snps <- read.table(paste0(out_prefix,".clumped.snps"),
                      header = F, stringsAsFactors = F, sep = "\t")
   bim <- data.table::fread(paste(genotype_prefix,"bim",sep="."),
                            data.table = F, stringsAsFactors = F)
   snps_bim <- dplyr::inner_join(snps, bim, by=c("V1"="V2"))
   rm(bim)

   # Make this the same format as the reference file
   snps_bim[,2] <- paste("chr", sprintf("%02d", snps_bim[,2]), sep="")

   # Write the matching records to disk
   write.table(snps_bim[,c(2, 4, 1, 5:6)],
               file = paste0(out_prefix, ".clumped.snps.target_list"),
               append = F, quote = F, sep = "\t", col.names = F, row.names = F)
   if (verbose) cat("DONE.\n")


   if (verbose) cat("Extracting the corresponding regions in the genome...")
   system(paste(TABIX, "-s1 -b2 -e2", reference_alleles,
                "-R", paste0(out_prefix, ".clumped.snps.target_list"), "|",
                AWK, "'{chr=$1; sub(/chr/, \"\", chr); print chr \"\\t\" $1 \"\\t\" $2 \"\\t\" $3;}'", "|",
                SORT, "-k 1n,1 -k 3n,3", "|",
                AWK, "'{print $2 \"\\t\" $3 \"\\t\" $4;}'",
                ">" , paste0(out_prefix, ".clumped.snps.alleles_from_reference")
   ))
   if (verbose) cat("DONE.\n")


   # If execution managed to reach this line, then everything went well.
   return (TRUE)
}
