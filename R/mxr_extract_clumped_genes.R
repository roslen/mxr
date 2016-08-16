#' Extract Genes from the Clumps.
#'
#' \code{mxr_extract_clumped_genes} extracts the genes that overlap from the
#' clumps.
#'
#' This function needs to be run after the \code{mxr_clump} function.
#'
#' @param out_prefix Path and prefix of the output files.
#' @param verbose (Optional) Show verbose output. (DEFAULT=FALSE)
#' @return TRUE if the PLINK run completed successfully. FALSE, otherwise.
#'
#' @export
mxr_extract_clumped_genes <- function(out_prefix = "", verbose = FALSE) {

   if (!file.exists(paste0(out_prefix,".clumped.ranges"))) stop(paste0(out_prefix,
                                                                       ".clumped.ranges",
                                                                       " does not exist."))
   OUTPUT_DIRECTORY <- dirname(out_prefix)
   OUTPUT_PREFIX <- basename(out_prefix)
   PLINK2 <- findApplication("plink2")
   SED <- findApplication("sed")
   TAIL <- findApplication("tail")
   AWK <- findApplication("awk")
   XARGS <- findApplication("xargs")
   SORT <- findApplication("sort")
   UNIQ <- findApplication("uniq")
   TABIX <- findApplication("tabix")

   if (verbose) cat("Creating a list of the associated or nearby genes...")
   cmd_line <- paste(SED, "'/^$/d'", paste0(out_prefix,".clumped.ranges"), "|",
                     TAIL, "-n+2", "|",
                     AWK, "'{print  $7}'", "|",
                     SED, "'s/[][]//g'", "|",
                     XARGS, "|",
                     SED, "'s/ /,/g'", "|",
                     SED, "'s/,/\\n/g'", "|",
                     AWK, "'{loci=$0; sub(/LOC_Os/,\"\",loci);
                               split(loci,a,\"g\");
                               print $0 \"\\t\" a[1] \"\\t\" a[2]}'", "|",
                     SORT, "-k 2n,2 -k 3n,3", "|",
                     UNIQ, "|",
                     AWK, "'{print $1}'",
                     ">", paste0(out_prefix,".clumped.ranges.genes")
   )
   # cat(paste0(cmd_line, "\n"))
   system(cmd_line)
   if (verbose) cat("DONE.\n")

   # If execution managed to reach this line, then everything went well.
   return (TRUE)
}
