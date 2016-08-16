#' Extract SNPs from the Clumps.
#'
#' \code{mxr_extract_clumped_snps} creates a list of SNPs that were identified
#' to be part of the clumps.
#'
#' This function needs to be run after the \code{mxr_clump} function.
#'
#' @param out_prefix Path and prefix of the output files.
#' @param verbose (Optional) Show verbose output. (DEFAULT=FALSE)
#' @return TRUE if the PLINK run completed successfully. FALSE, otherwise.
#'
#' @export
mxr_extract_clumped_snps <- function(out_prefix = "", verbose = FALSE) {

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

   if (verbose) cat("Extracting the significant SNPs...")
   cmd_line <- paste(SED, "'/^$/d'", paste0(out_prefix,".clumped"), "|",
                     TAIL, "-n+2", "|",
                     AWK, "'{print $3 \",\" $12}'", "|",
                     SED, "'s/,NONE//g'", "|",
                     SED, "'s/(1)//g'", "|",
                     XARGS, "|",
                     SED, "-e 's/ /,/g'", "|",
                     SED, "'s/,/\\n/g'", "|",
                     AWK, "'{split($0, a,\"_\");
                               print $0 \"\\t\" a[2] \"\\t\" a[3];}'", "|",
                     SORT, "-k 2n,2 -k 3n,3", "|",
                     AWK, "'{print $1}'",
                     ">", paste0(out_prefix, ".clumped.snps")
   )
   # cat (paste0(cmdline,"\n"))
   system(cmd_line)
   if (verbose) cat("DONE.\n")

   # If execution managed to reach this line, then everything went well.
   return (TRUE)
}
