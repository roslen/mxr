#' Create the Haploview Input Files.
#'
#' \code{mxr_create_haploview_inputs} creates the input files necessary to run
#' Haploview.
#'
#' The output of this file is required to run Haploview. It recodes in haploview
#' format the plink file, and then creates the batch input file that would
#' subsequently be fed into haploview.
#'
#' @param genotype_prefix The fullpath and prefix of the bim/bed/fam plink
#'   files.
#' @param snps_list The list of clumped SNPs.
#' @param out_prefix Path and prefix of the output files.
#' @param verbose (Optional) Show verbose output. (DEFAULT=FALSE)
#' @return TRUE if the PLINK run completed successfully. FALSE, otherwise.
#'
#' @export
mxr_create_haploview_inputs <- function(genotype_prefix, snps_list = "",
                                        out_prefix = "", verbose = FALSE) {

   if (!file.exists(paste0(out_prefix,".clumped"))) stop(paste0(out_prefix,
                                                                ".clumped",
                                                                " does not exist."))

   OUTPUT_DIRECTORY <- dirname(out_prefix)
   # OUTPUT_PREFIX <- basename(out_prefix)
   # SED <- findApplication("sed")
   # TAIL <- findApplication("tail")
   # XARGS <- findApplication("xargs")
   # UNIQ <- findApplication("uniq")
   # TABIX <- findApplication("tabix")

   PLINK2 <- findApplication("plink2")
   FIND <- findApplication("find")
   SORT <- findApplication("sort")
   PASTE <- findApplication("paste")
   AWK <- findApplication("awk")


   # Create haploview input files in preparation for generating the LD plots.
   if (verbose) cat("Recoding genotype file into haploview format...")
   system(paste(PLINK2, "--bfile", genotype_prefix,
                "--extract", snps_list,
                "--recode HV",
                "--out", paste0(out_prefix,".hv"),
                ifelse(!verbose,
                       paste("2>",paste0(out_prefix,".hv.errlog"),"1> /dev/null"),
                       "")
   ))
   if (verbose) cat("DONE.\n")


   # generate haploview batch input files
   if (verbose) cat("Creating haploview input files...")
   system(paste(FIND, OUTPUT_DIRECTORY,
                "-maxdepth 1",
                "-name '*.ped'",
                "-o -name '*.info'", "|",
                SORT, "|",
                PASTE, "-s -d ' \\n'", "|",
                AWK, "'{print $2 \"\\t\" $1;}'",
                ">", paste0(out_prefix,".hv.batch")
                ))
   if (verbose) cat("DONE.\n")

   # If execution managed to reach this line, then everything went well.
   return (TRUE)
}

