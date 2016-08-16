#' Generate the Haploview Plot.
#'
#' \code{mxr_create_haploview_plot} generates the haploview plot using Gabriel's
#' algorithm.
#'
#' For really large LD plots, it may be necessary to set the minimum memory
#' footprint to greater than 2GB. The default memory setting for this function
#' was already increased to 2GB from the default Haploview setting of 512MB.
#'
#' @param path_to_haploview The fullpath to the Haploview.jar.
#' @param ped_file The fullpath to the PED file in haploview format.
#' @param info_file The fullpath to the corresponding INFO file of the PED file.
#' @param max_distance The window size in kb of the LD block.
#' @param memory_mb The memory to set aside for the haploview run.
#'   (DEFAULT:2048MB)
#' @param out_prefix The output prefix for the generated files.
#' @param verbose (Optional) Show verbose output. (DEFAULT=FALSE)
#' @return TRUE if the PLINK run completed successfully. FALSE, otherwise.
#'
#' @export
mxr_create_haploview_plot <- function(path_to_haploview="",
                                      ped_file = "",
                                      info_file = "",
                                      max_distance = 200,
                                      memory_mb = 2048,
                                      out_prefix = "",
                                      verbose = FALSE) {

   # if (!file.exists(paste0(out_prefix,".clumped"))) stop(paste0(out_prefix,
   #                                                              ".clumped",
   #                                                              " does not exist."))

   JAVA <- findApplication("java")
   HAPLOVIEW <- path_to_haploview
   system(paste(JAVA, "-jar",
                HAPLOVIEW, "-n",
                "-pedfile", ped_file,
                "-info", info_file,
                "-skipcheck",
                "-blockoutput GAB",
                "-dprime",
                "-out", out_prefix,
                "-svg",
                "-png",
                "-pairwiseTagging",
                "-maxDistance", max_distance
                ))

   # If execution managed to reach this line, then everything went well.
   return (TRUE)
}
