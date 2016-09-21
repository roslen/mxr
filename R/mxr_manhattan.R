#' Generate the manhattan plot.
#'
#' \code{mxr_manhattan} generates the manhattan plot of the GWAS results.
#'
#' Map the needed columns in the GWAS results to the manhattan function.
#'
#' @param gwas_results The fullpath to the GWAS results. The file must have a
#'   header line indicating the snp_id, chromosome, pvalue, and base-pair
#'   position.
#' @param chr Header name of chromosome column.
#' @param bp Header for base-pair position.
#' @param p Header for the p-value column.
#' @param snp Header for the SNP ID column.
#' @param suggestive_line Default -log10(1e-5). Set to FALSE to disable.
#' @param genomewideline Default -log10(1e-7). Set to FALSE to disable.
#' @param highlight Character vector of SNPs to highlight in the manhattan plot.
#' @param ylim Maximum y-axis value.
#' @param file Fullpath to the where to store the filename.
#' @param width Width (in inches) of the figure.
#' @param height Height (in inches) of the figure.
#' @param unit Default: "in" (inches).
#' @param res Resolution of the image.
#' @param mar Vector of margins.
#' @param mgp The mgp vector in plot.
#'
#' @export
mxr_manhattan <- function(gwas_results = "",
                      chr = "CHR",
                      bp = "BP",
                      p = "P",
                      snp = "SNP",
                      suggestive_line = -log10(1e-5),
                      genomewide_line = -log10(1e-7),
                      highlight = NULL,
                      ylim = 20,
                      file = "",
                      width = 5,
                      height = 5,
                      unit = "in",
                      res = 300,
                      mar = c(5,4,4,0)+.1,
                      mgp = c(3,.8,0),
                      cex.axis = 0.8,
                      main = "",
                      y_tickmarks = "",
                      y_tickmark_labels = "") {


   if (!file.exists(gwas_results)) stop(paste(gwas_results, "does not exist."))

   # Load the data
   cat(paste0("Loading GWAS data..."))
   data <- data.table::fread(gwas_results, header = T, sep = "\t", stringsAsFactors = F,
                             data.table = F)
   cat("DONE.\n")


   # Set parameters based on the data
   genomewideline <- genomewide_line
   suggestiveline <- suggestive_line
   ylim <-ylim
   xlim <-ylim

   highlight <- highlight

   cat("Generating manhattan plot...")
   png(filename = file, width = width, height = height, units = unit, res=res)
   par(mar = mar, mgp = mgp)
   manhattan(data,
             highlight=highlight,
             ylim=c(0,ylim), las=1, bty="l", cex.axis = cex.axis,
             yaxt = "n", xlab = "", ylab = "",
             main = main,
             suggestiveline = suggestiveline, genomewideline = genomewideline)

   axis(side=2, at = y_tickmarks, labels = y_tickmark_labels, las=1,
        mgp = mgp, cex.axis = cex.axis)

   mtext("Chromosome", side=1, line=2.5)
   mtext(expression(-log[10](italic(p))), side=2, line=2.5)

   res <- dev.off()
   cat("DONE.\n")

}
