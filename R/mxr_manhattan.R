# This is a utility function
qq <- function (pvector, ...)
{
   if (!is.numeric(pvector)) stop("Input must be numeric.")
   #pvector <- pvector[!is.na(pvector) & !is.nan(pvector) & !is.null(pvector) &
   #                      is.finite(pvector) & pvector < 1 & pvector > 0]
   o = -log10(sort(pvector, decreasing = FALSE))
   e = -log10(ppoints(length(pvector)))
   def_args <- list(pch = 20,
                    xlim = c(0, max(e)),
                    ylim = c(0, max(o)),
                    xlab = expression(Expected ~ ~-log[10](italic(p))),
                    ylab = expression(Observed ~ ~-log[10](italic(p))))
   dotargs <- list(...)
   tryCatch(
      do.call("plot",
              c(list(x = e, y = o),
                def_args[!names(def_args) %in% names(dotargs)],
                dotargs)),
      warn = stop)
   abline(0, 1, col = "red")
} # qq


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
#' @param xlim Maximum x-axis value (for the qq plot).
#' @param ylim Maximum y-axis value.
#' @param file_prefix Fullpath and prefix for manhattan and qq plot file names.
#' @param width Width (in inches) of the figure.
#' @param height Height (in inches) of the figure.
#' @param unit Default: "in" (inches).
#' @param res Resolution of the image.
#' @param mar Vector of margins.
#' @param mgp The mgp vector in plot.
#' @param manhattan Set to TRUE to generate Manhattan plot.
#' @param qq Set to TRUE to generate QQ plot.
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
                      xlim = 10,
                      ylim = 20,
                      file_prefix = "",
                      width = 5,
                      height = 5,
                      unit = "in",
                      res = 300,
                      mar = c(5,4,4,0)+.1,
                      mgp = c(3,.8,0),
                      cex.axis = 0.8,
                      main = "",
                      y_tickmarks = "",
                      y_tickmark_labels = "",
                      manhattan = T,
                      qq = T) {


   if (!file.exists(gwas_results)) stop(paste(gwas_results, "does not exist."))

   # Load the data
   cat(paste0("Loading GWAS data..."))
   data <- data.table::fread(gwas_results, header = T, sep = "\t", stringsAsFactors = F,
                             data.table = F)
   cat("DONE.\n")


   # Set parameters based on the data
   genomewideline <- genomewide_line
   suggestiveline <- suggestive_line
   highlight <- highlight

   if (manhattan) {
      cat("Generating manhattan plot...")
      png(filename = paste0(file_prefix,".manhattan.png"),
          width = width, height = height, units = unit, res=res)
      par(mar=c(5,4,4,0)+.1, mgp=c(3,.8,0))
      qqman::manhattan(data,
                       chr = chr,
                       bp = bp,
                       p = p,
                       snp = snp,
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


   #par(mar = mar, mgp = mgp)
   if (qq) {
      # Calculate inflation factor
      # Reference: http://genometoolbox.blogspot.com/2014/08/how-to-calculate-genomic-inflation.html
      chisq <- qchisq(1 - data[, c(p)], 1)
      inflation_factor <- median(chisq)/qchisq(0.5,1)

      cat(paste0("Generating qq plot..."))
      png(filename = paste0(file_prefix,".qq.png"),
         width = width, height = height, units = unit, res=res)
      par(mar=c(5,4,4,0)+.1, mgp=c(3,.8,0))
      qq(data[, c(p)],
                xlim = c(0, xlim), ylim = c(0, ylim), las=1, bty="l",
                cex.axis = cex.axis, yaxt="n",
                main = main, xlab = "", ylab = "")
      # abline(0,1)
      axis(side=2, at = y_tickmarks, labels = y_tickmark_labels,
           las=1, mgp = mgp, cex.axis = cex.axis)
      mtext(expression(Expected ~ ~-log[10](italic(p))), side=1, line=2.5)
      mtext(expression(Observed ~ ~-log[10](italic(p))), side=2, line=2.5)
      mtext(paste0("Inflation factor: ", sprintf("%1.2f",inflation_factor)), side=1, line=3.5)

      res <- dev.off()
      cat("DONE.\n")
   }
}


