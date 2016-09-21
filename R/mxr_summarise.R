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
#' @param gwas_result The output of emmax2. (*.ps)
#' @param genotype_bim The .bim file of the genotype data.
#' @param clump_file The .clumped file.
#' @param clump_ranges_file The .clumped.ranges file.
#' @param clumped_genes The .clumped.genes file outputted by the plink file.
#' @param clumped_snps The .clumped.snps file.
#' @param reference_annotation The excel file containing the annotation file
#'   that will be used to annotate the genes.
#' @param variant_function Output of annovar.
#' @param exonic_variant_function Output of annovar.
#' @param output_prefix Path and prefix of the output files.
#' @param verbose (Optional) Show verbose output. (DEFAULT=FALSE)
#' @return TRUE if the PLINK run completed successfully. FALSE, otherwise.
#'
#' @export
mxr_summarise <- function(trait = "",
                          gwas_result = "",
                          genotype_bim = "",
                          clump_file = "",
                          clump_ranges_file = "",
                          clumped_genes = "",
                          clumped_snps = "",
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

   wb <- openxlsx::createWorkbook()
   openxlsx::modifyBaseFont(wb, fontSize=10)


   ### 2. Create the worksheets under this workbook

   # *.clumped
   openxlsx::addWorksheet(wb, sheetName="SNPs and linked SNPs")

   # *.ranges
   openxlsx::addWorksheet(wb, sheetName="SNPs and associated genes")

   # *.variant_function.consolidated
   # Without header
   openxlsx::addWorksheet(wb, sheetName="Annotation of SNPs")

   # *.clumped.ranges.genes.annotation
   openxlsx::addWorksheet(wb, sheetName="annotation of genes")
   #cat("DONE.\n")

   # create the tag SNPs
   openxlsx::addWorksheet(wb, sheetName="tag SNPs")


   ### 3. Write the data into the worksheets
   # sheet 1
   cat("Compiling SNPs and linked SNPs...")
   data <- read.table(clump_file, header = T, stringsAsFactors = F)
   openxlsx::writeData(wb = wb,
                       sheet = "SNPs and linked SNPs",
                       x = data,
                       colNames = TRUE,
                       rowNames = FALSE)
   openxlsx::setColWidths(wb,
                          sheet="SNPs and linked SNPs",
                          cols=1:(dim(data)[2]), widths="auto")
   cat("DONE.\n")

   # sheet 2
   cat("Compiling SNPs and associated genes...")
   data <- read.table(clump_ranges_file, header = T, stringsAsFactors = F)
   openxlsx::writeData(wb = wb,
                       sheet = "SNPs and associated genes",
                       x = data,
                       colNames = TRUE,
                       rowNames = FALSE)
   openxlsx::setColWidths(wb,
                          sheet="SNPs and associated genes",
                          cols=1:(dim(data)[2]), widths="auto")
   cat("DONE.\n")


   # sheet 3
   cat("Reading SNP annotation files...")
   data <- tryCatch(read.table(paste0(variant_function,".consolidated"),
                               header = F, sep = "\t", stringsAsFactors = F),
                    error=function(e) NULL)
   if (!is.null(data)) {
      data <- data[, -c(12)] # remove the duplicate "Notes" field
      colnames(data) <- c("snp_location", "locus_id", "chrom", "start", "end",
                          "ref", "alt", "notes", "line_no", "type_of_aa_change",
                          "transcript_and_aa_change")
   } else {
      # Read only the variant function since there were no exonic SNPs that were found.
      data <- tryCatch(read.table(variant_function,
                                  header = F, sep = "\t", stringsAsFactors = F),
                       error = function(e) NULL)
      colnames(data) <- c("snp_location", "locus_id",
                          "chrom", "start", "end", "ref", "alt", "notes")
   }
   cat("DONE.\n")


   cat("Converting chromosome codes to GWAS compatible format...")
   data$chrom <- gsub("Chr","", data$chrom)
   cat("DONE.\n")

   cat("Creating region file for SNP extraction...")
   write.table(data[, c("chrom", "start")],
               file = paste0(clumped_snps, ".region"),
               quote = F, sep = "\t", append = F, row.names = F, col.names = F)
   cat("DONE.\n")



   # cat("Putting bim details in to the gwas results...")
   # system(paste0(
   #    "tail -n+2 ", gwas_result,
   #    " | paste - ", genotype_bim,
   #    " > ", gwas_result, ".bim"
   # ))
   # cat("DONE.\n")


   cat("Indexing the detailed GWAS results...")
   system(paste0("bgzip -c ", gwas_result, ".bim > ",
                 gwas_result, ".bim.gz && ",
                 "tabix -s5 -b8 -e8 ", gwas_result, ".bim.gz"))
   cat("DONE.\n")


   # Extract the GWAS results that are in the clumped snps file
   cat("Extracting details of clumped SNPs...")
   system(paste0("tabix -s5 -b8 -e8 ", gwas_result, ".bim.gz -R ",
                 clumped_snps, ".region > ",
                 clumped_snps,".details"))


   ####>>>>>STOPPED HERE!!!<<<<<<####

   snp_details <- read.table(paste0(clumped_snps, ".details"),
                             sep = "\t", stringsAsFactors = F,
                             colClasses = c("character", "numeric", "numeric", "numeric",
                                            "character", "character", "numeric", "numeric",
                                            "character", "character"),
                             col.names=c("snp_id", "allele_effect", "std_error", "pval",
                                         "chrom", "snp_id", "cm", "bp", "a1", "a2"))
   cat("DONE.\n")

   # Add the snp details to the annotation data
   data <- dplyr::inner_join(data, snp_details[,c("snp_id",
                                                  "chrom",
                                                  "bp",
                                                  "allele_effect",
                                                  "pval")],
                             by=c("chrom"="chrom", "start"="bp"))
   # TODO: create friendly names for data

   # Writing the aggregated data to disk
   openxlsx::writeData(wb=wb, sheet="Annotation of SNPs", x=data, colNames=TRUE, rowNames=FALSE)
   openxlsx::setColWidths(wb, sheet="Annotation of SNPs", cols=1:(dim(data)[2]), widths="auto")


   # sheet 4
   # NOTE: Putting 'quote=""' ensures that any quotation marks ", or ' in any field is ignored.
   cat("Compiling gene annotations...")
   data <- read.table(paste0(clumped_genes, ".annotation"),
                      header = T, sep = "\t", stringsAsFactors = F, quote = "")
   openxlsx::writeData(wb = wb,
                       sheet = "annotation of genes",
                       x=data,
                       colNames=TRUE, rowNames=FALSE)
   openxlsx::setColWidths(wb, sheet="annotation of genes", cols=1:(dim(data)[2]), widths="auto")
   cat("DONE.\n")


   # sheet 5
   # save the TAG SNPs here
   cat("Compiling tag SNPs...")
   tag_snps <- data.frame(tag_snp = character(0), captured_alleles = character(0))
   tag_files <- list.files(path = PATH, pattern = "*.TAGS$", include.dirs = T)
   for (i in 1:length(tag_files)) {
      # Get the contents of the entire TAGS file into a vector
      all_data <- readLines(paste0(PATH, "/", tag_files[i]))

      # No of alleles is in line 1, second word
      alleles <- as.numeric(strsplit(all_data[1], " ")[[1]][2]) # numbers are hardcoded
      snps <- as.numeric(strsplit(all_data[3], " ")[[1]][2]) # numbers are also hardcoded
      for (j in 1:snps) {
         res <- strsplit(all_data[alleles+7+(j-1)], split="\t")
         tag_snps <- rbind(tag_snps,
                           data.frame(tag_snp=res[[1]][1], captured_alleles=res[[1]][2],
                                      stringsAsFactors=F))
      }
   }
   openxlsx::writeData(wb=wb, sheet="tag SNPs", x = tag_snps, colNames = TRUE, rowNames = FALSE)
   openxlsx::setColWidths(wb, sheet="tag SNPs", cols=1:(dim(tag_snps)[2]), widths="auto")
   cat("DONE.\n")


   # Save the workbook
   cat("Saving the summary...")
   openxlsx::saveWorkbook(wb, paste0(output_prefix, ".summary.xlsx"), overwrite = TRUE)
   cat("DONE.\n")


   return (TRUE)
}
