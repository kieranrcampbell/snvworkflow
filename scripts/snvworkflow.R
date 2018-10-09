
suppressPackageStartupMessages({
  library(tidyverse)
  library(GenomicRanges)
  library(matrixStats)
  library(aargh)
})


# input_clone_cnv <- "/cellassign/clonealign-analysis/analysis/alleleic-imbalance/build-snv-matrix/SA501-clone-cnv.csv"
# input_alt <- "/cellassign/clonealign-analysis/analysis/alleleic-imbalance/build-snv-matrix/SA501altmat.mtx"
# input_ref <- "/cellassign/clonealign-analysis/analysis/alleleic-imbalance/build-snv-matrix/SA501refmat.mtx"
# input_cov <- "/cellassign/clonealign-analysis/analysis/alleleic-imbalance/build-snv-matrix/SA501covmat.mtx"
# input_barcode_index <- "/cellassign/clonealign-analysis/analysis/alleleic-imbalance/build-snv-matrix/SA501-barcode_index.csv"
# input_region_index <- "/cellassign/clonealign-analysis/analysis/alleleic-imbalance/build-snv-matrix/SA501-region_index.csv"
# 
# output_cnv_clone <- "/cellassign/clonealign-analysis/data/SA501/allelic-imbalance/position-by-clone-mat.Rds"
# output_cov <- "/cellassign/clonealign-analysis/data/SA501/allelic-imbalance/scrnaseq-cov-counts.Rds"
# output_ref <- "/cellassign/clonealign-analysis/data/SA501/allelic-imbalance/scrnaseq-ref-counts.Rds"

# Should be a data frame with the columns
# chr start end copy_number clone
# input_df 


snvworkflow <- function(input_clone_cnv = "x",
                         input_alt = "x",
                         input_ref = "x",
                         input_cov = "x",
                         input_barcode_index = "x",
                         input_region_index = "x",
                         output_cnv_clone = "x",
                         output_cov = "x",
                         output_ref = "x") {

  input_df <- read_csv(input_clone_cnv)
  
  input_df_has_chr <- grepl("chr", input_df$chr[1])
  
  if(!input_df_has_chr) {
    input_df$chr <- paste0("chr", input_df$chr)
  }
  
  clones <- sort(unique(input_df$clone))
  
  
  # Create a GRanges object with meta data columns corresponding to clone copy number --------
  
  df_clone_list <- split(input_df, f = input_df$clone)
  
  df_clone_list <- df_clone_list[clones]
  
  cnv_grs <- lapply(df_clone_list, makeGRangesFromDataFrame, keep.extra.columns = TRUE)  
  
  cnv_grs_no_mcols <- lapply(cnv_grs, function(x) {mcols(x) <- NULL; x})
  
  all_ranges <- GenomicRanges::disjoin(Reduce(c, cnv_grs_no_mcols))
  
  for(i in seq_len(length(cnv_grs))) {
    mcols(all_ranges)[[ clones[i] ]] <- mcols(cnv_grs[[i]])$copy_number[subjectHits(findOverlaps(all_ranges, cnv_grs[[i]]))]
  }
  
  
  
  # Create a GRanges object for the SNV data --------------------------------
  
  alt <- Matrix::readMM(input_alt)
  ref <- Matrix::readMM(input_ref)
  cov <- Matrix::readMM(input_cov)
  
  keep_sites <- Matrix::colMeans(alt) > 0.01 & Matrix::colMeans(ref) > 0.01
  
  
  alt_fltr <- t(alt[,keep_sites]) # Now site by cell
  ref_fltr <- t(ref[,keep_sites])
  cov_fltr <- t(cov[,keep_sites])
  
  # Convert to normal matrices
  ref_fltr_m <- as.matrix(ref_fltr)
  cov_fltr_m <- as.matrix(cov_fltr)
  
  
  barcode_index <- read_csv(input_barcode_index) %>% arrange(index)
  region_index <- read_csv(input_region_index) %>% arrange(index)
  
  colnames(ref_fltr_m) <- colnames(cov_fltr_m) <- barcode_index$barcode
  rownames(ref_fltr_m) <- rownames(cov_fltr_m) <- region_index$region[keep_sites]
  
  # TODO: check if chr needs appended
  ref_df <- as_data_frame(ref_fltr_m) 
  
  
  ref_df <- ref_df %>% 
    dplyr::mutate(chrom = paste0("chr", sapply(strsplit(rownames(ref_fltr_m), ":"), `[`, 1)),
                  pos = as.numeric(sapply(strsplit(rownames(ref_fltr_m), ":"), `[`, 2)),
                  start = pos - 1) %>% 
    dplyr::rename(end = pos) %>% 
    dplyr::select(chrom, start, end, everything())
  
  ref_ranges <- GenomicRanges::makeGRangesFromDataFrame(ref_df, keep.extra.columns = TRUE)
  
  cov_df <- as_data_frame(cov_fltr_m) %>% 
    dplyr::mutate(chrom = paste0("chr", sapply(strsplit(rownames(cov_fltr_m), ":"), `[`, 1)),
                  pos = as.numeric(sapply(strsplit(rownames(cov_fltr_m), ":"), `[`, 2)),
                  start = pos - 1) %>% 
    dplyr::rename(end = pos) %>% 
    dplyr::select(chrom, start, end, everything())
  
  cov_ranges <- GenomicRanges::makeGRangesFromDataFrame(cov_df, keep.extra.columns = TRUE)
  
  olaps1 <- findOverlaps(ref_ranges, all_ranges)
  olaps2 <- findOverlaps(cov_ranges, all_ranges)
  
  stopifnot(all.equal(queryHits(olaps1), queryHits(olaps2)))
  stopifnot(all.equal(subjectHits(olaps1), subjectHits(olaps2)))
  
  cnv_final <- mcols(all_ranges)[subjectHits(olaps1),]
  cov_final <- mcols(cov_ranges)[queryHits(olaps1),]
  ref_final <- mcols(ref_ranges)[queryHits(olaps1),]
  
  cnv_final_mat <- as.matrix(cnv_final)
  
  # Keep only locations with variable LOH copy number 
  usable_locations <- rowVars(cnv_final_mat) > 0 & rowMins(cnv_final_mat) == 1 & rowMaxs(cnv_final_mat) == 2
  
  cnv_final <- cnv_final[usable_locations,]
  cov_final <- cov_final[usable_locations,]
  ref_final <- ref_final[usable_locations,]
  
  # TODO: add in location info
  
  saveRDS(as.matrix(cnv_final), output_cnv_clone)
  saveRDS(as.matrix(ref_final), output_ref)
  saveRDS(as.matrix(cov_final), output_cov)
  
  message("Done")
}

aargh(snvworkflow)





