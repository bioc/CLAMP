
library(dplyr)


#' Row-wise scaling (mean 0, sd 1)
#'
#' Standardize each row of a numeric matrix to have mean 0 and standard deviation 1.
#'
#' @param x A numeric matrix.
#' @return A matrix of the same shape, with each row scaled independently.
#' @export
tscale <- function(x) {
  row_means <- rowMeans(x)
  row_sds <- sqrt(rowMeans((x - row_means)^2))
  row_sds[row_sds == 0] <- 1  # avoid division by zero
  sweep(sweep(x, 1, row_means), 1, row_sds, "/")
}


#' Compare latent variable loadings against a target using correlation or other statistics
#'
#' This function compares two sets of latent variable loadings (`res1`, `res2`) with respect to a binary or continuous
#' target matrix. It computes the maximal association (e.g., correlation, t-statistic, or AUC) between each latent variable
#' and each column in the target, and returns a paired comparison plot and summary.
#' @export
compareBs<-function(res1, res2, target, method = "p", xlab="1", ylab="2", stat.method="t") {
  extract_B <- function(res) {
    if (class(res)[1] == "list") {
      return(as.matrix(res$B))
    } else if (class(res)[1] == "rsvd") {
      return(t(res$v))
    } else {
      return(as.matrix(res))
    }
  }

  B1 <- extract_B(res1)
  B2 <- extract_B(res2)
  noNA=!apply(target,1, function(x){any(is.na(x))})
  B1=B1[,noNA]
  B2=B2[,noNA]
  target=target[noNA,]

  if (method %in% c("s", "p")) {
    mat1 <- cor(t(B1), target, method = method)
    mat2 <- cor(t(B2), target, method = method)
    cor1 <- apply(mat1, 2, max, na.rm=T)
    cor2 <- apply(mat2, 2, max, na.rm=T)
    idx1 <- apply(mat1, 2, which.max)
    idx2 <- apply(mat2, 2, which.max)
  } else if (method == "a") {
    mat1 <- allAgainstAllAUCs(t(B1), target)
    mat2 <- allAgainstAllAUCs(t(B2), target)
    cor1 <- apply(mat1, 2, max)
    cor2 <- apply(mat2, 2, max)
    idx1 <- apply(mat1, 2, which.max)
    idx2 <- apply(mat2, 2, which.max)
  } else if (method == "t") {
    mat1 <- allAgainstAllTstats(t(B1), target)
    mat2 <- allAgainstAllTstats(t(B2), target)
    cor1 <- apply(mat1, 1, max)
    cor2 <- apply(mat2, 1, max)
    idx1 <- apply(mat1, 1, which.max)
    idx2 <- apply(mat2, 1, which.max)
  }
  Labels <- colnames(target)

  df <- data.frame(Cor1 = cor1, Cor2 = cor2, BestIdxRes1=idx1, BestIdxRes2=idx2,  Label = Labels)
  df$corMean=(df$Cor1+df$Cor2)/2
  if(stat.method=="t"){
    pval <- t.test(df$Cor2, df$Cor1, paired = TRUE, alternative = "greater")$p.value
  }else{
    pval <- wilcox.test(df$Cor2, df$Cor1, paired = TRUE, alternative = "greater")$p.value
  }
  method_label <- switch(method,
                         p = "Pearson correlation",
                         s = "Spearman correlation",
                         a = "AUC",
                         t = "T-statistic"
  )
  pl<-ggplot(df, aes(x = Cor1, y = Cor2, label = Label)) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red", linewidth = 1) +
    geom_text_repel() +
    labs(x = paste("Max", method_label, xlab),
         y = paste("Max", method_label, ylab))+
    annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.1, label = paste("p =", signif(pval, 3)))+
    #    annotate("text", x = Inf, y = -Inf, label = paste("p =", signif(pval, 3)),
    #            hjust = 1.1, vjust = -0.5, size = 4)+
    theme_minimal()
  return(list(plot=pl, df=df))

}


#' Print a concatenated message
#'
#' Wrapper around \code{message()} that pastes arguments together into a single string.
#'
#' @param ... Character strings to concatenate and print.
#'
#' @export
mymessage <- function(...) {
  message(paste(...))
}

#' Get maximum AUC per latent variable
#'
#' Summarizes a cross-validation results data frame to extract the highest AUC value
#' associated with each latent variable (LV).
#'
#' @param summary A data frame (e.g., from \code{crossVal()}) .
#' @param verbose Logical; if \code{TRUE}, prints counts of LVs exceeding AUC thresholds. Default is \code{FALSE}.
#'
#' @return A data frame with columns \code{LV index} and \code{max_AUC}.
#'
#' @export
getMaxAUC=function(summary, verbose=F){

  max_auc_per_lv <- summary %>%
    group_by(.data$`LV index`) %>%
    summarize(max_AUC = max(.data$AUC, na.rm = TRUE)) %>%
    ungroup()

  if(verbose){
    message(paste("There are", sum(max_auc_per_lv$max_AUC>0.70), " LVs with AUC>0.70"))
    message(paste("There are", sum(max_auc_per_lv$max_AUC>0.90), " LVs with AUC>0.90"))
  }
  max_auc_per_lv
}


#' Count number of latent variables exceeding AUC thresholds
#'
#' Given a summary data frame from cross-validation, reports the number of latent variables
#' with maximum AUC exceeding 0.7, 0.8, and 0.9.
#'
#' @param summary A data frame with \code{LV index} and \code{AUC} columns.
#'
#' @return A named numeric vector with counts for thresholds 0.7, 0.8, and 0.9.
#'
#' @export
getAUCstats=function(summary){
  out=getMaxAUC(summary)
  unlist(lapply(c(0.7,0.8, 0.9), function(x){sum(out$max_AUC>x)}))
}


#' Greedy maximum correspondence from correlation matrix
#'
#' Finds a one-to-one assignment (permutation) between rows and columns of a square correlation matrix
#' that maximizes the total correlation score, using a greedy algorithm.
#'
max_correspondence_greedy <- function(cor_mat) {
  n <- nrow(cor_mat)
  used_rows <- rep(FALSE, n)
  used_cols <- rep(FALSE, n)
  assignment <- integer(n)

  corr_entries <- as.data.frame(which(!is.na(cor_mat), arr.ind = TRUE))
  corr_entries$val <- cor_mat[cbind(corr_entries$row, corr_entries$col)]
  corr_entries <- corr_entries[order(-corr_entries$val), ]

  for (i in seq_len(nrow(corr_entries))) {
    r <- corr_entries$row[i]
    c <- corr_entries$col[i]
    if (!used_rows[r] && !used_cols[c]) {
      assignment[r] <- c
      used_rows[r] <- TRUE
      used_cols[c] <- TRUE
    }
  }

  list(permutation = assignment, sum = sum(cor_mat[cbind(1:n, assignment)]))
}

#' Download and read a GMT file from a URL
#'
#' Downloads a Gene Matrix Transposed (GMT) file from a specified URL, reads it into R as a list,
#' and removes the temporary file afterward.
#'
#' @param url A character string specifying the URL to a GMT file.
#'
#' @return A named list where each element is a character vector of gene names for a given gene set.
#' @export
#'
getGMT <- function(url, name = NULL, cache_dir = NULL, redownload = FALSE) {
  if (is.null(name)) {
    name <- sub(".*[=]", "", url)
    message("Auto-detected name: ", name)
  }
  if (is.null(cache_dir)) {
    cache_dir <- system.file("extdata", package = "PLIER2")
  }
  cache_file <- file.path(cache_dir, paste0(name, ".gmt"))
  if (!file.exists(cache_file) || redownload) {
    message("Downloading ", name, " from Enrichr...")
    if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)
    download.file(url, cache_file)
  } else {
    message("Using cached file for ", name)
  }
  read_gmt(cache_file)
}





#' Read a GMT file into a list
#'
#' Parses a local GMT file and returns a list of gene sets, with each gene set represented as
#' a character vector of unique gene names.
#'
#' @param filename A character string giving the path to a .gmt file.
#'
#' @return A named list where each element is a character vector of gene names.
#' @export
read_gmt=function (filename) {
  gmt = list()
  lines = readLines(filename)
  for (line in lines) {
    line = gsub("\"", "", trimws(line))
    sp = unlist(strsplit(line, "\t"))
    sp[3:length(sp)] = gsub(",.*$", "", sp[3:length(sp)])
    gmt[[sp[1]]] = sort(unique(sp[3:length(sp)]))
  }
  return(gmt)
}

#' Convert a list of GMT gene sets to a sparse matrix
#'
#' Converts a list of named gene sets (e.g., from \code{getGMT()}) into a sparse binary matrix
#' where rows are genes, columns are gene sets, and entries are 1 if the gene is in the set.
#'
#' @param gmtList A nested list of gene sets. Outer names are gene set names; each entry is a character vector of gene names.
#'
#' @return A sparse binary matrix with genes as rows and gene sets as columns.
#' @export
gmtListToSparseMat=function(gmtList){

  allnames=unlist(lapply(gmtList, names))
  #there are usually no duplicates
  stopifnot(all(table(allnames)==1))
  allGenes=unique(unlist(lapply(gmtList, unlist)))

  row_indices <- integer(0)
  col_indices <- integer(0)
  values <- integer(0)
  for (gmt in seq_along(gmtList)) {
    for (path in names(gmtList[[gmt]])) {
      thisPathGenes <- gmtList[[gmt]][[path]]
      iiGenes <- match(thisPathGenes, allGenes)
      iPath <- match(path, allnames)

      # Store indices and values
      row_indices <- c(row_indices, iiGenes)
      col_indices <- c(col_indices, rep(iPath, length(iiGenes)))
      values <- c(values, rep(1, length(iiGenes)))
    }
  }

  # Use sparseMatrix to create the matrix in one go

  pathMat=sparseMatrix(i = row_indices, j = col_indices, x = values,
                       dims = c(length(allGenes), length(allnames)))
  rownames(pathMat) = allGenes
  colnames(pathMat) = allnames
  pathMat
}

#' Find common row names between two matrices or data frames
#'
#' Returns the intersection of row names shared by two input objects.
#'
#' @param data1 A matrix, data frame, or similar object with row names.
#' @param data2 A matrix, data frame, or similar object with row names.
#'
#' @return A character vector of row names common to both inputs.
#'
#' @export
commonRows=function(data1, data2){
  intersect(rownames(data1), rownames(data2))
}


#' Clean a Filebacked Big Matrix (FBM) by log-transforming and handling NAs
#'
#' This function inspects an FBM to determine if log-transformation is needed (based on value range)
#' and whether NA values are present. If the maximum value is ≥ 100, it applies a log2(x + 1)
#' transformation in-place. If any NA values are detected, they are replaced with 0.
#'
#' @param fbm A \code{bigmemory::FBM} or \code{bigstatsr::FBM} object.
#'
#' @return A list with:
#' \describe{
#'   \item{\code{max_value}}{The maximum value encountered in the FBM (after log transformation if applied).}
#'   \item{\code{had_na}}{Logical indicating whether any NA values were found and filled.}
#' }
#'
#' @details
#' Modifies the FBM in place. Uses \code{bigstatsr::big_apply()} to process in parallel-safe chunks.
#'
#' @export
cleanFBM=function(fbm){
  # Step 1: Check for NA and max value
  max_value <- -Inf
  has_na <- FALSE

  big_apply(fbm, a.FUN = function(X, ind) {
    max_value <<- max(max_value, max(X[ind, ], na.rm = TRUE))
    if (anyNA(X[ind, ])) {
      has_na <<- TRUE
    }
    NULL  # No return, just updating global values
  }, ind = rows_along(fbm))

  # Step 2: Log2 transform if necessary
  if (max_value >= 100) {
    message("Applying log2 transformation")
    big_apply(fbm, a.FUN = function(X, ind) {
      X[ind, ] <- log2(X[ind, ] + 1)
    }, ind = rows_along(fbm))
  }
  else{
    message("Already on log scale")
  }

  # Step 3: Fill NAs with 0 if necessary
  if (has_na) {
    message("Filling NAs with 0")
    big_apply(fbm, a.FUN = function(X, ind) {
      X[ind, ][is.na(X[ind, ])] <- 0
      NULL
    }, ind = rows_along(fbm))
  }
  else{
    message("No NA values found")
  }

  return(list(max_value = max_value, had_na = has_na))
}


#' Compute row-wise sum and sum of squares for a Filebacked Big Matrix
#'
#' Efficiently computes row sums and row sum-of-squares for a \code{bigstatsr::FBM} using
#' column-wise chunking, suitable for large datasets that cannot be loaded fully into memory.
#'
#' @param fbm A \code{bigstatsr::FBM} object.
#' @param chunk_size Number of columns to process at a time. Default is 1000.
#'
#' @return A list with two numeric vectors:
#' \describe{
#'   \item{\code{row_sums}}{Sum of each row.}
#'   \item{\code{row_sums_sq}}{Sum of squares of each row.}
#' }
#'
#' @export
computeRowStatsFBM <- function(fbm, chunk_size = 1000) {
  n_rows <- nrow(fbm)
  n_cols <- ncol(fbm)

  # Initialize vectors to store row sums and row sums of squares
  row_sums <- numeric(n_rows)
  row_sums_sq <- numeric(n_rows)

  # Iterate over columns in chunks
  for (start_col in seq(1, n_cols, by = chunk_size)) {
    end_col <- min(start_col + chunk_size - 1, n_cols)

    # Extract the current chunk of columns
    col_chunk <- fbm[, start_col:end_col]

    # Update row sums and sums of squares
    row_sums <- row_sums + rowSums(col_chunk)
    row_sums_sq <- row_sums_sq + rowSums(col_chunk^2)
  }

  # Compute row means and variances
  row_means <- row_sums / n_cols
  row_variances <- (row_sums_sq / n_cols) - (row_means^2)

  return(list(row_means = row_means, row_variances = row_variances))
}

#' Filter rows of a Filebacked Big Matrix based on mean and variance
#'
#' Filters an FBM based on row-level mean and variance thresholds, returning a new FBM
#' with only the selected rows.
#'
#' @param fbm A \code{bigstatsr::FBM} object.
#' @param rowStats A list with numeric vectors \code{row_means} and \code{row_variances}.
#' @param mean_cutoff Optional minimum mean threshold; rows with means below this are removed.
#' @param var_cutoff Optional minimum variance threshold; rows with variances below this are removed.
#' @param backingfile A character string specifying the filename (without extension) for the new FBM. Default is \code{"filtered_fbm"}.
#'
#' @return A list with:
#' \describe{
#'   \item{\code{fbm_filtered}}{A new FBM object containing only filtered rows.}
#'   \item{\code{kept_rows}}{Indices of rows retained in the filtering step.}
#' }
#' @details
#' This function creates a new FBM and copies over only the rows that pass the filtering criteria.
#' The original FBM is unchanged.
#' @export
filterFBM<- function(fbm, rowStats, mean_cutoff = NULL, var_cutoff = NULL, backingfile = "filtered_fbm") {
  row_means <- rowStats$row_means
  row_variances <- rowStats$row_variances

  # Determine rows to keep based on cutoffs
  keep_rows <- rep(TRUE, length(row_means))  # Default: keep all rows

  if (!is.null(mean_cutoff)) {
    keep_rows <- keep_rows & (row_means >= mean_cutoff)
  }

  if (!is.null(var_cutoff)) {
    keep_rows <- keep_rows & (row_variances >= var_cutoff)
  }

  # Number of rows to keep
  n_kept <- sum(keep_rows)

  if (n_kept == 0) {
    stop("No rows meet the filtering criteria.")
  }

  # Create a new FBM with the filtered data
  fbm_filtered <- FBM(n_kept, ncol(fbm), backingfile = backingfile)
  fbm_filtered[] <- fbm[keep_rows, ]

  return(list(fbm_filtered = fbm_filtered, kept_rows = which(keep_rows)))
}

#' Z-score a filtered expression matrix for PLIER2
#'
#' Centers each gene to mean 0 and scales to unit variance.
#'
#' @param Y_filtered Numeric matrix (genes × samples) returned by preprocessPLIER2
#' @param rowStats   Data frame with numeric columns `mean` and `variance`,
#'                   row-named to match `rownames(Y_filtered)`
#'
#' @return Numeric matrix of the same dimensions as `Y_filtered`, with each
#'         row centered and scaled.
#' @export
zscorePLIER2 <- function(Y_filtered, rowStats) {
  # 1) Input validation
  if (!is.matrix(Y_filtered) || !is.numeric(Y_filtered)) {
    stop("`Y_filtered` must be a numeric matrix (genes × samples).")
  }
  if (!is.data.frame(rowStats) ||
      !all(c("mean", "variance") %in% colnames(rowStats))) {
    stop("`rowStats` must be a data.frame with columns 'mean' and 'variance'.")
  }
  # 2) Align rowStats to Y_filtered
  if (!all(rownames(Y_filtered) %in% rownames(rowStats))) {
    stop("Row names of `Y_filtered` and `rowStats` do not match.")
  }
  rowStats <- rowStats[rownames(Y_filtered), , drop = FALSE]
  # 3) Ensure numeric
  mu  <- as.numeric(rowStats$mean)
  var <- as.numeric(rowStats$variance)
  if (any(is.na(mu)) || any(is.na(var))) {
    stop("Missing values detected in 'mean' or 'variance'.")
  }
  if (any(var <= 0)) {
    stop("All variances must be positive; zero or negative found.")
  }
  # 4) Compute standard deviation
  sd  <- sqrt(var)
  # 5) Center and scale
  #    subtract mu from each row, then divide by sd
  Y_centered <- sweep(Y_filtered, 1L, mu,  "-")
  Y_scaled   <- sweep(Y_centered, 1L, sd,  "/")
  # 6) Return
  return(Y_scaled)
}


#' Preprocess a bigstatsr FBM for PLIER2
#'
#' Makes a writable copy of the input FBM, cleans it (log2 transform if needed, fill NAs),
#' filters rows by mean/variance, and returns the filtered FBM plus stats and indices.
#'
#' @param fbm A bigstatsr::FBM (genes × samples), possibly read-only.
#' @param mean_cutoff Numeric or NULL. Minimum row mean to keep (NULL = no mean filter).
#' @param var_cutoff  Numeric or NULL. Minimum row variance to keep (NULL = no var filter).
#' @param backingfile Character or NULL. Base name for the *copy* FBM and filtered FBM on disk.
#'                    If NULL, defaults to paste0(fbm$backingfile, "_preproc") and "_filtered".
#' @return A list with:
#'   \item{fbm_filtered}{The filtered FBM (writable).}
#'   \item{rowStats}{List with row_means & row_variances for fbm_filtered.}
#'   \item{kept_rows}{Integer vector of original row indices that were retained.}
#' @export
preprocessPLIER2FBM <- function(fbm,
                                mean_cutoff = NULL,
                                var_cutoff  = NULL,
                                backingfile = NULL) {
  # 1. Choose base names
  base_bk <- if (is.null(backingfile)) paste0(fbm$backingfile, "_preproc") else backingfile

  # 2. Make a writable copy
  fbm_copy <- FBM(
    nrow        = nrow(fbm),
    ncol        = ncol(fbm),
    backingfile = base_bk,
    create_bk   = TRUE
  )
  # copy all data
  fbm_copy[] <- fbm[]

  # 3. Clean in-place (log2 if needed, fill NAs)
  cleanFBM(fbm_copy)

  # 4. Compute row stats on cleaned copy
  rs_all <- computeRowStatsFBM(fbm_copy)

  # 5. Filter rows, writing to a new filtered FBM
  filt_bk <- paste0(base_bk, "_filtered")
  filter_res <- filterFBM(
    fbm_copy,
    rowStats    = rs_all,
    mean_cutoff = mean_cutoff,
    var_cutoff  = var_cutoff,
    backingfile = filt_bk
  )
  fbm_filtered <- filter_res$fbm_filtered
  kept_rows    <- filter_res$kept_rows

  # 6. Subset stats to kept rows
  stats_filt <- list(
    row_means     = rs_all$row_means[kept_rows],
    row_variances = rs_all$row_variances[kept_rows]
  )

  list(
    fbm_filtered = fbm_filtered,
    rowStats     = stats_filt,
    kept_rows    = kept_rows
  )
}

#' Z-score a filtered FBM in-place
#'
#' Standardizes each row of an FBM using provided row means and variances.
#'
#' @param fbm_filtered A bigstatsr::FBM produced by preprocessPLIER2FBM().
#' @param rowStats A list with row_means and row_variances from that FBM.
#' @param chunk_size Columns per block (default 1000).
#' @export
zscorePLIER2FBM <- function(fbm_filtered, rowStats, chunk_size = 1000) {
  message("Applying Z-score transformation")
  means <- rowStats$row_means
  sds   <- sqrt(rowStats$row_variances)
  sds[sds == 0] <- 1

  for (start in seq(1, ncol(fbm_filtered), by = chunk_size)) {
    end <- min(start + chunk_size - 1, ncol(fbm_filtered))
    mat <- fbm_filtered[, start:end]
    mat <- sweep(mat, 1, means, FUN = "-")
    mat <- sweep(mat, 1, sds,   FUN = "/")
    fbm_filtered[, start:end] <- mat
  }

  invisible(NULL)
}

#' Preprocess an expression matrix for PLIER2
#'
#' Filters genes by mean expression and variance, returning the filtered matrix
#' and per-gene statistics.
#'
#' @param Y Numeric matrix of gene expression (rows = genes, cols = samples)
#' @param mean_cutoff Numeric. Minimum row-mean required to keep a gene (default 0).
#' @param var_cutoff  Numeric. Minimum row-variance required to keep a gene (default 0).
#'
#' @return A list with components:
#'   - **fbm_filtered**: filtered matrix (genes × samples)
#'   - **rowStats**: data.frame with columns `mean` and `variance` for each kept gene
#'   - **kept_rows**: integer vector of the original row indices that were kept
#'
#' @export
preprocessPLIER2 <- function(Y, mean_cutoff = 0, var_cutoff = 0) {
  if (!is.matrix(Y) || !is.numeric(Y)) {
    stop("`Y` must be a numeric matrix (genes × samples).")
  }
  # Compute per‐gene statistics
  row_mean <- rowMeans(Y, na.rm = TRUE)
  row_var  <- apply(Y, 1, stats::var,  na.rm = TRUE)

  rowStats <- data.frame(
    mean     = row_mean,
    variance = row_var,
    stringsAsFactors = FALSE
  )
  rownames(rowStats) <- rownames(Y)

  # Identify genes passing both thresholds
  keep <- which(rowStats$mean >= mean_cutoff & rowStats$variance >= var_cutoff)
  if (length(keep) == 0) {
    stop("No genes passed the mean/variance filters.")
  }

  # Subset matrix and stats
  Y_filtered       <- Y[keep, , drop = FALSE]
  rowStats_filtered <- rowStats[keep, , drop = FALSE]

  return(list(
    Y_filtered = Y_filtered,
    rowStats   = rowStats_filtered,
    kept_rows  = keep
  ))
}


#' visualize the top genes contributing to the LVs
#'
#' @param plierRes the result returned by PLIER
#' @param data the data to be displayed in a heatmap, typically the z-scored input data (or some subset thereof)
#' @param priorMat the same gene by geneset binary matrix that was used to run PLIER
#' @param top the top number of genes to use
#' @param index the subset of LVs to display
#' @param allLVs plot even the LVs that have no pathway association
#' @param ... Additional arguments to be passed to pheatmap, such as a column annotation data.frame (annotation_col). See ?pheatmap for details.
#' @export
plotTopZ=function(plierRes, data, priorMat, top=10, index=NULL, allLVs=F,...){
  data=data[rownames(plierRes$Z),]
  priorMat=priorMat[rownames(plierRes$Z),]
  ii=which(colSums(plierRes$U)>0)
  if(!allLVs){
    if(! is.null(index)){
      ii=intersect(ii,index)
    }
  }
  else{
    ii=index
  }

  tmp=apply(-plierRes$Z[, ii, drop=F],2,rank)
  nn=character()
  nncol=character()
  nnpath=character()
  nnindex=double()
  for (i in 1:length(ii)){
    nn=c(nn,nntmp<-names(which(tmp[,i]<=top)))
    nncol=c(nncol, rep(rownames(plierRes$B)[ii[i]], length(nntmp)))
    nnpath=c(nnpath,rowSums(priorMat[nntmp,plierRes$U[,ii[i]]>0, drop=F])>0)
    nnindex=c(nnindex,rep(ii[i], length(nntmp)))

  }
  names(nncol)=nn
  nncol=strtrim(nncol, 30)

  nnrep=names(which(table(nn)>1))
  if(length(nnrep)>0){
    nnrep.im=match(nnrep,nn)
    nn=nn[-nnrep.im]
    nncol=nncol[-nnrep.im]
    nnpath=nnpath[-nnrep.im]
    nnindex=c(nnindex,rep(ii[i], length(nntmp)))

  }
  nnpath[nnpath=="TRUE"]="inPathway"
  nnpath[nnpath=="FALSE"]="notInPathway"

  nncol=as.data.frame(list(nncol,nnpath))

  names(nncol)=c("pathway", "present")
  ll=c(inPathway="black", notInPathway="beige")

  anncol=list(present=ll)
  toPlot=tscale(data[nn,])



  maxval=max(abs(toPlot))

  pheatmap(toPlot, breaks=seq(-maxval, maxval, length.out = 99),color=colorpanel(100, "green", "white", "red"),annotation_row=nncol, show_colnames = F, annotation_colors = anncol, ...)
}





library(ComplexHeatmap)
library(circlize)

plotTopZ_Complex = function(plierRes, data, priorMat, top = 10, top.pathway = 5,
                            index = NULL, allLVs = FALSE, Zheat=FALSE) {
  data = data[rownames(plierRes$Z), ]
  priorMat = priorMat[rownames(plierRes$Z), ]

  ii = which(colSums(plierRes$U) > 0)
  if (!allLVs) {
    if (!is.null(index)) ii = intersect(ii, index)
  } else {
    ii = index
  }

  tmp = apply(-plierRes$Z[, ii, drop = FALSE], 2, rank)
  nn = unique(unlist(apply(tmp, 2, function(x) names(which(x <= top)))))
  nn=unique(sort(nn))

  data_sub = t(scale(t(data[nn, , drop = FALSE])))

  # build annotation for genes: inPathway or not
  nnpath = sapply(seq_along(ii), function(i) {
    gene_idx = match(nn, rownames(priorMat))
    col_idx = which(plierRes$U[, ii[i]] > 0)
    Matrix::rowSums(priorMat[gene_idx, col_idx, drop = FALSE]) > 0


  })
  nnpath = rowSums(nnpath) > 0
  gene_annot = rowAnnotation(present = nnpath,
                             col = list(present = c("TRUE" = "black", "FALSE" = "beige")))

  # build gene × top pathway binary matrix
  top_pathways = unique(unlist(lapply(ii, function(i) {
    names(sort(plierRes$U[, i], decreasing = TRUE))[1:top.pathway]
  })))
  gene_idx=match(nn, rownames(priorMat))
  path_idx=match(top_pathways, colnames(priorMat))
  pathway_mat = priorMat[gene_idx, path_idx, drop = FALSE]
  pathway_mat = as.matrix(pathway_mat > 0)
  pathway_mat=pathway_mat[, colSums(pathway_mat)>0]

  ht1 = Heatmap(data_sub,
                name = "expression",
                show_row_names = TRUE,
                show_column_names = FALSE,
                cluster_rows = TRUE,

                cluster_columns = TRUE,
                width = unit(7, "cm"),
                row_dend_width = unit(0, "mm"))

  col_fun = circlize::colorRamp2(c(0, 1), c("white", "black"))

  ht2 = Heatmap(pathway_mat+1-1,
                name = "in pathway",
                #      col = c("TRUE" = "red", "FALSE" = "white"),
                show_row_names = TRUE,
                col = col_fun,
                show_column_names = TRUE,
                cluster_rows = FALSE,
                cluster_columns = FALSE,
                column_names_rot = 45,
                row_names_gp = gpar(fontsize = 8),
                column_names_gp = gpar(fontsize = 8),
                row_dend_width = unit(0, "mm"),
                width = unit(6, "cm")
  )
  if(Zheat){
    gene_idx=match(nn, rownames(plierRes$Z))
    z_sub = scale(plierRes$Z[gene_idx, index, drop = FALSE], center = F)

    z_col_fun = circlize::colorRamp2(c(0, max(z_sub, na.rm = TRUE)),
                                     c("#e5f5e0", "#31a354"))

    ht_z = Heatmap(z_sub,
                   name = "Z",
                   col = z_col_fun,
                   cluster_rows = TRUE,
                   cluster_columns = TRUE,
                   show_row_names = TRUE,
                   show_column_names = TRUE,
                   column_names_rot = 90,
                   column_names_gp = gpar(fontsize = 10),
                   width = unit(1, "cm"))

  }
  lv_labels = colnames(plierRes$Z)[max.col(scale(plierRes$Z,center = F), ties.method = "first")]
  names(lv_labels) = rownames(plierRes$Z)
  gene_lv = lv_labels[rownames(data_sub)]

  set.seed(1)
  row_annot = rowAnnotation(LV = gene_lv,
                            #    col = list(LV = lv_colors),
                            show_annotation_name = FALSE)
  if(!Zheat){
    draw(row_annot+ht1 + ht2 , row_dend_side = "left")
  }
  else{
    draw(row_annot+ht_z+ht1 + ht2 , row_dend_side = "left")
  }
}

library(Matrix)
library(matrixStats)

corQuantile <- function(expr, pathways, quant = 0.75) {
  expr <- as.matrix(expr)
  pathways <- as(pathways, "dgCMatrix")

  result <- numeric(ncol(pathways))
  names(result) <- colnames(pathways)

  for (j in seq_len(ncol(pathways))) {
    genes_in_path <- which(pathways[, j] != 0)
    if (length(genes_in_path) >= 2) {
      sub_expr <- expr[genes_in_path, , drop = FALSE]
      cmat <- cor(t(sub_expr))
      result[j] <- quantile(cmat[lower.tri(cmat)], quant)
    } else {
      result[j] <- NA
    }
  }

  return(result)
}


getCorrelationMat<-function(zscore_data){
  cor_mat=tcrossprod(zscore_data)/(ncol(zscore_data)+1)
}



upper_outlier_threshold <- function(x, coef = 2) {
  qs <- quantile(x, probs = c(0.25, 0.75), names = FALSE, type = 7, na.rm = TRUE)
  iqr <- qs[2] - qs[1]
 qs[2] + coef * iqr

}



pathway_summary_with_eigengene <- function(expr, pathways, quant = 0.75, gene_cor_thresh = 0.5) {
  expr <- as.matrix(expr)

  pathways <- as(pathways, "dgCMatrix")

  quant_vec <- numeric(ncol(pathways))

  strong_gene_list <- vector("list", ncol(pathways))
  names(quant_vec)  <- names(strong_gene_list) <- colnames(pathways)

  for (j in seq_len(ncol(pathways))) {
    genes <- which(pathways[, j] != 0)
    if (length(genes) >= 2) {
      # quantile of gene-gene correlation
      cvals <- cor_mat[genes, genes]
      quant_vec[j] <- quantile(cvals[lower.tri(cvals)], quant)

      # compute eigengene (1st PC of gene expression)
      e <- expr[genes, , drop = FALSE]
      svd_res <- svd(scale(e, center = TRUE, scale = FALSE), nu = 0, nv = 1)
      eigengene <- svd_res$v[, 1]

      # enforce positive orientation
      avg_cor <- mean(cor(eigengene, t(e)))
      if (avg_cor < 0) eigengene <- -eigengene



      # identify genes with correlation > threshold
      gene_corrs <- cor(eigengene, t(e))
      strong_genes <- rownames(e)[which(gene_corrs > gene_cor_thresh)]
      strong_gene_list[[j]] <- strong_genes

    } else {
      quant_vec[j] <- NA
      strong_gene_list[[j]] <- character(0)
    }
  }

  data.frame(
    pathway = colnames(pathways),
    cor_quantile = quant_vec,
    strong_genes = I(strong_gene_list),
    row.names = NULL
  )
}

combinePlierBaseResults=function(plierBase1, plierBase2){
  ncol1=ncol(plierBase1$Z)
  ncol2=ncol(plierBase2$Z)
  if(!all (rownames(plierBase2$Z) %in% rownames(plierBase1$Z))){
    message("genes in plierBase2 must be a subset of genes in plierBase1")
  }
  padZ=matrix(0,nrow=nrow(plierBase1$Z), ncol=ncol2)
  plierBase1$Z=cbind(plierBase1$Z, padZ)
  plierBase1$Z[rownames(plierBase2$Z), (ncol1+1):(ncol1+ncol2)]=plierBase2$Z
  plierBase1$B=rbind(plierBase1$B, plierBase2$B)
  colnames(plierBase1$Z)<-rownames(plierBase1$B)<-paste("LV", 1:ncol(plierBase1$Z))
plierBase1
}


findSplineMax <- function(x, y, n = 1000, spar = NULL) {
  fit <- smooth.spline(x, y, spar = spar)
  grid <- seq(min(x), max(x), length.out = n)
  pred <- predict(fit, grid)
  max_idx <- which.max(pred$y)
  list(x = pred$x[max_idx], y = pred$y[max_idx])
}

squashZscore=function(zdata, maxScore=2){
  2*tanh(zdata/2)
}


allAgainstAllTstats <- function(B, target) {
  B <- as.matrix(B)  # samples × features
  target <- as.matrix(target)  # samples × targets

  n_pos <- colSums(target == 1)
  n_neg <- colSums(target == 0)

  mu1 <- crossprod(target == 1, B) / n_pos     # targets × features
  mu0 <- crossprod(target == 0, B) / n_neg

  var1 <- crossprod(target == 1, B^2) / n_pos - mu1^2
  var0 <- crossprod(target == 0, B^2) / n_neg - mu0^2

  se <- sqrt(var1 / n_pos + var0 / n_neg)
  t_stats <- (mu1 - mu0) / se

  return(t_stats)  # targets × features
}

allAgainstAllAUCs=function(B, target){
  B=as.matrix(B)
  ranks <- colRanks(B, ties.method = "average")  # rows = samples, cols = features

  n_pos <- colSums(target == 1)
  n_neg <- colSums(target == 0)
  show(dim(target))
  show(dim(ranks))
  pos_mean_rank <-  ranks %*% (target == 1)
  pos_mean_rank=sweep(pos_mean_rank,2,n_pos*(n_pos+1)/2, "-")
  #  / n_pos  # targets × features
  auc_matrix <- sweep(pos_mean_rank, 2, (n_pos) *n_neg, "/")  # targets × features
}

#' Compute counts-per-million (CPM) for PLIER2 pipelines
#'
#' @param counts A numeric matrix or data.frame of raw counts (genes × samples).
#' @return A numeric matrix of CPM values (same dimensions), ready for PLIER2 input.
#' @examples
#' mat <- matrix(1:12, nrow = 3)
#' cpmPLIER2(mat)
#' @export
cpmPLIER2 <- function(counts) {
  mat <- if (is.data.frame(counts)) as.matrix(counts) else counts
  stopifnot(is.numeric(mat), length(dim(mat)) == 2)
  lib_sizes <- colSums(mat, na.rm = TRUE)
  if (any(lib_sizes == 0)) {
    warning("Some samples have zero total counts – CPM will be Inf/NaN.")
  }
  sweep(mat, 2, lib_sizes, "/") * 1e6
}


#' Compute CPM on a file-backed matrix for PLIER2 (in-place)
#'
#' @param fbm_counts A bigstatsr::FBM of raw counts (genes × samples).
#' @param block_size Numeric; rows per block (default 1000).
#' @return Invisibly returns the modified FBM (now holding CPM values).
#' @examples
#' # fbm <- FBM(1000, 10, backingfile="raw", backingpath=".")
#' # … fill fbm with counts …
#' cpmPLIER2FBM(fbm)
#' @export
cpmPLIER2FBM <- function(fbm_counts, block_size = 1000) {
  if (!inherits(fbm_counts, "FBM")) {
    stop("`fbm_counts` must be a bigstatsr::FBM object.")
  }
  n_r <- nrow(fbm_counts)
  n_c <- ncol(fbm_counts)
  lib_sizes <- numeric(n_c)

  # 1) compute library sizes
  for (rs in seq(1, n_r, by = block_size)) {
    rows <- rs:min(rs + block_size - 1, n_r)
    lib_sizes <- lib_sizes + colSums(fbm_counts[rows, , drop = FALSE])
  }

  # 2) in-place CPM
  for (rs in seq(1, n_r, by = block_size)) {
    rows <- rs:min(rs + block_size - 1, n_r)
    block <- fbm_counts[rows, , drop = FALSE]
    block <- sweep(block, 2, lib_sizes, "/") * 1e6
    fbm_counts[rows, ] <- block
  }

  invisible(fbm_counts)

}
