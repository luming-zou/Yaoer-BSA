#!/usr/bin/env Rscript

# Setting options to avoid scientific notation for better readability
options(scipen = 20)

# Reading filtered depth data
dp <- read.table("filtered_depth.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)

# Calculating SNP-Index for Sample S01
dp_mut <- ifelse(dp[, 6] > dp[, 5], dp[, 6], dp[, 5])
index_s01 <- dp_mut / (dp[, 5] + dp[, 6])

# Calculating SNP-Index for Sample S02
dp_mut <- ifelse(dp[, 6] > dp[, 5], dp[, 8], dp[, 7])
index_s02 <- dp_mut / (dp[, 7] + dp[, 8])

# Calculating Delta SNP-Index
d_index <- index_s01 - index_s02

# Creating a matrix of depth values for subsequent analysis
depth <- as.matrix(dp[, 5:8])

# Creating a data frame for SNP-Index information
index <- data.frame(
  CHROM = dp[, 1],
  POS = dp[, 2],
  INDEX_r = index_s01,
  INDEX_S = index_s02,
  D_INDEX = d_index
)

# Function to calculate delta SNP-Index in a given window
delta_index <- function(chr, start, end) {
  index_temp1 <- index_temp[index_temp[, 2] > start & index_temp[, 2] <= end, ]
  if (dim(index_temp1)[1] >= 0) {
    m <- sum(index_temp1$D_INDEX) / (end - start)  # average SNP-Index
    return(m)
  } else {
    return(0)
  }
}

# Reading in the reference genome index file
fai <- read.table("REFERENCE.fa.fai", header = FALSE, stringsAsFactors = FALSE)
index_all <- NULL

# Looping through each chromosome
for (i in 1:dim(fai)[1]) {
  D_INDEX_CHR <- NULL
  chr <- fai[i, 1]
  print(chr)
  index_temp <- index[index$CHROM == chr, ]

  # Window and step size for sliding window analysis
  window <- 1e3
  step <- 7e2

  # Producing windows
  length <- fai[i, 2]
  if (length <= window) {
    windows <- matrix(c(1, length), ncol = 2)
  } else {
    starts <- seq(1, length, step)
    ends <- starts + window
    windows <- cbind(starts, ends)
    windows <- windows[ends <= length, ]
    windows <- matrix(windows, ncol = 2)
    windows <- rbind(windows, c(max(windows[, 1]) + step, length))
  }

  # Calculating delta SNP-Index for each window
  for (j in 1:dim(windows)[1]) {
    D_INDEX_CHR <- rbind(D_INDEX_CHR, delta_index(chr, windows[j, 1], windows[j, 2]))
  }

  # Creating a data frame for delta SNP-Index information for the chromosome
  index_chr <- data.frame(
    CHROM = chr,
    STARTS = windows[, 1],
    ENDS = windows[, 2],
    D_INDEX = D_INDEX_CHR[, 1]
  )

  # Combining data for all chromosomes
  index_all <- rbind(index_all, index_chr)
}

# Outputting delta SNP-Index to a CSV file
write.table(index_all, "delta_snp_index.csv", row.names = FALSE, sep = "\t", quote = FALSE)
