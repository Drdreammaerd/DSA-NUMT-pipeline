library(optparse)
library(Biostrings)
library(data.table)
library(stringr)
library(rtracklayer)
library(GenomicRanges)

options(scipen=999)

option_list <- list(
  make_option(c("-b", "--bed_file"), type="character", help="Path to bed file (or TSV)"), 
  make_option(c("-c", "--chain"), type="character", help="Path to chain file"),
  make_option(c("-o", "--output"), type="character", help="output file")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$bed_file)) stop("Missing -b")
if (is.null(opt$chain)) stop("Missing -c")
if (is.null(opt$output)) stop("Missing -o")

bed_file_path <- opt$bed_file
chain_path <- opt$chain
output_path <- opt$output

# Read in the TSV
bed_data_raw <- fread(bed_file_path, header = TRUE, sep = "\t")
if (nrow(bed_data_raw) == 0) {
  fwrite(bed_data_raw, output_path, sep = "\t", col.names = TRUE)
  quit(save="no")
}

# Generate an internal ID
bed_data_raw$INTERNAL_ID <- 1:nrow(bed_data_raw)

bed_data <- bed_data_raw[, c("contig", "locus_start", "locus_end", "INTERNAL_ID"), with = FALSE]
colnames(bed_data) <- c("chr", "start", "end", "ID")

# Import the chain file
chain <- rtracklayer::import(chain_path, format = "chain")

# Create a GRanges object from the BED file
grange_bed <- makeGRangesFromDataFrame(bed_data, keep.extra.columns = TRUE)

# Perform liftover
lifted <- liftOver(grange_bed, chain)

result_list <- list()
for (i in seq_along(lifted)) {
  gr <- lifted[[i]]
  if (length(gr) == 0) next
  
  original_id <- grange_bed$ID[i]
  original_size <- width(grange_bed[i])
  
  chr_groups <- split(gr, seqnames(gr))
  for (chr_name in names(chr_groups)) {
    grp <- chr_groups[[chr_name]]
    if (length(grp) == 0) next
    span_start <- min(start(grp))
    span_end <- max(end(grp))
    lifted_size <- span_end - span_start
    
    if (lifted_size < original_size * 0.5) {
      center <- as.integer((span_start + span_end) / 2)
      half_size <- as.integer(ceiling(original_size / 2))
      span_start <- max(1, center - half_size)
      span_end <- center + half_size
    }
    
    result_list[[length(result_list) + 1]] <- data.table(
      hg38_chr = as.character(chr_name),
      hg38_start = span_start,
      hg38_end = span_end,
      INTERNAL_ID = original_id
    )
  }
}

if (length(result_list) > 0) {
  lifted_bed <- rbindlist(result_list)
  # Merge with original data
  merged_data <- merge(lifted_bed, bed_data_raw, by="INTERNAL_ID", all.x=TRUE)
  merged_data[, INTERNAL_ID := NULL]
  
  # Ensure column order: hg38_chr, hg38_start, hg38_end, ... then original cols
  orig_cols <- setdiff(colnames(merged_data), c("hg38_chr", "hg38_start", "hg38_end"))
  final_cols <- c("hg38_chr", "hg38_start", "hg38_end", orig_cols)
  merged_data <- merged_data[, ..final_cols]
} else {
  # Empty output with correct header
  orig_cols <- colnames(bed_data_raw)
  orig_cols <- setdiff(orig_cols, "INTERNAL_ID")
  final_cols <- c("hg38_chr", "hg38_start", "hg38_end", orig_cols)
  merged_data <- data.table(matrix(ncol = length(final_cols), nrow = 0))
  colnames(merged_data) <- final_cols
}

# Write the output to a file
fwrite(merged_data, output_path, sep = "\t", col.names = TRUE, row.names = FALSE)
