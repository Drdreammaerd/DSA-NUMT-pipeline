library(optparse)
library(Biostrings)
library(data.table)
library(stringr)
library(rtracklayer)
library(GenomicRanges)

options(scipen=999)

option_list <- list(
  make_option(c("-b", "--bed_file"),
              type="character",
              help="Path to bed file"), 
  make_option(c("-t", "--threads"),
              type="integer",
              help="number of threads"),
  make_option(c("-c", "--chain"),
              type="character",
              help="Path to chain file"),
  make_option(c("-o", "--output"),
              type="character",
              help="output file"),
  make_option(c("-i", "--id_col"),
              type="integer", default=4,
              help="Column number for ID (default: 4)")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Validate required arguments
if (is.null(opt$bed_file)) {
  stop("Missing -b or --bed_file argument. Please provide the input bed file path.")
}
if (is.null(opt$chain)) {
  stop("Missing -c or --chain argument. Please provide the input chain file path.")
}
if (is.null(opt$output)) {
  stop("Missing -o or --output argument. Please provide the output bed file path.")
}

bed_file_path <- opt$bed_file
chain_path <- opt$chain
output_path <- opt$output
id_col <- opt$id_col

# Read in the BED file
bed_data_raw <- fread(bed_file_path, header = FALSE, sep = "\t")

# Check if id_col is valid
if (id_col > ncol(bed_data_raw)) {
  stop(paste("The ID column specified (", id_col, ") is out of bounds for the input file which has", ncol(bed_data_raw), "columns."))
}
bed_data <- bed_data_raw[, c(1:3, id_col), with = FALSE]

# Define BED file column names
bed_columns <- c("chr", "start", "end","ID")
colnames(bed_data) <- bed_columns

# Import the chain file
chain <- rtracklayer::import(chain_path, format = "chain")

# Create a GRanges object from the BED file
grange_bed <- makeGRangesFromDataFrame(bed_data, keep.extra.columns = TRUE)

# Perform liftover
lifted <- liftOver(grange_bed, chain)

# Group-aware span: for each original region, take overall span per target chromosome
# This prevents fragmentation from splitting a single NUMT into tiny pieces
result_list <- list()
for (i in seq_along(lifted)) {
  gr <- lifted[[i]]
  if (length(gr) == 0) next
  
  original_id <- grange_bed$ID[i]
  original_size <- width(grange_bed[i])
  
  # Group fragments by target chromosome
  chr_groups <- split(gr, seqnames(gr))
  for (chr_name in names(chr_groups)) {
    grp <- chr_groups[[chr_name]]
    if (length(grp) == 0) next
    span_start <- min(start(grp))
    span_end <- max(end(grp))
    lifted_size <- span_end - span_start
    
    # Expand step: if lifted region is < 50% of original size, 
    # expand symmetrically around center to match original size
    if (lifted_size < original_size * 0.5) {
      center <- as.integer((span_start + span_end) / 2)
      half_size <- as.integer(ceiling(original_size / 2))
      span_start <- max(1, center - half_size)
      span_end <- center + half_size
    }
    
    result_list[[length(result_list) + 1]] <- data.table(
      seqnames = chr_name,
      start = span_start,
      end = span_end,
      ID = original_id
    )
  }
}

if (length(result_list) > 0) {
  lifted_bed <- rbindlist(result_list)
} else {
  lifted_bed <- data.table(seqnames = character(), start = integer(), 
                           end = integer(), ID = character())
}

# Write the output to a file
fwrite(lifted_bed, output_path, sep = "\t", col.names = FALSE, row.names = FALSE)
