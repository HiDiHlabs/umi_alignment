# library(optparse)

# # Define the option parser
# option_list <- list(
#     make_option(
#         c("-i", "--input"), type = "character", default = NULL,
#         help = "Input file path", metavar = "FILE"),
#     make_option(
#         c("-o", "--output"), type = "character", default = NULL,
#         help = "Output png"),
#     # make_option(
#     #     c("-r", "--rds"), type = "character", default = NULL,
#     #     help = "Output rds"),
#     make_option(
#         c("-b", "--binSize"), type = "integer", default = 50,
#         help = "BinSize")
# )

# # Parse the command-line arguments
# opt_parser <- OptionParser(usage = "Usage: %prog [options]", option_list = option_list)
# opt <- parse_args(opt_parser)

# if (is.na(opt$input)){
#     stop("Input must be provided")
# }
# if (is.na(opt$output)) {
#     stop("Output must be provided")
# }
# # print(opt$input)

# # print(opt$output)
# # print(opt$binSize)

# # print(opt$rds)

library(ACE)
library(QDNAseq)
library(future)
library(ggplot2)

future::plan("multisession", workers=snakemake@threads)

input <- snakemake@input[['bam']]
binSize<-snakemake@params[['binsize']]
output<-snakemake@output[['plot']]

# print(input)
# print(binSize)
# print(output)


bins <-getBinAnnotations(binSize=binSize)


readCounts <- binReadCounts(bins, bamfiles=input, chunkSize=10e6)


readCountsFiltered <- applyFilters(readCounts, residual=TRUE, blacklist=TRUE)

readCountsFiltered <- estimateCorrection(readCountsFiltered)
png(filename=output, width=3000, height=1000, res=300)


plot(readCountsFiltered, ylim=c(10, 20))
dev.off()

