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

