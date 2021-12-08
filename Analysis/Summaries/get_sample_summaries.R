# get summaries for each MB sample against all DNA samples, one chr
library(stringr)

snpcalls_dir <- "SnpCalls"
readcounts_dir <- "Readcounts"
out_dir <- "Summaries"
chr <- "1"

# load the corresponding imputed SNPs for that chromosome
load(paste0(snpcalls_dir, "/imp_snp_", chr, ".RData")) # loads as `imp_snps`
# snpinfo$pos_bp <- round(snpinfo$pos_Mbp * 1e6)

files <- list.files(readcounts_dir, pattern=paste0("_", chr, "_"))
file.chunks.df <- str_split_fixed(files, "_", n=6)
samples <- str_extract(files, paste0("(.*)(?=_", chr, "_)"))
mb_ids <- paste(file.chunks.df[,1], file.chunks.df[,2], file.chunks.df[,3], sep="-")

# check that all IDs with imp_snps data
stopifnot(all(mb_ids %in% rownames(imp_snps)))

# initialize output
sample_results <- pair_results <- vector("list", length(files))
names(sample_results) <- names(pair_results) <- samples

for(indi in seq_along(files)) {
    cat("file", indi, "of", length(files), "\n") 
    file <- files[indi]
    sample <- samples[indi]
    mb_id <- mb_ids[indi]

    # read the read counts
    readcounts <- readRDS(file.path(readcounts_dir, file))

    # find the reads' SNP positions in the snpinfo table
    snpinfo_row <- match(readcounts$pos, snpinfo$pos)
    
    # should all have been found
    # stopifnot(!any(is.na(snpinfo_row)))

    # remove NAs
    snpinfo_row <- snpinfo_row[!is.na(snpinfo_row)]
    
    # get SNP IDs since that's what's in the columns of imp_snp
    imp_snps_col <- snpinfo$snp_id[snpinfo_row]
    
    # subset to SNPs present in the columns of imp_snps
    imp_snps_col <- imp_snps_col[imp_snps_col %in% colnames(imp_snps)]

    # create object to contain the results for the single samples
    sample_results[[indi]] <- array(0, dim=c(nrow(imp_snps), 3, 2))
    dimnames(sample_results[[indi]]) <- list(rownames(imp_snps), c("AA", "AB", "BB"), c("A", "B"))

    for(i in 1:nrow(imp_snps)) {
        g <- imp_snps[i, imp_snps_col]
        for(j in 1:3) {
            sample_results[[indi]][i,j,1] <- sum(readcounts$count1[!is.na(g) & g==j])
            sample_results[[indi]][i,j,2] <- sum(readcounts$count2[!is.na(g) & g==j])
        }
    }

    # create object to contain the results for sample pairs
    pair_results[[indi]] <- array(0, dim=c(nrow(imp_snps), 3, 3, 2))
    dimnames(pair_results[[indi]]) <- list(rownames(imp_snps), c("AA", "AB", "BB"),
                                          c("AA", "AB", "BB"), c("A", "B"))

    g0 <- imp_snps[mb_id, imp_snps_col]
    for(i in 1:nrow(imp_snps)) {
        g <- imp_snps[i, imp_snps_col]
        for(j in 1:3) {
            for(k in 1:3) {
                pair_results[[indi]][i,j,k,1] <- sum(readcounts$count1[!is.na(g0) & g0==j & !is.na(g) & g==k])
                pair_results[[indi]][i,j,k,2] <- sum(readcounts$count2[!is.na(g0) & g0==j & !is.na(g) & g==k])
            }
        }
    }

} # loop over MB samples

saveRDS(sample_results, file.path(out_dir, paste0("sample_results_chr", chr, ".rds")))
saveRDS(pair_results, file.path(out_dir, paste0("pair_results_chr", chr, ".rds")))
