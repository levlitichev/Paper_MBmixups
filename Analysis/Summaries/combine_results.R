in_dir <- "Summaries"
out_dir <- "Summaries"

pair_results <- sample_results <- NULL
for(chr in c(1,12,13)) {
# for(chr in 1:19) {
    cat(chr, "\n")

    file <- file.path(in_dir, paste0("pair_results_chr", chr, ".rds"))
    if(file.exists(file)) {
        x <- readRDS(file)
        if(is.null(pair_results)) {
            pair_results <- x
        } else {
            stopifnot(length(x) == length(pair_results))
            x <- x[names(pair_results)]

            for(i in seq_along(pair_results)) {
                stopifnot(dim(x[[i]]) == dim(pair_results[[i]]),
                          all(rownames(x[[i]]) == rownames(pair_results[[i]])))
                pair_results[[i]] <- pair_results[[i]] + x[[i]]
            }
        }
    } else {cat(file, "does not exist")}

    file <- file.path("Summaries", paste0("sample_results_chr", chr, ".rds"))
    if(file.exists(file)) {
        x <- readRDS(file)
        if(is.null(sample_results)) {
            sample_results <- x
        } else {
            stopifnot(length(x) == length(sample_results))
            x <- x[names(sample_results)]

            for(i in seq_along(sample_results)) {
                stopifnot(dim(x[[i]]) == dim(sample_results[[i]]),
                          all(rownames(x[[i]]) == rownames(sample_results[[i]])))
                sample_results[[i]] <- sample_results[[i]] + x[[i]]
            }
        }
    }
}

saveRDS(pair_results, file=file.path(out_dir, "pair_results_allchr.rds"))
saveRDS(sample_results, file=file.path(out_dir, "sample_results_allchr.rds"))
