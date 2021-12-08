library(RSQLite)
library(qtl2)

# get SNP calls for one chromosome
chr <- "1"

# data related to the SNP arrays are in /z/Proj/attie/kbroman/AttieDO/
# dir <- "/z/Proj/attie/kbroman/AttieDO"

sql_file <- "cc_variants.sqlite"
out_dir <- "SnpCalls"

#########
### Accessing the genotype probabilities
########
# genotype probabilities (as arrays n_ind x 36 x n_markers) are in .../CalcGenoProb/attieDO_probs*.rds
# where * = 1, 2, ..., 19, X

# read in chromosome
# pr <- readRDS( file.path(dir, paste0("CalcGenoProb/attieDO_probs", chr, ".rds")) )

# I am using allele probabilities, rather than genotype probabilities because this is all I could find
# creates `apr` object, names of this object are the chromosomes
# under each chromosome is an array with size 946 x 8 x 5000-ish
load("prob.8state.allele.qtl2_191121.Rdata")

# we also need the physical map of the markers in these probabilities
# pmap <- readRDS( file.path(dir, "CalcGenoProb/attieDO_pmap.rds") )

# pmap has columns marker, chr, and pos, where pos is in Mb
tmp.pmap <- read.csv("gigamuga_pmap_190510.csv")

# convert position from Mb to bp
tmp.pmap$pos <- round(tmp.pmap$pos*1e6)

# convert pmap into the form expected by index_snps
pmap <- list()
for (ii in unique(tmp.pmap$chr)) {
  tmp.pmap.this.chr <- tmp.pmap[tmp.pmap$chr == ii, ]
  pmap.this.chr <- tmp.pmap.this.chr$pos
  names(pmap.this.chr) <- tmp.pmap.this.chr$marker
  pmap[[ii]] <- pmap.this.chr
}

#########
### Accessing the founder SNP genotypes
########
# genotypes for *all* SNPs in the genome for the 8 founder lines are in an SQLite database
#    .../DerivedData/ccfoundersnps.sqlite
# sql_file <- file.path(dir, "DerivedData/ccfoundersnps.sqlite")

db <- dbConnect(SQLite(), sql_file) # connect to database

# now, to grab ALL of the SNPs on the chr
snpinfo <- dbGetQuery(db, paste0("SELECT * FROM variants WHERE chr=='", chr, "'"))
dbDisconnect(db)

# the alleles column is like "A|G" where A is the major allele and G is the minor allele
# some "complex" SNPs will be like "A|G/C" where there are multiple possiblities for the minor allele
# let's drop those
complex <- grepl("/", snpinfo$alleles, fixed=TRUE)
snpinfo <- snpinfo[!complex,]


#########
### Getting imputed SNP genotype probabilities for these SNPs for all mice
########
# we first need to create indexes connecting the SNP locations to the physical map for the genotype probabilities
#    - this returns only a portion of the snps...within the region defined by markers in pmap
#    - also, uses an indexing system for focus on subset of snps with distinct patterns
# this adds some columns to the "snps" object
snpinfo <- index_snps(pmap, snpinfo)

# Now convert the genotype probabilities to SNP probabilities
snp_pr <- genoprob_to_snpprob(apr, snpinfo)

# we can use maxmarg() in the qtl2 package to get imputed SNP genotypes
# also drop the surrounding list business and just make it a matrix
imp_snps <- maxmarg(snp_pr)[[1]]

# CRITICAL CHANGE!
# Because I used allele probabilities, rather than genotype probabilities,
# the current object will fail in get_sample_summaries.R.
# To get around this, I need to recode imp_snps:
# -- 2 becomes 3 (homozygous for minor allele)
# -- NA becomes 2 (heterozygous)
# -- 1 stays as 1 (homozygous for major allele)
imp_snps[imp_snps==2] <- 3
imp_snps[is.na(imp_snps)] <- 2

# pull out the alleles from the SNP table
alleles <- strsplit(snpinfo$alleles, "\\|")
snpinfo$allele1 <- sapply(alleles, "[", 1)
snpinfo$allele2 <- sapply(alleles, "[", 2)

# save as .RData file
save(imp_snps, snpinfo, file=file.path(out_dir, paste0("imp_snp_", chr, ".RData")))
