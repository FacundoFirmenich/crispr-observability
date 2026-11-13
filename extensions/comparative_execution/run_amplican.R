args <- commandArgs(trailingOnly=TRUE)
w <- args[1]
library(amplican)
stopifnot(as.character(packageVersion('amplican')) == '1.32.1')
sink(file.path(w,'evidence','R_SESSION.txt'));sessionInfo();print(formals(amplicanPipeline));sink()
out <- file.path(w,'evidence','native');dir.create(out)
amplicanPipeline(file.path(w,'amplican.csv'),file.path(w,'inputs'),out,
 knit_reports=FALSE,average_quality=0,min_quality=0,fastqfiles=1,
 normalize=FALSE,continue=TRUE)
