args <- commandArgs(trailingOnly=TRUE)
library(jsonlite)
files <- sort(list.files(args[1], pattern="_template-object.rds$", full.names=TRUE))
stopifnot(length(files)==9)
rows <- lapply(files, function(f) {
  x <- readRDS(f)
  flat <- as.character(unlist(x, recursive=TRUE, use.names=FALSE))
  known_types <- c("no-changes","delin","ins-out","ins-in","dels-out","dels-in","subs")
  type_at <- which(flat %in% known_types)
  stopifnot(length(type_at)==1)
  prefix <- as.numeric(flat[seq_len(type_at-1)])
  ids <- flat[grepl("^SRR[0-9]+[.]",flat)]
  list(run=sub("_template-object.rds$","",basename(f)),
       component_lengths=as.integer(lengths(x)),
       native_first_count=as.numeric(x[[1]]),
       native_second_value=as.character(x[[2]]),
       type_position=type_at, type=flat[type_at],
       count_prefix=prefix, sum_count_prefix=sum(prefix),
       linked_id_count=length(ids), unique_linked_id_count=length(unique(ids)),
       scalar_interface_valid=(type_at==2))
})
write_json(rows,args[2],pretty=TRUE,auto_unbox=TRUE)
stopifnot(all(vapply(rows,function(x)x$sum_count_prefix==x$linked_id_count || x$sum_count_prefix==0,logical(1))))
