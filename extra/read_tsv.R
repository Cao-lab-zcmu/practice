read_tsv <- function(path){
  file <- data.table::fread(input=path, sep="\t", header=T, quote="", check.names=F)
  return(file)
}
write_tsv <-
  function(x, filename){
    write.table(x, file = filename, sep = "\t", col.names = T, row.names = F, quote = F)
  }
