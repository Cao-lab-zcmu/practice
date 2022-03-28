library(XML)

xml <- xmlToDataFrame("sites.xml")

write.table(xml, file = paste0("data",".tsv"), quote = FALSE, append = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)









