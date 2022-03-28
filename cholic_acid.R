cat("Cholic acid formula: C24-H40-O5\n")
cat("Cholic acid-d4 formula: C24-H36-D4-O5\n")
cat("---Ion mode---\n")

cat("Cholic acid-d4 [M-H]-: C24-H35-D4-O5\n")
d4_neg <- 12*24 + 1.007825*35 + 2.014102*4 + 15.994915*5
cat(paste0("Cholic acid-d4 in negtive ion mode [M-H]-: ", d4_neg, "\n"))

cat("Cholic acid-d4 [M-2*H2O+H]+: C24-H33-D4-O3\n")
d4_pos <- 12*24 + 1.007825*33 + 2.014102*4 + 15.994915*3
cat(paste0("Cholic acid-d4 in positive ion mode [M-H2O+H]+: ", d4_pos, "\n"))


