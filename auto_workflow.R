library(usethis)
library(devtools)
library(progress)
## data reformat
library(pbapply)
library(data.table)
library(plyr)
library(dplyr)
library(MSnbase)
library(igraph)
## visualize
library(ggplot2)
library(ggraph)
library(ggsci)
library(stringr)
library(grid)
## suggest
library(ChemmineR)
library(ChemmineOB)
library(rsvg)
library(grImport2)
library(ggimage)
## extra
library(gt)
library(classyfireR)
library(aplot)
library(ggalluvial)
##
load_all("~/MCnebula/R")
load_all("~/extra/R")
initialize_mcnebula(".")
# collate_structure(exclude_element = c("P", "S", "B", "Si"), fc = NA)
collate_structure()
build_classes_tree_list()
# collate_ppcp(min_possess = 50, max_possess_pct = 0.1, filter_via_struc_score = NA)
collate_ppcp()
# generate_parent_nebula(min_tanimoto = 0.5)
generate_parent_nebula()
generate_child_nebulae()
visualize_parent_nebula()
visualize_child_nebulae()

