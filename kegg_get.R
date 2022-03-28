library(KEGGREST)
library(org.Mm.eg.db)
library(clusterProfiler)
#BiocManager::install("org.Hs.eg.db")
#org <- keggList('organism')
### human <- hsa
### mice <- mmu
### rats <- rno
# save
# load
species="mmu"
pathway <- keggLink("pathway", species)
pathway <- unique(pathway)
# capture.output(mmu_test,file="test")
pathway.database <- vector(mode = "list", length = length(pathway))
for(i in 1:length(pathway)){
 	cat(i," >>> ",pathway[[i]],"\n")
 	pathway.database[[i]] <- keggGet(dbentries = pathway[[i]]) }
############ compound message
compound.id <- keggList('compound')
compound.id <- names(compound.id)
compound.database <- vector(mode = 'list', length = length(compound.id))
for(i in 1:length(compound.id)){
 	cat(i," >>> ",compound.id[i],"\n")
 	compound.database[[i]] <- keggGet(dbentries = compound.id[i]) }
#### search <- keggGet(entries = object)
#### 
## export DBLINKS
sink("~/operation/re_fecal_neg/kegg/dblink.list")
for(i in 1:length(compound.database)){
  cat("BEGIN_COMPOUND\n")
  cat(compound.database[[i]][[1]]$ENTRY, "\n")
  cat("DBLINKS\n")
  for(j in 1:length(compound.database[[i]][[1]]$DBLINKS)){
    cat(compound.database[[i]][[1]]$DBLINKS[j], "\n")
  }
  cat("NAME\n")
  for(j in 1:length(compound.database[[i]][[1]]$NAME)){
    cat(compound.database[[i]][[1]]$NAME[j], "\n")
  }
  cat("END_COMPOUND\n")
  cat("\n")
}
sink()

####################################################
## kegg api
png <- keggGet("mmu01100", "image")
library(png)
writePNG(png, "test.png")
####################################################
query=keggGet(c("mmu:01100"))

