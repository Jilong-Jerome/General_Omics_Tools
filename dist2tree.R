library(ape)
library(dplyr)
library(readr)
library(tidyr)
args <- commandArgs(T)
distance_all <- read_tsv(args[1],col_names = F)
for (i in c(1:(NROW(distance_all)/64))) {
  distance <- distance_all[(64*(i-1)+1):(64*i),]
  distance_cut <- distance[,2:65]
  names(distance_cut) <- distance[,1] %>% unlist()
  tre <- nj(as.dist(distance_cut))
  write.tree(tre, file = paste(args[1],".newick",sep = ""), append = TRUE,
             digits = 10, tree.names = FALSE)
}
