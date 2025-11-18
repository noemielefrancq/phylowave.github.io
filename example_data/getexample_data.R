## Get exmaple data for phylowave website
library(ape)
library(stringr)
setwd('~/Dropbox/Projects/202312_Index_paper_methods_Sc2_H3N2_BP_TB/Phylowave_Learning-fitness-dynamics-pathogens-in-phylogenies/1_Data/')
path_to_save = '~/Dropbox/Web_page/phylowave.github.io/example_data/'

## SC2 
tree = read.nexus('1_1_SARS_CoV_2/Tree_SARSCoV2_global_alltime_nextstrain_20230414.nexus')
# tree$tip.label = sapply(tree$tip.label, function(x)str_replace(x, pattern = '\\(', replacement = '-'))
# tree$tip.label = sapply(tree$tip.label, function(x)str_replace(x, pattern = '\\)', replacement = ''))
metadata = read.csv('1_1_SARS_CoV_2/List_all_isolates_SARSCoV2_global_20230414.csv')
metadata = metadata[,-1]
metadata = metadata[,-2]
metadata$Isolate = tree$tip.label
metadata$Collection_time = as.numeric(sapply(tree$tip.label, function(x)tail(str_split(x, pattern = '/')[[1]],2)[1]))
metadata$Country = sapply(tree$tip.label, function(x)str_split(x, pattern = '/')[[1]][2])
metadata$Nextstrain_clade = sapply(tree$tip.label, function(x)tail(str_split(x, pattern = '/')[[1]],1))
write.tree(tree, paste0(path_to_save, 'Example_SC2.newick'))
write.csv(metadata, paste0(path_to_save, 'Example_SC2.csv'), row.names = F, quote = F)

## H3N2 
tree = read.nexus('1_2_H3N2/Tree_H3N2_ha_global_12y_nextstrain_20230414.nexus')
metadata = read.csv('1_2_H3N2/List_all_isolates_H3N2_global_20230414.csv')
metadata = metadata[,-1]
metadata$Isolate = tree$tip.label
write.tree(tree, paste0(path_to_save, 'Example_H3N2.newick'))
write.csv(metadata, paste0(path_to_save, 'Example_H3N2.csv'), row.names = F, quote = F)

## Pertussis 
tree = read.nexus('1_3_Pertussis/Tree_Bpertussis_France_20230901.tree')
metadata = read.csv('1_3_Pertussis/List_all_isolates_Bpertussis_France.csv')
metadata = metadata[,-1]
metadata$Isolate = tree$tip.label
metadata$Source = unlist(lapply(metadata$Source, function(x)paste0(str_split(x, ',')[[1]], collapse = '')))
write.tree(tree, paste0(path_to_save, 'Example_BP.newick'))
write.csv(metadata, paste0(path_to_save, 'Example_BP.csv'), row.names = F, quote = F)

## TB
tree = read.tree('1_4_TB/Tree_Mtuberculosis_Samara_Russia_20230901.tree')
tree$node.label = NULL
metadata = read.csv('1_4_TB/List_all_isolates_TB_Samara_Russia_20231216.csv')
metadata = metadata[,-1]
tree$tip.label[which(tree$tip.label == "NC_000962.3")] = 'H37Rv'
write.tree(tree, paste0(path_to_save, 'Example_TB.newick'))
write.csv(metadata, paste0(path_to_save, 'Example_TB.csv'), row.names = F, quote = F)

# which(is.na(match(tree$tip.label, metadata$Isolate)))
# tree$tip.label[which(is.na(match(tree$tip.label, metadata$Isolate)))
# 
# 
# 
# str_replace(tree$tip.label[3000], pattern = '\\(', replacement = '-')
# tree$tip.label = sapply(tree$tip.label, function(x)str_replace(x, pattern = '\\(', replacement = '-'))
# tree$tip.label = sapply(tree$tip.label, function(x)str_replace(x, pattern = '\\)', replacement = ''))



Q = matrix(0.26, nrow = 2, ncol = 2)
diag(Q) = -0.26
rowSums(Q)
expm::expm(Q*1) 



