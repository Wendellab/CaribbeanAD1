setwd(getwd())

library(gridExtra)
library(tidyverse)
library(ggrepel)
library(grid)
library(RColorBrewer)
library(ggrepel)
require(gridExtra)


cbbPalette <- brewer.pal(n = 12, name = "Paired")
#c("#000000", "#E69F00", "#56B4E9", "#009E73",  "#0072B2", "#D55E00", "#CC79A7")
###################black,  orange,    sky blue,  blulsh green, yellow, blue, vermilion, reddish purple 
scale_colour_manual(values=cbbPalette)

file1 <- list.files(pattern = "\\.eigenvec$")
file2 <- list.files(pattern = "\\.eigenval$")
file3 <- list.files(pattern = "\\.distmatrix$")
file4 <- list.files(pattern = "\\.mds$")



#######################################
pca <- read_table(file1, col_names = FALSE)
eigenval <- scan(file2)

pca <- pca[,-1]
pve <- data.frame(PC = 1:length(eigenval), pve = eigenval/sum(eigenval)*100)

a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
#a + ylab("Percentage variance explained") + theme_light()

# set names
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))

pca$population <- stringr::str_extract(pca$ind, "[^_]*_[^_]*")
#pca$species <- stringr::str_extract(pca$ind, "[^_]*")
pca$population <- gsub('OPB.*', 'OPB', pca$population)
pca$population <- gsub('MK.*', 'MK', pca$population)
pca$population <- gsub('AD1_GD_G.*', 'AD1_GD_G', pca$population)



b <- ggplot(pca, aes(PC1, PC2, col = population)) + geom_point(size = 2) +
  geom_text_repel(data=pca, max.overlaps=Inf,
                  aes(label= ind),  size = 4, segment.alpha =0,
                  show.legend = F)

b <- b + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))




b1 <- ggplot(pca, aes(PC1, PC3, col = population)) + geom_point(size = 2) +
  geom_text_repel(data=pca, max.overlaps=Inf,
                  aes(label= ind),  size = 4, segment.alpha =0,
                  show.legend = F) 


b1 <- b1 + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC3 (", signif(pve$pve[3], 3), "%)"))

write.table(pca[c(22, 1:21)], file = "GD_n90.eigenvec.3d",  row.names=F, col.names=F, quote = F, sep = "\t")
write.table(pve[2], file = "GD_n90.eigenval.3d",  row.names=F, col.names=F, quote = F, sep = "\t")

#########################################################################
#########################################################################

library(ape)
library(ggtree)
library(treeio)

tree_dis <- as.matrix(read.table(file3, row.names = 2)[-1])
tree_distree <- nj(tree_dis)

tree_distree2 <- as_tibble(tree_distree) %>% 
  mutate(Population = str_extract(label, "[^_]*_[^_]*")) %>%
  mutate(Population = gsub("OPB.*", "OPB", Population)) %>%
  mutate(Population = gsub("MK.*", "MK", Population)) %>%
  mutate(Population = gsub('AD1_GD_G.*', 'AD1_GD_G', Population)) %>%
  mutate(Population = gsub('AD1_Wild.*', 'AD1_Wild', Population)) %>%
  as.treedata()

njtree <- ggtree(tree_distree2, aes(color=Population), size = 0.9, layout="equal_angle", show.legend = F) +
  geom_tippoint(aes(shape=Population), size=4) +
  geom_text_repel(aes(label=label, color = Population), max.overlaps=Inf, size = 2, show.legend = F) +
  scale_shape_manual(values=seq(0,25))
  theme(#legend.position =  c(0.95, 0.2),
        legend.key = element_rect(fill = "white", colour = "grey30", linetype="dotted")) 

write.nexus(tree_distree, file = "nj.tre", translate = TRUE)


##############################################################################
##############################################################################
#mds2 <- read.table(file4, header=T)


#mds2$population <- stringr::str_extract(mds2$IID, "[^_]*_[^_]*")

#b3 <- ggplot(mds2, aes(C1, C2, col = population)) + geom_point(size = 2) +
#  geom_text_repel(data=mds2, max.overlaps=Inf,
#                  aes(label= IID),  size = 4, segment.alpha =0,
#                  show.legend = F)


##############################################################################
##############################################################################

pdf("PCA_tree.pdf", width = 15, height = 15)
b
b1
njtree
#b3
dev.off()


