library(LEA)
library(tidyverse)
library(hues)


setwd(getwd())

file1 <- list.files(pattern = "\\.ped$")

ped2geno(file1)
ped2lfmm(file1)

file2 <- list.files(pattern = "\\.lfmm$")
pc=pca(file2,scale=TRUE)
tw=tracy.widom(pc)


pdf("LEA1_pca_subset.pdf", width = 10, height = 10)
#plotcross-entropycriterionforallrunsinthesnmfproject 
plot(tw$percentage,pch =19,col="darkblue",cex=.8)
dev.off()


file3 <- list.files(pattern = "\\.geno$")
project=NULL 
project=snmf(file3, K=1:10, entropy= TRUE, repetitions=10, project= "new")

pdf("LEA2_crvalues.pdf", width = 10, height = 10)
#plotcross-entropycriterionforallrunsinthesnmfproject 
plot(project,col="blue",pch=19,cex=1.2)
dev.off()

summary.info = summary(project)$crossEntropy

pdf("LEA3_crvalues.pdf", width = 10, height = 10)
as.data.frame(t(summary.info)) %>%
  tibble::rownames_to_column("K") %>%
  mutate(K = factor(as.numeric(substring(K, 4)))) %>%
  ggplot() +
  geom_hline(aes(yintercept = min(summary.info[2,]), color = "red")) +
  geom_pointrange(aes(x = K, ymin=min, y = mean, ymax=max)) +
  labs(title = "Cross-entropy versus K", x = "Number of ancestral populations (K)", y = "Cross-entropy") +
  theme_minimal() +
  theme(legend.position = "none",
        panel.grid.major.x=element_blank(),
        panel.grid.minor.x=element_blank())
dev.off()


lapply(1:10, cross.entropy, object=project) %>%
  sapply(which.min) %>%
  as.data.frame(nm="run") %>%
  rowid_to_column("K") %>%
  mutate(cross.entropy = mapply(cross.entropy, K=K, run=run,
                                MoreArgs = list(object=project))) %>%
  write.csv(.,file = 'LEA4_crvalues.txt', quote = F, row.names = F)

save.image()


