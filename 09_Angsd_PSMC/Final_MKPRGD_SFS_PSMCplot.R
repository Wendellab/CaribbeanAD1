setwd("C:/Users/weixuan/Desktop/GDPR_samples/00_FigureTables/07_Angsd/")
library(stringr)
library(dplyr)
library(grid)   # for the textGrob() function
library(ggpubr)
library(ggplot2)
library(scales)
library(reshape2)
library(pheatmap)
library(correlation)
library(scales)
library(ggcorrplot)
library(directlabels)
options(scipen=999)

#######################################################################

safe_colorblind_palette <- c("LR1" = "#F4A582",
                             "LR2" = "#D6604D",
                             "Cultivar" = "#FDDBC7",
                             "Wild" = "#B2182B",
                             "PR_CR" = "#01665E",
                             "PR_PCTW" = "#59A09A",
                             "PR_Ph" = "#80CDC1",
                             "PR_PR325" = "#35978F",
                             "PR_TBT" = "#C7EAE5",
                             "MK"= "black",
                             "GD" = "#4D9221",
                             "GD2" = "#7FBC41")
#######################################################################
###################################################################

psmc <- read.table("psmc_MKPRGD2.txt", header = F) %>%
  mutate(sample = str_extract(V6,  "(?<=_)[^_.]+")) %>%
  mutate(seq = ave(V1,V6,FUN = seq_along) ) #%>%
  #filter(!sample %in% c("CR", "PR325"))

unique(psmc$sample)

#https://www.geeksforgeeks.org/numbering-rows-within-groups-of-dataframe-in-r/
#https://stackoverflow.com/questions/46330352/means-multiple-columns-by-multiple-groups
psmc2 <- psmc %>%
  group_by(sample, seq) %>% 
  mutate(sample = gsub("CR", "PR_CR", sample), 
         sample = gsub("Ph", "PR_Ph", sample), 
         sample = gsub("PR325", "PR_PR325", sample), 
         sample = gsub("YUC", "GD", sample) ) %>% 
  summarise_at(vars("V1", "V2"), mean) %>%
  as.data.frame()

psmc2$sample <- factor(as.factor(psmc2$sample), 
                       levels = c("PR_Ph", "GD", "MK", "PR_PR325","PR_CR")) 


plotpsmc <- ggplot() + 
  geom_step(data=psmc, aes(x = V1, y = V2, group=V6), linewidth  = 0.2, color = 'black', alpha = 0.1) +
  geom_step(data=psmc2, aes(x = V1, y = V2, group = sample, color = sample), linewidth  = 2) +
  geom_dl(data=psmc2,aes(x = V1, y = V2, label = sample), method = list(cex = 1, dl.trans(x = x + 0.2), "last.points")) +
  ylim(0, 1.5)+
  scale_color_manual(values= safe_colorblind_palette) +
  xlab("Years (g=2, μ=4.56 x 10-9)") +
  ylab("Effective population size (x104)") +
  scale_x_log10(limits = c(1000,120000)) +
  theme_classic() +
  theme(legend.position =  "top",
        axis.title=element_text(size=12)) 


plotpsmc2 <- ggplot() + 
  geom_step(data=psmc, aes(x = V1, y = V2, group=V6), linewidth  = 0.2, color = 'black', alpha = 0.1) +
  geom_step(data=psmc2, aes(x = V1, y = V2, group = sample, color = sample), linewidth  = 2) +
  geom_dl(data=psmc2,aes(x = V1, y = V2, label = sample), method = list(cex = 1, dl.trans(x = x + 0.2), "last.points")) +
  ylim(0,30)+
  scale_color_manual(values= safe_colorblind_palette) +
  xlab("Years (g=2, μ=4.56 x 10-9)") +
  ylab("Effective population size (x104)") +
  scale_x_log10(limits = c(1000,10000000)) +
  #xlim(5000,10000000) +
  theme_classic() +
  theme(legend.position =  "top",
        axis.title=element_text(size=12)) 

plotpsmc$labels$group <- "Populations"
plotpsmc$labels$colour <- "Populations"
plotpsmc

plotpsmc2$labels$group <- "Populations"
plotpsmc2$labels$colour <- "Populations"
###################################################################
#SFS plot below with tWatterson caculated from Angsd###############
###################################################################

#https://thomasadventure.blog/posts/ggplot2-percentage-scale/
sfs.GD <- as.data.frame(scan("PR_GD.sfs")[c(2:22)]) %>%
  mutate(filename = "GD")  %>%
  mutate(site = as.numeric(rownames(.))) %>%
  rename("frequency" = 1) %>% 
  mutate(freq = (frequency / sum(frequency))*100) %>%
  mutate(theta = ifelse(site == max(site),
                        24 * (1 / site),
                        24 * (1 / site + 1 / (42 - site)))) %>%   
  mutate(thetafreq = theta / sum(theta) *100)  



sfs.CR <- as.data.frame(scan("PR_CR.sfs")[c(2:16)]) %>%
  mutate(filename = "PR_CR")  %>%
  mutate(site = as.numeric(rownames(.))) %>%
  rename("frequency" = 1) %>% 
  mutate(freq = (frequency / sum(frequency))*100) %>%
  mutate(theta = ifelse(site == max(site),
                        111 * (1 / site),
                        111 * (1 / site + 1 / (30 - site)))) %>%   
  mutate(thetafreq = theta / sum(theta) *100)  


sfs.PR325 <- as.data.frame(scan("PR_PR325.sfs")[c(2:16)]) %>%
  mutate(filename = "PR_PR325")  %>%
  mutate(site = as.numeric(rownames(.))) %>%
  rename("frequency" = 1) %>% 
  mutate(freq = (frequency / sum(frequency))*100) %>%
  mutate(theta = ifelse(site == max(site),
                        96.7 * (1 / site),
                        96.7 * (1 / site + 1 / (30 - site)))) %>%    
  mutate(thetafreq = theta / sum(theta) *100) 

#sfs.Ph <- as.data.frame(scan("PR_Ph.sfs")[c(2:6)]) %>%
#  mutate(filename = "PR_Ph")  %>%
#  mutate(site = as.numeric(rownames(.))) %>%
#  rename("frequency" = 1) %>% 
#  mutate(freq = (frequency / sum(frequency))*100) %>%
#  mutate(theta=10624986* (1/site + 1/(10 - site))) %>% 
#  mutate(thetafreq = theta / sum(theta) *100) 

sfs.MK <- as.data.frame(scan("PR_MK.sfs")[c(2:26)]) %>%
  mutate(filename = "MK")  %>%
  mutate(site = as.numeric(rownames(.))) %>%
  rename("frequency" = 1) %>% 
  mutate(freq = (frequency / sum(frequency))*100)  %>% 
  mutate(theta = ifelse(site == max(site),
                        79.7 * (1 / site),
                        79.7 * (1 / site + 1 / (50 - site)))) %>%  
  mutate(thetafreq = theta / sum(theta) *100) 
rm(sfs.Ph)

theta=79.7* (1/25 )

sum(as.data.frame(scan("PR_GD.sfs")[c(2:22)]))
sum(as.data.frame(scan("PR_CR.sfs")[c(2:16)]))
sum(as.data.frame(scan("PR_PR325.sfs")[c(2:16)]))
sum(as.data.frame(scan("PR_Ph.sfs")[c(2:6)]))
sum(as.data.frame(scan("PR_MK.sfs")[c(2:25)]))
#####################

sfs_names <- ls(pattern = "^sfs\\.")
# Retrieve the data frames and bind them
sfs_all <- bind_rows(mget(sfs_names))

sfs_all$filename <- factor(as.factor(sfs_all$filename), 
                           levels = c("MK", "GD", "PR_CR","PR_PR325"))

plotsfs <- ggplot(data = sfs_all, aes(x=site, y=freq, fill = filename)) +
  geom_bar(stat = "identity") +
  geom_line(aes(y = thetafreq), color = "red", size = 1) +
  facet_grid(~ filename, 
             scales = "free_x", switch = 'x') +
  #scale_y_continuous(labels = percent_format(scale = 1)) +
  ylab ("Proportion of SNP sites (%)") +
  xlab ("Minor allele frequency") +
  scale_fill_manual(values= safe_colorblind_palette) +
  theme_minimal() +
  theme(legend.position = "none",
        strip.placement = "outside",
        #axis.title = element_blank(),
        plot.title = element_text(size = 12))

plotsfs2 <- ggplot(data = sfs_all, aes(x=site, y=frequency, fill = filename)) +
  geom_bar(stat = "identity") +
  #geom_line(aes(y = thetafreq), color = "red", size = 1) +
  facet_grid(~ filename, 
             scales = "free_x", switch = 'x') +
  #scale_y_continuous(labels = percent_format(scale = 1)) +
  ylab ("SNP sites") +
  xlab ("Minor allele frequency") +
  scale_fill_manual(values= safe_colorblind_palette) +
  theme_minimal() +
  theme(legend.position = "none",
        strip.placement = "outside",
        #axis.title = element_blank(),
        plot.title = element_text(size = 12))


###################################################################
#Tajima's D plot below with facet #################################
###################################################################

library(dplyr)
library(ggridges)
library(forcats)

file_list <- list.files(path = "C:/Users/weixuan/Desktop/GDPR_samples/00_FigureTables/07_Angsd/", pattern = "*.pestPG", full.names = TRUE)

read_file_with_filename <- function(file) {
  df <- read.table(file, sep="\t",header=T)
  df$filename <- basename(file)  
  return(df)}

combined_df <- file_list %>%
  lapply(read_file_with_filename) %>%
  bind_rows()%>%
  mutate(filename = gsub ("PR_MK", "MK", filename))%>%
  mutate(filename = gsub ("PR_GD", "GD", filename))%>%
  mutate(filename = gsub (".thetasWindow.gz.pestPG", "", filename))%>%
  mutate(filename = fct_rev(fct_infreq(filename)))

combined_df$filename <- factor(as.factor(combined_df$filename), 
                     levels = c("PR_PR325","PR_CR","PR_Ph", "GD", "MK")) 


combined_df%>%
  group_by(filename) %>%
  summarise(mean_value = mean(Tajima, na.rm = TRUE))

tajamadplot_PR <-
  ggplot(combined_df %>%
           filter(!filename %in% c("PR_Ph")),
         aes(x = Tajima, y = filename, group = filename)) + 
  geom_density_ridges(aes(fill = filename), alpha = 0.5,
                      quantile_lines = T, quantile_fun = mean, show.legend = F) +
  geom_vline(xintercept =0, color = 'blue') +
  scale_fill_manual(values= safe_colorblind_palette) +
  xlab ("Tajima's D") +
  theme_minimal() +
  theme(axis.title.y = element_blank())

#############################################################
#############################################################
library(cowplot)



finalplot <-
  ggdraw() +
  draw_plot(tajamadplot_PR, x = 0, y = 0, width = 0.3, height = 1) +
  draw_plot(plotsfs, x = 0.3, y = 0.5, width = 0.7, height = 0.5) +
  draw_plot(plotpsmc, x = 0.3, y = 0, width = 0.7, height = 0.5)  +
  draw_plot_label(label = c("A", "B", "C"), size = 15,fontface = "bold",
                  x = c(0, 0.3, 0.3), y = c(1, 1, 0.5))


pdf("Fig7_MKPRGD_SFS_PSMC_final.pdf", width = 11, height = 8)
finalplot
dev.off()


finalplot2 <-
  ggdraw() +
  draw_plot(plotsfs2, x = 0, y = 0.5, width = 1, height = 0.5) +
  draw_plot(plotpsmc2, x = 0, y = 0, width = 1, height = 0.5) +
  draw_plot_label(label = c("A", "B"), size = 15,fontface = "bold",
                  x = c(0, 0), y = c(1, 0.5))

pdf("FigS8_MKPRGD_PSMC.pdf", width =10, height = 8)
finalplot2
dev.off()







##########################################

smc <- read.csv("../sfs_psmc/joint.csv", header = T, sep = ',')
#smc_wild <- read.csv("smc_plot_wild.csv", header = T, sep = ',')

smc_mk_wild <- rbind(smc_wild,smc)

colnames(smc) <- c("Population/Group", "Year", "Ne", "type" , "num")

breaks = c(2000, 3000,  5000, 7000, 10000, 50000)
labels = as.character(breaks)

smc.plot <- ggplot(smc, aes(x=Year, y=Ne)) +
#  geom_rect(aes(xmin=9000, xmax=11000, ymin=0, ymax=Inf), alpha = 0.01, color = "black",
#            fill = "grey")+ 
  geom_line(aes(color=`Population/Group`), size =2, alpha = 0.7) +
#  geom_point(aes(color=label)) + 
  scale_y_continuous(trans='log10',labels = scales::comma) + 
  
#  scale_x_continuous(trans='log10',labels = trans_format('log10',math_format(10^.x))) +
  scale_x_continuous(trans='log10', breaks = breaks) +
  scale_color_manual(values = c("#D55E00", "grey20"))+
  theme_bw() +  
  theme(legend.position =  c(0.13, 0.8),
        #legend.title=element_blank(),
        legend.key = element_rect(fill = "white", colour = "grey30", linetype="dotted")) 

smc.plot