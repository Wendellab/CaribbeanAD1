setwd("C:/Users/weixuan/Desktop/GDPR_samples/00_FigureTables/05_ROH_LD/")

library(dplyr)
library(ggplot2)
library(directlabels)
library(stringr)
#######################################################################

safe_colorblind_palette <- c("LR1" = "#F4A582",
                             "LR2" = "#D6604D",
                             "Cultivar" = "#FCAB86",
                             "Wild" = "#B2182B",
                             "PR_CR" = "#01665E",
                             "PR_PCTW" = "#59A09A",
                             "PR_Ph" = "#80CDC1",
                             "PR_PR325" = "#35978F",
                             "PR_TBT" = "#C7EAE5",
                             "MK"= "black",
                             "GD" = "#4D9221",
                             "GD2" = "#7FBC41")

safe_shape_palette <- c("MK"= 3,
                        "Wild" = 13,
                        "LR1" = 17,
                        "Cultivar" = 17,
                        "LR2" = 17,
                        "GD" = 19,
                        "GD2" = 19,
                        "PR_CR" = 15,
                        "PR_PCTW" = 15,
                        "PR_Ph" = 15,
                        "PR_PR325" = 15,
                        "PR_TBT" = 15)

#######################################################################
LD_dec <- read.table("LD_Fig_test.txt", header = T)
LD_dec$Mean_r.2 <- as.numeric(LD_dec$Mean_r.2)
LD_dec$Dist <- as.numeric(LD_dec$Dist)
LD_dec$Accession <- gsub("[.].*","",LD_dec$name_file)

LD_dec.avg <- LD_dec %>%
  group_by(Accession, Dist) %>%
  dplyr::summarize(Mean.Mean_r.2 = mean(Mean_r.2, na.rm=TRUE)) %>%
  as.data.frame()

LD_dec_plot <- ggplot(LD_dec.avg, aes(x=Dist/1000, y=Mean.Mean_r.2)) +
  geom_line(aes(color=Accession), size = 1)+
  scale_y_continuous(limits=c(0, 1), breaks = seq(0, 1, 0.2)) +
  geom_dl(aes(label = Accession,color=Accession), method = list(cex = 0.5, dl.trans(x = x + 0.2), "last.points")) +
  ylab(expression(r^{2})) +
  xlab("Distance (kb)") +
  theme_bw() +
  scale_x_continuous(limits=c(0, 550), breaks = seq(0,500,100)) +
  scale_color_manual(name = "Population/Group", values=safe_colorblind_palette)  +
  theme(legend.position = "none",
        legend.key = element_rect(fill = "white", colour = "grey30", linetype="dotted")) 
LD_dec_plot

#pdf("Fig4_GD_LD_n121.pdf", width = 7, height = 7)
#LD_dec_plot
#dev.off()

########################################################
########################################################
########################################################
library(detectRUNS)
load("ROH_n86.final.RData")

slidingRuns$id <- gsub("AD1_GD_G1A","AD1_GD2_G1A",slidingRuns$id)
slidingRuns$id <- gsub("AD1_GD_G23A","AD1_GD2_G23A",slidingRuns$id)
slidingRuns$id <- gsub("AD1_GD_G24A","AD1_GD2_G24A",slidingRuns$id)
slidingRuns$id <- gsub("AD1_GD_G25B","AD1_GD2_G25B",slidingRuns$id)
slidingRuns$id <- gsub("Pop","Site",slidingRuns$id) 

slidingRuns$group <- slidingRuns$id

slidingRuns$group <- gsub("AD1_","",slidingRuns$group) 
slidingRuns$group <- gsub("_peru","",slidingRuns$group)
slidingRuns$group <- gsub("^(Wild.*?_.*?)_(.*)$", "\\1", slidingRuns$group)

unique(slidingRuns$group)
#https://stackoverflow.com/questions/54064503/custom-factor-levels-in-a-concatenated-string
x6 <- as.factor(slidingRuns$group)
i1 <- match(sub("_[^_]+$*", "", levels(as.factor(slidingRuns$group))), 
            c("Cultivar", "PR_PCTW" ,"LR2", "LR1","PR_TBT", "GD", "GD2", "MK_Site1","MK_Site2","MK_Site3", 
              "PR_CR", "PR_PR325", "PR_Ph","Wild", "AD2_Wild","AD4_mus"))
slidingRuns$group <- factor(x6, levels = levels(x6)[seq_along(levels(x6))[order(i1)]])

unique(slidingRuns$id)




###########################

#https://forum.posit.co/t/multiple-ifelse-statement-with-in/101363
slidingRuns2 <- slidingRuns %>%
  mutate(category = ifelse(lengthBps < 1000000, "0.25-1 Mb", 
                           ifelse(lengthBps %in% 1000000:2000000, "1-2 Mb",
                                         ifelse(lengthBps > 2000000, ">2 Mb",NA))))

slidingRuns3 <- slidingRuns2 %>%
  group_by(category, group) %>% 
  dplyr::summarize(lengthBps_sum = sum(lengthBps, na.rm=TRUE)/2296245394) %>%
  as.data.frame()

category_fac <- as.factor(slidingRuns3$category)
slidingRuns3$category <- factor(category_fac, levels = c(">2 Mb", "1-2 Mb" ,"0.25-1 Mb"))


ROH_plot <- ggplot(slidingRuns3 %>%
                     filter(!grepl("AD2|AD4", group)), 
                   aes(fill=category, x=group, y=lengthBps_sum )) + 
  geom_bar(position="stack", stat="identity", width = 0.7) +
  scale_fill_manual(values=c("0.25-1 Mb" = "steelblue",
                             "1-2 Mb" = "#E69F00",
                             ">2 Mb" = "#009E73")) +
  geom_vline(xintercept = c(30.5, 51.5, 76.5, 111.5)) +
  #  xlab("Populations/Groups") +
  ylab("Propotion of total genome") +
  theme_classic()+
  theme(legend.position =  c(0.96, 0.8),
        #legend.box.background = element_rect(linetype="dotted", colour = "black"),
        #legend.title=element_blank(),
        axis.title.x  = element_blank(),
        legend.key = element_rect(fill = "white", colour = "grey30"),
        axis.text.x = element_blank()
        #axis.text.x = element_text(angle = -90, hjust= 0.001, vjust=0.5)
        )

ROH_plot$labels$fill <- "Length of ROH"
ROH_plot
#####################################################################
#####################################################################

slidingRuns4 <- slidingRuns2 %>%
  group_by(group) %>% 
  dplyr::summarize(Froh = sum(lengthBps, na.rm=TRUE)/2296245394) %>%
  mutate(group2 = gsub("_[^_]+$*", "", group))%>%
  mutate(group2 = gsub("_Site.*", "", group2))%>%
  filter(! group2 %in% c("PR_TBT", "AD2_Wild", "AD4_mus", "PR_PCTW", "GD2"))%>%
  as.data.frame()

group_ordered <- with(slidingRuns4,  reorder(group2,  Froh,  mean))
slidingRuns4$group2 <- factor(slidingRuns4$group2, levels = levels(group_ordered))


Roh_Froh_Plot <- ggplot(slidingRuns4, aes(x=group2, y=Froh, fill=group2)) + 
  geom_boxplot(show.legend = F, width = 0.7, outlier.shape = NA) + 
  geom_jitter(color="black", size=1, alpha=0.9) +
  geom_text(aes(label = formatC(after_stat(y), format = "f", digits  = 3), group = group2), size = 4,
            stat = 'summary', fun = mean,  nudge_y = 0.1 ) +
  stat_summary(fun.y=mean, geom="point", shape=18, size=3, color="blue", fill="blue") +
  scale_fill_manual(values=safe_colorblind_palette)  +
  xlab("Populations/Groups") +
  ylab(expression("F"[ROH]*" inbreeding")) +
  theme_classic()+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = -45, v = -1,  size = 9),
        panel.border = element_rect(colour = "black", fill=NA))


Roh_Froh_Plot

#####################################################################
####################################################################
####################################################################


heteperce  <- read.csv("MKGDPRYuan_n121.het.het", sep = '\t')
heteperce$INDV<- gsub("AD1_GD_G1A","AD1_GD2_G1A",heteperce$INDV)
heteperce$INDV<- gsub("AD1_GD_G23A","AD1_GD2_G23A",heteperce$INDV)
heteperce$INDV<- gsub("AD1_GD_G24A","AD1_GD2_G24A",heteperce$INDV)
heteperce$INDV<- gsub("AD1_GD_G25B","AD1_GD2_G25B",heteperce$INDV)

heteperce$group <- stringr::str_extract(heteperce$INDV, "[^_]*_[^_]*")
heteperce$group <- gsub("AD1_","",heteperce$group) 


heteperce2 <- heteperce %>%
  mutate(group = str_extract(INDV, "[^_]*_[^_]*_[^_]*"))%>%
  mutate(group = gsub("AD1_","", group)) %>%
  mutate(group2 = str_extract(INDV, "PR_[^_]*"))%>%
  mutate(group3 = str_extract(group, "[^_]*"))%>%
  mutate(group2 = coalesce(group2, group3))%>%
  mutate(percentage = (N_SITES-O.HOM.)/N_SITES) 


heteperce_mean <- heteperce2 %>%
  group_by(group2) %>% 
  dplyr::summarize(percent_mean = mean(percentage, na.rm=TRUE)) 

group_ordered <- with(heteperce_mean,  reorder(group2,  percent_mean,  mean))
heteperce2$group2 <- factor(heteperce2$group2, levels = levels(group_ordered))


heteperce_Plot <- ggplot(heteperce2, aes(x=group2, y=percentage, fill = group2)) + 
  geom_boxplot(show.legend = F, width = 0.7, outlier.shape = NA) + 
  geom_text(aes(label = formatC(after_stat(y), format = "f", digits  = 3), group = group2), size = 4,
            stat = 'summary', fun = mean,  nudge_y = 0.03 ) +
  stat_summary(fun.y=mean, geom="point", shape=18, size=3, color="blue", fill="blue") +
  scale_fill_manual(values=safe_colorblind_palette)  +
  geom_jitter(color="black", size=1, alpha=0.9) +
  xlab("Populations/Groups") +
  ylab("Proportion of heterozygosity sites") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = -45, v = -1, size = 9),
        panel.border = element_rect(colour = "black", fill=NA))
heteperce_Plot


#####################################################################
####################################################################
####################################################################
library(cowplot)
finalplot <- ggdraw() +
  draw_plot(ROH_plot, x = 0, y = 0.55, width = 1, height = 0.45) +
  draw_plot(LD_dec_plot , x = 0, y = 0, width = 0.333333, height = 0.5) +
  draw_plot(Roh_Froh_Plot , x = 0.333333, y = 0, width = 0.333333, height = 0.5) +
  draw_plot(heteperce_Plot, x = 0.66666, y = 0, width = 0.33333, height = 0.5) +
  
  draw_plot_label(label = c("A", "B", "C", "D"), size = 15,  fontface = "bold",
                  x = c(0, 0, 0.33333,0.66666), y = c(1, 0.5,  0.5, 0.5))

pdf("Fig6_ROH_n121v2.pdf", width = 15, height = 10)
finalplot
dev.off()



















#####################################################
#####################################################
#####################################################

ROH_total <- slidingRuns3 %>%
  mutate(id2 = str_extract(group, "[^_]*"))%>% 
  mutate(id3 = str_extract(group, "PR_[^_]*"))%>% 
  mutate(id = coalesce(id3, id2))%>% 
  group_by(id, category) %>% 
  dplyr::summarize(ROH_sum = mean(lengthBps_sum, na.rm=TRUE)) %>%
  as.data.frame()

group_ordered <- with(ROH_total,  reorder(id,  ROH_sum,  sum, decreasing = TRUE))
ROH_total$id <- factor(ROH_total$id, levels = levels(group_ordered))


ROH_plot2 <- ggplot(ROH_total, aes(fill=category, y=id, x=ROH_sum )) + 
  geom_bar(position="stack", stat="identity", width = 0.7) + 
  scale_fill_manual(values=c("0.25-1 Mb" = "steelblue",
                             "1-2 Mb" = "#E69F00",
                             ">2 Mb" = "#009E73")) +
  ylab("Populations/Groups") +
  xlab("Propotion of total genome") +
  theme_classic()+
  theme(legend.position =  c(0.8, 0.8),
        #legend.box.background = element_rect(linetype="dotted", colour = "black"),
        #legend.title=element_blank(),
        legend.key = element_rect(fill = "white", colour = "grey30"))
ROH_plot2$labels$fill <- "Length of ROH"
ROH_plot2
#########

slidingRuns4 <- slidingRuns2 %>%
  group_by(id) %>% 
  add_count(id) %>%
  group_by(id,n) %>% 
  dplyr::summarize(ROH_sum = sum(lengthBps, na.rm=TRUE)/2296245394) %>%
  mutate(id = gsub("AD1_", "", id))%>% 
  mutate(group2 = str_extract(id, "[^_]*"))%>% 
  mutate(group3 = str_extract(id, "PR_[^_]*"))%>% 
  mutate(group = coalesce(group3, group2))%>% 
  as.data.frame()

slidingRuns4 %>%
  group_by(group) %>% 
  dplyr::summarize(ROH_sum_mean = mean(ROH_sum, na.rm=TRUE), 
                   n_mean = mean(n, na.rm=TRUE))
  
cor(slidingRuns4$n, slidingRuns4$ROH_sum)



NR_ROH <- ggplot(slidingRuns4, aes(n, ROH_sum, col = group, shape = group)) + 
  geom_point(size=3) +
  xlab("Number of ROH") +
  ylab("Proportion of total genome") + 
  scale_color_manual(values= safe_colorblind_palette) +
  scale_shape_manual(values=safe_shape_palette) +
  theme_bw() +  
  theme(legend.position =  c(0.85, 0.25),
        legend.key = element_rect(fill = "white", colour = "grey30", linetype="dotted")  ) 

NR_ROH$labels$colour <- "Population/Group"
NR_ROH$labels$shape <- "Population/Group"

#####################################################################

finalplot <- ggdraw() +
  draw_plot(ROH_plot2, x = 0, y = 0, width = 0.5, height = 1) +
  draw_plot(NR_ROH , x = 0.5, y = 0, width = 0.5, height = 1) +
  draw_plot_label(label = c("a", "b"), size = 15,  fontface = "bold",
                  x = c(0, 0.5), y = c(1,1))

finalplot


pdf("FigS6_ROH_n121.pdf", width = 12, height = 6)
finalplot
dev.off()
