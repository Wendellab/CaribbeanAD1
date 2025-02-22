setwd("C:/Users/weixuan/Desktop/GDPR_samples/00_FigureTables/03_TreemixDsuite_n128/Treemixoutput/") # of course this needs to be adjusted

library(RColorBrewer)
library(R.utils)
source("C:/Users/weixuan/Desktop/GDPR_samples/00_FigureTables/03_TreemixDsuite_n128/plotting_funcs.R") # here you need to add the path

prefix="TreeMix.7"

pdf("Treemix.pdf", height = 20, width = 20)
par(mfrow=c(2,3))
for(edge in 1:5){
  plot_tree(cex=1,font =1, scale = F, ybar = 0.3, plus = 0.2, paste0(prefix,".",edge))
  title(paste("Three edges"))
}
dev.off()


library(ggplotify)
p <- as.ggplot(function() plot_tree(cex=1, font =1, scale = F, mbar = T,  disp = 0.01, ybar = 0.1, 
                                    plus = 0.05, arrow = 0.1, lwd = 3, "TreeMix.7.1") )

p

library(cowplot)
finalplot <- ggdraw() +
  draw_plot(p, x = 0, y = 0, width = .5, height = 1)  +
  draw_plot_label(label = c("a)", "b)"), size = 15,  fontface = "bold",
                  x = c(0, 0.5), y = c(1, 1))

setwd("C:/Users/weixuan/Desktop/GDPR_samples/00_FigureTables/03_TreemixDsuite_n128/") 

pdf("Fig3_GDPRMK_Geneflow_n128.pdf", width = 10, height = 5)
finalplot
dev.off()

