setwd("C:/Users/weixuan/Desktop/GDPR_samples/00_FigureTables/01_AD1AD2AD4_n145/")

library(mapdata)
library(usmap)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(sp)
library(ggspatial)
library(ggmap)
library(ggsn)
library(dplyr)
library(maps)
library(maptools)
library(ggplot2)
library(ggrepel)
library(grid)
library(cowplot)

##################################################################
############## Caribbean region ##################################
##################################################################

register_google("AIzaSyAFMCC4ztGaZzBKvTMo4ByO5pBZVPhR5tQ") #for the ggmap package

gpsmk <- read.csv("CottonSamplingMap3.csv")

gpsmk2 <- gpsmk %>%
  filter(Source == "Field") %>%
  filter(Population == "Puerto Rico" | Population == "Mound Key" | Population == "Guadalupe")  %>%
  filter(Group != "Landrace")

#########################

mapdot <- gpsmk %>%
  filter(Source == "Field") %>%
  filter(Region == "Puerto Rico")  %>%
  mutate(maploc = gsub("AD1_", "", Pop_Site)) %>%
  add_count(maploc, name = 'id_occurence')  %>%
  mutate(maploc2 = paste0( maploc, " (n = ", id_occurence, ")")) %>%
  select(maploc2, lat, long) %>% 
  distinct

##################################################################
##################################################################
##################################################################
# Get the map from Google
map1 <- get_map(c(left = -98, bottom = 4, right = -55, top = 31), 
               zoom = 4, 
               source = "google", 
               maptype = "satellite")

# Plot the map
bigmap <- ggmap(map1) +
  xlim(lon_range) +
  ylim(lat_range) +
  coord_sf(xlim = c(-98, -55), ylim = c(4, 31), expand = FALSE, crs = "WGS84")+
  xlab("Longitude") + ylab("Latitude") +
  geom_point(data = gpsmk2, mapping = aes(x = long, y = lat-2,  shape = Source), color = "yellow", size =4, show.legend = F) +
  annotate(geom = "text", x = -84, y = 25, label = "Mound Key", fontface = "bold", color = "white", size = 4) +
  annotate(geom = "text", x = -65, y = 17, label = "Puerto Rico", fontface = "bold", color = "white", size = 4) +
  annotate(geom = "text", x = -60, y = 15, label = "Guadeloupe", fontface = "bold", color = "white", size = 4) +
  scale_x_continuous(labels = scales::label_number(accuracy = 0.01)) +
  scale_y_continuous(labels = scales::label_number(accuracy = 0.01))+
  scalebar( x.min = -94, x.max = -90, y.min = 6, y.max = 8,
            box.fill = c("yellow", "white"), st.color = "white",
            dist = 200, dist_unit = "km", st.size =3, height =0.3, st.dist = 0.2,
            transform = TRUE, model = "WGS84") +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.0, "in"), pad_y = unit(0.15, "in"),
                         style = north_arrow_fancy_orienteering(text_col = 'white',
                                                                line_col = 'white')) +
  #theme_void() +
  theme(panel.background = element_rect(fill = "white"),
        plot.background = element_blank(),   # Remove overall background
        axis.text = element_text(family = 'serif'),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.1),
        panel.grid.major = element_line(color = "grey90", linetype = "dashed",  size = 0.2))


##################################################################
##################################################################
##################################################################

lon_range <- c(-61.192, -61.187)  # Longitude (xlim)
lat_range <- c(16.2421, 16.2450)     # Latitude (ylim)

# Get the map from Google
map2 <- get_map(location = c(lon = -61.180, lat = 16.245), 
               zoom = 15, 
               source = "google", 
               maptype = "satellite")

GD_samplingplot <- ggmap(map2) +
  coord_sf(xlim = c(-61.193, -61.169), ylim = c(16.2415, 16.257), expand = FALSE, crs = "WGS84") +
  #  xlab("Longitude") + ylab("Latitude") +
  annotate(geom = "text", x = -61.180, y = 16.252, label = "GD (n = 21)", fontface = "bold", color = "white", size = 3) +
  annotate(geom = "text", x = -61.185, y = 16.245, label = "GD2 (n = 4)", fontface = "bold", color = "white", size = 3) +
  scalebar( x.min = -61.191, x.max = -61.186, y.min = 16.2430, y.max = 16.246,
            box.fill = c("yellow", "white"), st.color = "white",
            dist = 200, dist_unit = "m", st.size =3, height =0.14, st.dist = .12,
            transform = TRUE, model = "WGS84") +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.0, "in"), pad_y = unit(0.2, "in"),
                         style = north_arrow_fancy_orienteering(text_col = 'white',
                                                                line_col = 'white')) +
  theme_void() +
  theme(panel.background = element_blank(),
        plot.background = element_blank(),   # Remove overall background
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = c(0.86, 0.7),
        legend.background = element_rect(fill= "transparent", size=.5, linetype="dotted"))

GD_samplingplot

##################################################################
##################################################################
##################################################################

lon_range <- c(-67.25, -66.8)  # Longitude (xlim)
lat_range <- c(17.8, 18.1)     # Latitude (ylim)

# Get the map from Google
map3 <- get_map(location = c(lon = mean(lon_range), lat = mean(lat_range)), 
               zoom = 11, 
               source = "google", 
               maptype = "satellite")

ggmap(map)

PR_samplingplot <- ggmap(map3) +
  coord_sf(xlim = c(-67.23, -66.81) , ylim = c(17.84, 18.1), expand = FALSE, crs = "WGS84") +
  geom_point(data = mapdot3, mapping = aes(x = long, y = lat), color = "yellow", size = 3, show.legend = F) +
  geom_text_repel(data = mapdot, aes(x = long, y = lat, label=maploc2), 
                  min.segment.length = unit(0, 'lines'),  color = "white", 
                  max.overlaps=Inf, size = 3, show.legend = F ) +
  scalebar( x.min = -67.18, x.max = -67.1, y.min = 17.86, y.max = 17.89,
            box.fill = c("yellow", "white"), st.color = "white",
            dist = 4, dist_unit = "km", st.size =3, height =0.3, st.dist = 0.2,
            transform = TRUE, model = "WGS84") +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.0, "in"), pad_y = unit(0.2, "in"),
                         style = north_arrow_fancy_orienteering(text_col = 'white',
                                                                line_col = 'white')) +
  theme_void() +
  theme(panel.background = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.1))

  


##################################################################
##################################################################
##################################################################
finalplot <- ggdraw() +
  draw_plot(bigmap , x = 0, y = 0.05, width = 0.6, height = 1) +
  draw_plot(PR_samplingplot , x = 0.6, y = 0.5, width = 0.397, height = 0.5) +
  draw_plot(GD_samplingplot, x = 0.6, y = 0, width = 0.397, height = 0.7) 

finalplot
#  draw_plot(bigmap, x = 0, y = 0.5, width = 0.6, height = 0.5) +
#  draw_plot(PR_samplingplot , x = 0.6, y = 0.75, width = 0.4, height = 0.25) +
#  draw_plot(GD_samplingplot, x = 0.6, y = 0.5, width = 0.4, height = 0.25) +
  
#  draw_plot(pca_plot, x = 0, y = 0, width = 0.6, height = 0.5) +
#  draw_plot(GDnjtree, x = 0.6, y = 0, width = 0.4, height = 0.5) +
#  draw_plot_label(label = c("b)", "c)"), size = 20,  fontface = "bold",
#                  x = c(0, 0.5), y = c(1, 1))

pdf("GD_PR_map2.pdf", width = 15, height = 10)
finalplot
dev.off()
