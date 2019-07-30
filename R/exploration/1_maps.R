# Creating global map to highlight differences in approach

# install.packages("ggmap")

library(ggmap)

mp <- NULL
mapWorld <- borders("world", colour="gray70", fill="gray70") # create a layer of borders
mp <- ggplot() + mapWorld

mp + 
  theme_classic(base_size = 16) +
  coord_cartesian(expand = 0) + 
  theme(aspect.ratio = 1/2)
