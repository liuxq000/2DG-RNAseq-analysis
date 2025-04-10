library(fmsb)

df <- read_xlsx("radar_data.xlsx",sheet = 1)
rownames(df) <- df$...1
df <- df[,-1]

pdf("radar_gene.pdf",width=10)
par(mfrow=c(2,3))
radarchartcirc(df[c(1,2,3),],
               cglty = 1,       # Grid line type
               cglcol = "lightgray", # Grid line color
               pcol = "#299D8F", # Color for each line
               pfcol = alpha("#299D8F", 0.3),    # Color for each fill polygon
               plwd = 2,        # Width for each line
               plty = 1,
               axistype = 1,
               title = "Sham",
               axislabcol = "gray")

radarchartcirc(df[c(1,2,4),],
               cglty = 1,       # Grid line type
               cglcol = "lightgray", # Grid line color
               pcol = "#55B7E6", # Color for each line
               pfcol = alpha("#55B7E6", 0.3),    # Color for each fill polygon
               plwd = 2,        # Width for each line
               plty = 1,
               axistype = 1,
               title = "CLP",
               axislabcol = "gray")

radarchartcirc(df[c(1,2,5),],
               cglty = 1,       # Grid line type
               cglcol = "lightgray", # Grid line color
               pcol = "#B395BD", # Color for each line
               pfcol = alpha("#B395BD", 0.3),    # Color for each fill polygon
               plwd = 2,        # Width for each line
               plty = 1,
               axistype = 1,
               title = "CLP2DG",
               axislabcol = "gray")

radarchartcirc(df[c(1,2,6),],
               cglty = 1,       # Grid line type
               cglcol = "lightgray", # Grid line color
               pcol = "#C59D94", # Color for each line
               pfcol = alpha("#C59D94", 0.3),    # Color for each fill polygon
               plwd = 2,        # Width for each line
               plty = 1,
               axistype = 1,
               title = "ShamLPS",
               axislabcol = "gray")

radarchartcirc(df[c(1,2,7),],
               cglty = 1,       # Grid line type
               cglcol = "lightgray", # Grid line color
               pcol = "#EA8379", # Color for each line
               pfcol = alpha("#EA8379", 0.3),    # Color for each fill polygon
               plwd = 2,        # Width for each line
               plty = 1,
               axistype = 1,
               title = "CLPLPS",
               axislabcol = "gray")

radarchartcirc(df[c(1,2,8),],
               cglty = 1,       # Grid line type
               cglcol = "lightgray", # Grid line color
               pcol = "#F09739", # Color for each line
               pfcol = alpha("#F09739", 0.3),    # Color for each fill polygon
               plwd = 2,        # Width for each line
               plty = 1,
               axistype = 1,
               title = "CLP2DGLPS",
               axislabcol = "gray")
dev.off()

##############################################################################################################

library(ggplot2)
library(dplyr)

set.seed(123) 

gene_data <- read_xlsx("gene_anno.xlsx",sheet = 1)
plot_data <- gene_data %>%
  mutate(
    angle = seq(0, 2*pi*(1 - 1/nrow(.)), length.out  = nrow(.)),
    angle = angle + pi/2,  
    radius = 5,
    x = radius * cos(angle),
    y = radius * sin(angle)
  )
class_colors <- c("Cytokines"="#DC0000","ROS"="#7E6148","NETs"="#3C5488","Chemokines"="#a2d2a4")

pdf("gene_point.pdf")
ggplot(plot_data, aes(x = x, y = y)) +
  geom_point(
    aes(color = Class),   
    size = 6,              
    shape = 19             
  ) +
  scale_color_manual(values = class_colors) + 
  coord_fixed(ratio = 1) +                     
  theme_void() +                               
  theme(
    legend.position  = "right",                
    plot.margin  = unit(c(1,1,1,1), "cm")      
  )
dev.off()

