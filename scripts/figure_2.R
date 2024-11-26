# This script generates figure 2 from PAPER TITLE
# It first generates the following sub-pannels:
# a) global map of study locations
# b) bar plot of samples across tissues and sex
# c) histogram of the age distribution
# d) bar plot of tissues in which the different plastic polymers were found
# e) boxplots of the concentrations (particles/g) of microplastics in the different tissues
# f) boxplots of the concentrations (ug/g) of microplastics in the different tissues
# it then combines the sub-pannels into a single figure

# Last updated: November 14 2024

#Import the packages
library(ggplot2)
library(tidyr)
library(dplyr)
library(terra)
library(sf)
library(tidyterra)
library(gridExtra)
library(viridis)

# Load data
data <- read.csv("data/MPs_human_tissues.csv", header = TRUE)
data[which(data$Tissue == "Olfactory_bulb"),"Tissue"] <- "Olfactory bulb"
data$Tissue <- as.factor(data$Tissue)
data$Sex <- as.factor(data$Sex)

#---------------------------------------------------------------------
# Figure 2A - study location map with HFI
#---------------------------------------------------------------------

# Import the human footprint index raster from: https://www.frontiersin.org/articles/10.3389/frsen.2023.1130896/full
#Data are available here: https://source.coop/repositories/vizzuality/hfp-100/description
HFI <- rast("~/Dropbox/UBC/Projects/microplastics_brazil/data/environmental_data/hfp_2021_100m_v1-2_cog.tif")

#Reproject and process the HFI data
world_sf <- st_as_sf(rworldmap::getMap(resolution = "low"))
world_sf <- subset(world_sf, continent != "Antarctica")
world <- vect(lwgeom::st_transform_proj(world_sf, crs = crs(HFI)))
HFI <- mask(HFI, world) #Note: this step is slow
plot(HFI)


# Count the number of patients in each city
city_data <- aggregate(Species ~ City, data = data, FUN = "length")

# Manually assigned latitude and longitude for each city
city_data <- city_data %>%
  mutate(
    lat = case_when(
      City == "Beijing" ~ 39.9042,
      City == "Foggia" ~ 41.4622,
      City == "Kota Bharu" ~ 4.2105,
      City == "Miami" ~ 25.7617,
      City == "Shanghai" ~ 31.2304,
      City == "Sao Paulo" ~ -23.5505,
      City == "Guangdong" ~ 23.3790,
      City == "Albuquerque" ~ 35.0844 
    ),
    lon = case_when(
      City == "Beijing" ~ 116.4074,
      City == "Foggia" ~ 15.5446,
      City == "Kota Bharu" ~ 101.9758,
      City == "Miami" ~ -80.1918,
      City == "Shanghai" ~ 121.4737,
      City == "Sao Paulo" ~ -46.6333,
      City == "Guangdong" ~ 113.7633,
      City == "Albuquerque" ~ -106.6504
    )
  )

#convert to sf and reproject
locations <- st_as_sf(city_data,
                      coords = c("lon", "lat"),
                      crs = 4326) %>% 
  st_transform(dst = crs(HFI))


a <- 
  ggplot() +
  ggtitle("a") +
  geom_spatraster(data = HFI, maxcell = 5e+06,
                  alpha = 0.7) +
  scale_fill_gradient(low = "#eeeeee",
                      high = "black",
                      na.value = NA,
                      guide = FALSE) +
  #Add locations of study sites
  #geom_sf(data = world, size = 0.2, fill = "transparent", linewidth = 0.1, col = "black") +
  geom_sf(data = locations, aes(size = city_data$Species), fill = "#fca311", shape = 21) +
  scale_size_continuous(name = "Number of patients") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.position = "inside",
        legend.position.inside = c(0.15,0.3),
        legend.title = element_text(size=6, family = "sans", face = "bold"),
        legend.text = element_text(size=6, family = "sans", face = "bold"),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(0.2, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        plot.title = element_text(hjust = .1, vjust = -2.5, size = 10, family = "sans", face = "bold"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks = element_blank(),
        strip.background=element_blank(),
        plot.margin = unit(c(0.4,0.6,0.4,-0.2), "cm")) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) 


#---------------------------------------------------------------------
# Figure 2B - tissue type bar chart 
#---------------------------------------------------------------------

#vector of names to order the figure elements
ORDER <- names(sort(table(data$Tissue)))

b <-
  ggplot(data = data) +
  ggtitle("b") +
  geom_bar(aes(Tissue, fill =Sex),  alpha = 0.9, size = 0.1, width = 0.8) +
  scale_fill_manual(values = c("#184E77", "#de425b")) +
  coord_flip() +
  scale_x_discrete(limits = ORDER) +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(linewidth = 0.2),
        axis.line.x = element_line(linewidth = 0.2),
        axis.ticks.x = element_line(linewidth = 0.2),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size=6, family = "sans", face = "bold"),
        axis.text.y = element_text(size=5, family = "sans", face = "bold", color = "black"),
        axis.text.x  = element_text(size=5, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 10, family = "sans", face = "bold"),
        legend.position = "inside",
        legend.position.inside = c(0.7,0.2),
        legend.title = element_blank(),
        legend.text = element_text(size=6, family = "sans", face = "bold"),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(0.2, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(1.2,0.1,1.3,-0.5), "cm")) +
  guides(fill = guide_legend(byrow = TRUE)) +
  ylab("Number of patients") +
  scale_y_continuous(expand = expansion(mult = c(0,.1), add = 0))


#---------------------------------------------------------------------
# Figure 2c - histogram of ages 
#---------------------------------------------------------------------

c <- 
  ggplot(data, aes(x = Age)) +
  geom_histogram(position = "identity", fill = "#fca311") +
  ggtitle("c") +
  xlab("Age (years)") +
  ylab ("Number of patients") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(linewidth = 0.2),
        axis.line.x = element_line(linewidth = 0.2),
        axis.ticks.x = element_line(linewidth = 0.2),
        axis.ticks.y = element_line(linewidth = 0.2),
        axis.title.y = element_text(size=6, family = "sans", face = "bold"),
        axis.title.x = element_text(size=6, family = "sans", face = "bold"),
        axis.text.y = element_text(size=5, family = "sans"),
        axis.text.x  = element_text(size=5, family = "sans", colour = "black"),
        strip.text.x = element_text(size=9, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 10, family = "sans", face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size=6, family = "sans", face = "bold"),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(0.2, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm")) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) 


#---------------------------------------------------------------------
# Figure 2d - histogram of ages 
#---------------------------------------------------------------------

#vector of names to order the figure elements
ORDER <- rev(names(sort(table(data$Most_prevalent_MP))))

d <-  
  ggplot(data, aes(fill=Tissue, x=Most_prevalent_MP)) + 
  geom_bar(position="stack", stat="count") +
  scale_fill_viridis(discrete = TRUE) +
  ggtitle("d") +
  xlab("Most prevalent polymer") +
  ylab ("Number of patients") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(linewidth = 0.2),
        axis.line.x = element_line(linewidth = 0.2),
        axis.ticks.x = element_line(linewidth = 0.2),
        axis.ticks.y = element_line(linewidth = 0.2),
        axis.title.y = element_text(size=6, family = "sans", face = "bold"),
        axis.title.x = element_text(size=6, family = "sans", face = "bold"),
        axis.text.y = element_text(size=5, family = "sans"),
        axis.text.x  = element_text(size=5, family = "sans", colour = "black"),
        strip.text.x = element_text(size=9, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 10, family = "sans", face = "bold"),
        legend.position = "inside",
        legend.position.inside = c(0.7,0.7),
        legend.title = element_blank(),
        legend.text = element_text(size=6, family = "sans", face = "bold"),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(0.2, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm")) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0), limits = ORDER) 

#---------------------------------------------------------------------
# Figure 2e - boxplot of concentration by tissue type
#---------------------------------------------------------------------

#vector of names to order the figure elements
ORDER <- aggregate(MPs_concentration_particles_g ~ Tissue, data = data, FUN = "mean")
ORDER <- ORDER[rev(order(ORDER$MPs_concentration_particles_g)),"Tissue"]


e <- 
  ggplot(data[!is.na(data$MPs_concentration_particles_g),],
         aes(x = Tissue, y = MPs_concentration_particles_g)) +
  geom_boxplot(aes(fill = Tissue),
               outlier.shape = NA,
               color = "black",
               alpha = 0.7,
               linewidth = 0.2) +
  geom_jitter(aes(color = Tissue), width = 0.2, size = 0.4, alpha = 0.7, shape = 16) +
  stat_summary(fun = mean,
               geom = "point",
               shape = 23,
               size = 0.8,
               fill = "black",
               color = "black") +
  
  ggtitle("e") +
  xlab("Tissue") +
  ylab ("MP concentration (particles/g)") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(linewidth = 0.2),
        axis.line.x = element_line(linewidth = 0.2),
        axis.ticks.x = element_line(linewidth = 0.2),
        axis.ticks.y = element_line(linewidth = 0.2),
        axis.title.y = element_text(size=6, family = "sans", face = "bold"),
        axis.title.x = element_text(size=6, family = "sans", face = "bold"),
        axis.text.y = element_text(size=5, family = "sans"),
        axis.text.x  = element_text(size=5, family = "sans", colour = "black"),
        strip.text.x = element_text(size=9, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 10, family = "sans", face = "bold"),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size=6, family = "sans", face = "bold"),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(0.2, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm"))  +
  scale_fill_viridis(discrete = TRUE) +
  scale_color_viridis(discrete = TRUE) +
  scale_y_log10(limits = c(0.1,100000), breaks = c(0.1,1,10,100,1000,10000,100000)) +
  scale_x_discrete(guide = guide_axis(angle = 90), limits = ORDER) +
  annotation_logticks(sides= "l",
                      short = unit(0.04, "cm"),
                      mid = unit(0.04, "cm"),
                      long = unit(0.1, "cm"),
                      linewidth = 0.2)



#---------------------------------------------------------------------
# Figure 2f - boxplot of concentration by tissue type
#---------------------------------------------------------------------

#vector of names to order the figure elements
ORDER <- aggregate(MPs_concentration_weight ~ Tissue, data = data, FUN = "mean")
ORDER <- ORDER[rev(order(ORDER$MPs_concentration_weight)),"Tissue"]

f <- 
  ggplot(data[!is.na(data$MPs_concentration_weight),],
         aes(x = Tissue, y = MPs_concentration_weight)) +
  geom_boxplot(aes(fill = Tissue),
               outlier.shape = NA,
               color = "black",
               alpha = 0.7,
               linewidth = 0.2) +
  geom_jitter(aes(color = Tissue), width = 0.2, size = 0.4, alpha = 0.7, shape = 16) +
  stat_summary(fun = mean,
               geom = "point",
               shape = 23,
               size = 0.8,
               fill = "black",
               color = "black") +
  ggtitle("f") +
  xlab("Tissue") +
  ylab ("MP concentration (ug/g)") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(linewidth = 0.2),
        axis.line.x = element_line(linewidth = 0.2),
        axis.ticks.x = element_line(linewidth = 0.2),
        axis.ticks.y = element_line(linewidth = 0.2),
        axis.title.y = element_text(size=6, family = "sans", face = "bold"),
        axis.title.x = element_text(size=6, family = "sans", face = "bold"),
        axis.text.y = element_text(size=5, family = "sans"),
        axis.text.x  = element_text(size=5, family = "sans", colour = "black"),
        strip.text.x = element_text(size=9, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 10, family = "sans", face = "bold"),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size=6, family = "sans", face = "bold"),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(0.2, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm"))  +
  
  # Color Scales for Fill and Color
  scale_fill_viridis(discrete = TRUE) +
  scale_color_viridis(discrete = TRUE) +
  scale_y_log10(limits = c(0.1,100000), breaks = c(0.1,1,10,100,1000,10000,100000)) +
  scale_x_discrete(guide = guide_axis(angle = 90), limits = ORDER) +
  annotation_logticks(sides= "l",
                      short = unit(0.04, "cm"),
                      mid = unit(0.04, "cm"),
                      long = unit(0.1, "cm"),
                      linewidth = 0.2)



#---------------------------------------------------------------------
# Combine and save
#---------------------------------------------------------------------


TOP <-
  grid.arrange(a,b,
               ncol=2,
               nrow=1,
               widths=c(6,3))


MID <-
  grid.arrange(c,d,
               ncol=2,
               nrow=1) 

BOT <-
  grid.arrange(e,f,
               ncol=2,
               nrow=1,
               widths=c(6,3)) 

fig2 <-
  grid.arrange(TOP,
               MID,
               BOT,
               ncol=1,
               nrow=3,
               heights=c(3,1.5,2))


#Save the figures
ggsave(fig2,
       width = 6.86, height = 7, units = "in",
       dpi = 600,
       bg = "transparent",
       file="figures/Figure_2.png")
