# This script generates figure 2 from PAPER TITLE
# It first generates the following sub-pannels:
# a) Relationship with age for the spectroscopy data
# b) Relationship with age for the Pyrolysis-GC data
# c) Sex-based differences in MPs within each tissue for the spectroscopy data
# d) Sex-based differences in MPs within each tissue for the Pyrolysis-GC data
# it then combines the 4 sub-panels into a single figure

# Last updated: November 14 2024

#Import the packages
library(ggplot2)
library(gridExtra)
library(mgcv)
library(viridis)

# Load data
data <- read.csv("data/MPs_human_tissues.csv", header = TRUE)
data[which(data$Tissue == "Olfactory_bulb"),"Tissue"] <- "Olfactory bulb"
data$Tissue <- as.factor(data$Tissue)
data$Sex <- as.factor(data$Sex)
data$Article_number <- as.factor(data$Article_number)
data$Method_analysis <- as.factor(data$Method_analysis)


#---------------------------------------------------------------------
# Figure 4a - Relationship with age for the spectroscopy data
#---------------------------------------------------------------------

#Relationship with age for the spectroscopy data
FIT <- gam(MPs_concentration_particles_g ~  Age + s(Age, Tissue, bs = "re") +
             s(Article_number, bs = "re") +
             s(Tissue, bs = "re") +
             s(Method_analysis, bs = "re"),
           family = tw(),
           method = "REML",
           data = data[!is.na(data$MPs_concentration_particles_g),])


TISSUES <- unique(data[!is.na(data$MPs_concentration_particles_g),"Tissue"])

#Empty list for storing results
res <- list()

#Loop over all of the tissue types
for(i in 1:length(TISSUES)){
  
  #New dataset for predicting
  NEW_DATA <- data.frame(Age = 1:100,
                         Tissue = TISSUES[i],
                         Article_number = factor("new"),
                         Method_analysis = factor("new"))
  
  #Make predictions for that tissue type
  NEW_DATA$mp_hat <- predict(FIT, newdata = NEW_DATA, type = "response", se.fit = FALSE, exclude = c("s(Article_number)","s(Method_analysis)"))
  
  res[[i]] <-  NEW_DATA
}

res <- do.call(rbind, res)

#New dataset for predicting
NEW_DATA <- data.frame(Age = 1:100,
                       Tissue = factor("new"),
                       Article_number = factor("new"),
                       Method_analysis = factor("new"))

#Make predictions for that tissue type
mp_hat <- predict(FIT, newdata = NEW_DATA, type = "link", se.fit = TRUE, exclude = c("s(Article_number)","s(Method_analysis)"))
NEW_DATA$mp_hat <- exp(mp_hat$fit)
NEW_DATA$mp_hat_min <- exp(mp_hat$fit - mp_hat$se.fit*1.96)
NEW_DATA$mp_hat_max <- exp(mp_hat$fit + mp_hat$se.fit*1.96)

a <- 
ggplot() +
  geom_point(data = data[!is.na(data$MPs_concentration_particles_g),], aes(x = Age, y = MPs_concentration_particles_g, col = Tissue), size = 0.5) +
  geom_line(data = res, aes(x = Age, y = mp_hat, col = Tissue), linetype = "dotted", linewidth = 0.2) +
  geom_line(data = NEW_DATA, aes(x = Age, y = mp_hat), linewidth = 0.8, col = "black") +
  geom_ribbon(data = NEW_DATA, aes(x = Age, ymin = mp_hat_min, ymax = mp_hat_max), alpha = 0.2, fill = "grey70") +
  ggtitle("a") +
  xlab("Age (years)") +
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
        legend.title = element_blank(),
        legend.text = element_text(size=4, family = "sans", face = "bold"),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(0.2, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.position = "right",
        legend.position.inside = c(0.1,0.85),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm")) +
  scale_color_viridis(discrete = TRUE) +
  scale_y_log10(limits = c(0.1,100000), breaks = c(0.1,1,10,100,1000,10000,100000)) +
  #scale_x_continuous(expand = c(0.01,0)) +
  annotation_logticks(sides= "l",
                      short = unit(0.04, "cm"),
                      mid = unit(0.04, "cm"),
                      long = unit(0.1, "cm"),
                      linewidth = 0.2)



#---------------------------------------------------------------------
# Figure 4b - Relationship with age for males in the spectroscopy data
#---------------------------------------------------------------------

mps_males <- data[!is.na(data$MPs_concentration_particles_g),]
mps_males <- mps_males[mps_males$Sex == "Male",]

#Relationship with age for the spectroscopy data
FIT <- gam(MPs_concentration_particles_g ~  Age + s(Age, bs = "re") +
             s(Article_number, bs = "re") +
             s(Tissue, bs = "re") +
             s(Method_analysis, bs = "re"),
           family = tw(),
           method = "REML",
           data = mps_males)

#New dataset for predicting
NEW_DATA <- data.frame(Age = 1:100,
                       Sex = factor("Male"),
                       Tissue = factor("new"),
                       Article_number = factor("new"),
                       Method_analysis = factor("new"))

#Make predictions for that tissue type
mp_hat <- predict(FIT, newdata = NEW_DATA, type = "link", se.fit = TRUE, exclude = c("s(Article_number)","s(Method_analysis)","s(Tissue)"))
NEW_DATA$mp_hat <- exp(mp_hat$fit)
NEW_DATA$mp_hat_min <- exp(mp_hat$fit - mp_hat$se.fit*1.96)
NEW_DATA$mp_hat_max <- exp(mp_hat$fit + mp_hat$se.fit*1.96)

b <- 
  ggplot() +
  geom_point(data = mps_males, aes(x = Age, y = MPs_concentration_particles_g), col = "#de425b", size = 0.3) +
  geom_line(data = NEW_DATA, aes(x = Age, y = mp_hat), linewidth = 0.6, col = "#de425b") +
  geom_ribbon(data = NEW_DATA, aes(x = Age, ymin = mp_hat_min, ymax = mp_hat_max), alpha = 0.2, fill = "#de425b") +
  ggtitle("b") +
  xlab("Age (years)") +
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
        legend.title = element_blank(),
        legend.text = element_text(size=4, family = "sans", face = "bold"),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(0.2, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.position = "right",
        legend.position.inside = c(0.1,0.85),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm")) +
  scale_color_viridis(discrete = TRUE) +
  scale_y_log10(limits = c(0.1,100000), breaks = c(0.1,1,10,100,1000,10000,100000)) +
  #scale_x_continuous(expand = c(0.01,0)) +
  annotation_logticks(sides= "l",
                      short = unit(0.04, "cm"),
                      mid = unit(0.04, "cm"),
                      long = unit(0.1, "cm"),
                      linewidth = 0.2)



#---------------------------------------------------------------------
# Figure 4c - Relationship with age for females in the spectroscopy data
#---------------------------------------------------------------------

mps_females <- data[!is.na(data$MPs_concentration_particles_g),]
mps_females <- mps_females[mps_females$Sex == "Female",]

#Relationship with age for the spectroscopy data
FIT <- gam(MPs_concentration_particles_g ~  Age + s(Age, bs = "re") +
             s(Article_number, bs = "re") +
             s(Tissue, bs = "re") +
             s(Method_analysis, bs = "re"),
           family = tw(),
           method = "REML",
           data = mps_females)

#New dataset for predicting
NEW_DATA <- data.frame(Age = 1:100,
                       Sex = factor("Female"),
                       Tissue = factor("new"),
                       Article_number = factor("new"),
                       Method_analysis = factor("new"))

#Make predictions for that tissue type
mp_hat <- predict(FIT, newdata = NEW_DATA, type = "link", se.fit = TRUE, exclude = c("s(Article_number)","s(Method_analysis)","s(Tissue)"))
NEW_DATA$mp_hat <- exp(mp_hat$fit)
NEW_DATA$mp_hat_min <- exp(mp_hat$fit - mp_hat$se.fit*1.96)
NEW_DATA$mp_hat_max <- exp(mp_hat$fit + mp_hat$se.fit*1.96)

c <- 
  ggplot() +
  geom_point(data = mps_females, aes(x = Age, y = MPs_concentration_particles_g), col = "#184E77", size = 0.3) +
  geom_line(data = NEW_DATA, aes(x = Age, y = mp_hat), linewidth = 0.6, col = "#184E77") +
  geom_ribbon(data = NEW_DATA, aes(x = Age, ymin = mp_hat_min, ymax = mp_hat_max), alpha = 0.2, fill = "#184E77") +
  ggtitle("c") +
  xlab("Age (years)") +
  ylab ("MP concentration (particles/g)") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(linewidth = 0.2),
        axis.line.x = element_line(linewidth = 0.2),
        axis.ticks.x = element_line(linewidth = 0.2),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size=6, family = "sans", face = "bold"),
        axis.text.y = element_blank(),
        axis.text.x  = element_text(size=5, family = "sans", colour = "black"),
        strip.text.x = element_text(size=9, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 10, family = "sans", face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size=4, family = "sans", face = "bold"),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(0.2, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.position = "right",
        legend.position.inside = c(0.1,0.85),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm")) +
  scale_color_viridis(discrete = TRUE) +
  scale_y_log10(limits = c(0.1,100000), breaks = c(0.1,1,10,100,1000,10000,100000)) #+
#scale_x_continuous(expand = c(0.01,0)) +
# annotation_logticks(sides= "l",
#                     short = unit(0.04, "cm"),
#                     mid = unit(0.04, "cm"),
#                     long = unit(0.1, "cm"),
#                     linewidth = 0.2)


#---------------------------------------------------------------------
# Figure 4d - Relationship with age for the py-GC data
#---------------------------------------------------------------------

#Relationship with age for the spectroscopy data
FIT <- gam(MPs_concentration_weight ~  Age + s(Age, Tissue, bs = "re") +
             s(Article_number, bs = "re") +
             s(Tissue, bs = "re"),
           family = tw(),
           method = "REML",
           data = data[!is.na(data$MPs_concentration_weight),])


TISSUES <- unique(data[!is.na(data$MPs_concentration_weight),"Tissue"])

#Empty list for storing results
res <- list()

#Loop over all of the tissue types
for(i in 1:length(TISSUES)){
  
  #New dataset for predicting
  NEW_DATA <- data.frame(Age = 1:100,
                         Tissue = TISSUES[i],
                         Article_number = factor("new"))
  
  #Make predictions for that tissue type
  NEW_DATA$mp_hat <- predict(FIT, newdata = NEW_DATA, type = "response", se.fit = FALSE, exclude = c("s(Article_number)"))
  
  res[[i]] <-  NEW_DATA
}

res <- do.call(rbind, res)

#New dataset for predicting
NEW_DATA <- data.frame(Age = 1:100,
                       Tissue = factor("new"),
                       Article_number = factor("new"))

#Make predictions for that tissue type
mp_hat <- predict(FIT, newdata = NEW_DATA, type = "link", se.fit = TRUE, exclude = c("s(Article_number)"))
NEW_DATA$mp_hat <- exp(mp_hat$fit)
NEW_DATA$mp_hat_min <- exp(mp_hat$fit - mp_hat$se.fit*1.96)
NEW_DATA$mp_hat_max <- exp(mp_hat$fit + mp_hat$se.fit*1.96)

d <- 
  ggplot() +
  geom_point(data = data[!is.na(data$MPs_concentration_weight),], aes(x = Age, y = MPs_concentration_weight, col = Tissue), size = 0.5) +
  geom_line(data = res, aes(x = Age, y = mp_hat, col = Tissue), linetype = "dotted", linewidth = 0.2) +
  geom_line(data = NEW_DATA, aes(x = Age, y = mp_hat), linewidth = 0.8, col = "black") +
  geom_ribbon(data = NEW_DATA, aes(x = Age, ymin = mp_hat_min, ymax = mp_hat_max), alpha = 0.2, fill = "grey70") +
  ggtitle("d") +
  xlab("Age (years)") +
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
        legend.title = element_blank(),
        legend.text = element_text(size=4, family = "sans", face = "bold"),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(0.2, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.position = "right",
        legend.position.inside = c(0.1,0.15),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm")) +
  scale_color_viridis(discrete = TRUE) +
  scale_y_log10(limits = c(0.1,100000), breaks = c(0.1,1,10,100,1000,10000,100000)) +
  #scale_x_continuous(expand = c(0.01,0))  +
  annotation_logticks(sides= "l",
                      short = unit(0.04, "cm"),
                      mid = unit(0.04, "cm"),
                      long = unit(0.1, "cm"),
                      linewidth = 0.2)



#---------------------------------------------------------------------
# Figure 4e - Relationship with age for males in the py-GC data
#---------------------------------------------------------------------

mps_males <- data[!is.na(data$MPs_concentration_weight),]
mps_males <- mps_males[mps_males$Sex == "Male",]

#Relationship with age for the spectroscopy data
FIT <- gam(MPs_concentration_weight ~  Age + s(Age, bs = "re") +
             s(Article_number, bs = "re") +
             s(Tissue, bs = "re"),
           family = tw(),
           method = "REML",
           data = mps_males)

#New dataset for predicting
NEW_DATA <- data.frame(Age = 1:100,
                       Sex = factor("Male"),
                       Tissue = factor("new"),
                       Article_number = factor("new"),
                       Method_analysis = factor("new"))

#Make predictions for that tissue type
mp_hat <- predict(FIT, newdata = NEW_DATA, type = "link", se.fit = TRUE, exclude = c("s(Article_number)","s(Method_analysis)","s(Tissue)"))
NEW_DATA$mp_hat <- exp(mp_hat$fit)
NEW_DATA$mp_hat_min <- exp(mp_hat$fit - mp_hat$se.fit*1.96)
NEW_DATA$mp_hat_max <- exp(mp_hat$fit + mp_hat$se.fit*1.96)

e <- 
  ggplot() +
  geom_point(data = mps_males, aes(x = Age, y = MPs_concentration_weight), col = "#de425b", size = 0.3) +
  geom_line(data = NEW_DATA, aes(x = Age, y = mp_hat), linewidth = 0.6, col = "#de425b") +
  geom_ribbon(data = NEW_DATA, aes(x = Age, ymin = mp_hat_min, ymax = mp_hat_max), alpha = 0.2, fill = "#de425b") +
  ggtitle("e") +
  xlab("Age (years)") +
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
        legend.title = element_blank(),
        legend.text = element_text(size=4, family = "sans", face = "bold"),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(0.2, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.position = "right",
        legend.position.inside = c(0.1,0.85),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm")) +
  scale_color_viridis(discrete = TRUE) +
  scale_y_log10(limits = c(0.1,100000), breaks = c(0.1,1,10,100,1000,10000,100000)) +
  #scale_x_continuous(expand = c(0.01,0)) +
  annotation_logticks(sides= "l",
                      short = unit(0.04, "cm"),
                      mid = unit(0.04, "cm"),
                      long = unit(0.1, "cm"),
                      linewidth = 0.2)



#---------------------------------------------------------------------
# Figure 4f - Relationship with age for females in the py-GC data
#---------------------------------------------------------------------

mps_females <- data[!is.na(data$MPs_concentration_weight),]
mps_females <- mps_females[mps_females$Sex == "Female",]

#Relationship with age for the spectroscopy data
FIT <- gam(MPs_concentration_weight ~  Age + s(Age, bs = "re") +
             s(Article_number, bs = "re") +
             s(Tissue, bs = "re"),
           family = tw(),
           method = "REML",
           data = mps_females)

#New dataset for predicting
NEW_DATA <- data.frame(Age = 1:100,
                       Sex = factor("Female"),
                       Tissue = factor("new"),
                       Article_number = factor("new"),
                       Method_analysis = factor("new"))

#Make predictions for that tissue type
mp_hat <- predict(FIT, newdata = NEW_DATA, type = "link", se.fit = TRUE, exclude = c("s(Article_number)","s(Method_analysis)","s(Tissue)"))
NEW_DATA$mp_hat <- exp(mp_hat$fit)
NEW_DATA$mp_hat_min <- exp(mp_hat$fit - mp_hat$se.fit*1.96)
NEW_DATA$mp_hat_max <- exp(mp_hat$fit + mp_hat$se.fit*1.96)

f <- 
  ggplot() +
  geom_point(data = mps_females, aes(x = Age, y = MPs_concentration_weight), col = "#184E77", size = 0.3) +
  geom_line(data = NEW_DATA, aes(x = Age, y = mp_hat), linewidth = 0.6, col = "#184E77") +
  geom_ribbon(data = NEW_DATA, aes(x = Age, ymin = mp_hat_min, ymax = mp_hat_max), alpha = 0.2, fill = "#184E77") +
  ggtitle("f") +
  xlab("Age (years)") +
  ylab ("MP concentration (ug/g)") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(linewidth = 0.2),
        axis.line.x = element_line(linewidth = 0.2),
        axis.ticks.x = element_line(linewidth = 0.2),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size=6, family = "sans", face = "bold"),
        axis.text.y = element_blank(),
        axis.text.x  = element_text(size=5, family = "sans", colour = "black"),
        strip.text.x = element_text(size=9, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 10, family = "sans", face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size=4, family = "sans", face = "bold"),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(0.2, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.position = "right",
        legend.position.inside = c(0.1,0.85),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm")) +
  scale_color_viridis(discrete = TRUE) +
  scale_y_log10(limits = c(0.1,100000), breaks = c(0.1,1,10,100,1000,10000,100000)) #+
  #scale_x_continuous(expand = c(0.01,0)) +
  # annotation_logticks(sides= "l",
  #                     short = unit(0.04, "cm"),
  #                     mid = unit(0.04, "cm"),
  #                     long = unit(0.1, "cm"),
  #                     linewidth = 0.2)



#---------------------------------------------------------------------
# Combine and save
#---------------------------------------------------------------------

TOP <-
  grid.arrange(a,b,c,
               ncol=3,
               nrow=1,
               widths = c(3,1.5,1.5))

BOT <-
  grid.arrange(d,e,f,
               ncol=3,
               nrow=1,
               widths = c(3,1.5,1.5))

fig4 <-
  grid.arrange(TOP,
               BOT,
               ncol=1,
               nrow=2)


#Save the figures
ggsave(fig4,
       width = 6.86, height = 4, units = "in",
       dpi = 600,
       bg = "transparent",
       file="figures/Figure_4.png")

