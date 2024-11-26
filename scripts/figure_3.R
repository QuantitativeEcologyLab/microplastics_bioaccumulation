# This script generates figure 2 from PAPER TITLE
# It first generates the following sub-pannels:
# a) Tukey's post-hoc tests for the spectroscopy data
# b) Tukey's post-hoc tests for the Pyrolysis-GC data
# c) Sex-based differences in MPs within each tissue for the spectroscopy data
# d) Sex-based differences in MPs within each tissue for the Pyrolysis-GC data
# it then combines the 4 sub-panels into a single figure

# Last updated: November 14 2024

#Import the packages
library(ggplot2)
library(gridExtra)
library(mgcv)
library(emmeans)
library(tidyr)
library(dplyr)

# Load data
data <- read.csv("data/MPs_human_tissues.csv", header = TRUE)
data[which(data$Tissue == "Olfactory_bulb"),"Tissue"] <- "Olfactory bulb"
data$Tissue <- as.factor(data$Tissue)
data$Sex <- as.factor(data$Sex)
data$Article_number <- as.factor(data$Article_number)
data$Method_analysis <- as.factor(data$Method_analysis)

#---------------------------------------------------------------------
# Figure 3a - Tukey's post-hoc tests for the spectroscopy data
#---------------------------------------------------------------------

#Test for differences in MP concentrations between tissue types
fit <- gam(MPs_concentration_particles_g ~ Tissue + s(Article_number, bs = "re") + s(Method_analysis, bs = "re"),
           family = tw(link = "log"),
           data = data[!is.na(data$MPs_concentration_particles_g),],
           method = "REML")

summary(fit)

#Obtain Estimated Marginal Means for 'Tissue'
emm_tissue <- emmeans(fit, ~ Tissue)

#Conduct Tukey post-hoc pairwise comparisons
pairwise_comparisons <- pairs(emm_tissue, adjust = "tukey")

# Convert the pairwise comparisons to a data frame
comparison_df <- as.data.frame(pairwise_comparisons)

# Calculate 95% Confidence Intervals and Determine Significance Levels
comparison_df <- comparison_df %>%
  mutate(
    CI_lower = estimate - qt(0.975, df) * SE,
    CI_upper = estimate + qt(0.975, df) * SE,
    # Define significance based on p-value thresholds
    Significance = case_when(
      p.value < 0.001 ~ "p < 0.001",
      p.value < 0.01  ~ "p < 0.01",
      p.value < 0.05  ~ "p < 0.05",
      TRUE            ~ "ns"  # 'ns' stands for not significant
    )
  )


# Generate the figure
a <-
  ggplot(comparison_df, aes(x = reorder(contrast, estimate), y = estimate, color = Significance)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_point(size = 0.4) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), linewidth = 0.2, width = 0.3) +
  scale_color_manual(values = c("p < 0.001" = "darkred", 
                                "p < 0.01"  = "red", 
                                "p < 0.05"  = "orange", 
                                "ns"        = "black"),
                     name = "Significane") +
  ggtitle("a") +
  ylab ("Difference in log(MP concentration)") +
  coord_flip() +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(linewidth = 0.2),
        axis.line.x = element_line(linewidth = 0.2),
        axis.ticks.x = element_line(linewidth = 0.2),
        axis.ticks.y = element_line(linewidth = 0.2),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size=6, family = "sans", face = "bold"),
        axis.text.y = element_text(size=3, family = "sans"),
        axis.text.x  = element_text(size=5, family = "sans", colour = "black"),
        strip.text.x = element_text(size=9, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 10, family = "sans", face = "bold"),
        legend.title = element_text(size=6, family = "sans", face = "bold"),
        legend.text = element_text(size=6, family = "sans", face = "bold"),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(0.2, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.position = "inside",
        legend.position.inside = c(0.15,0.9),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm"))


#---------------------------------------------------------------------
# Figure 3b - Tukey's post-hoc tests for the pyrolysis-GC data
#---------------------------------------------------------------------

#Test for differences in MP concentrations between tissue types
fit <- gam(MPs_concentration_weight ~ Tissue + s(Article_number, bs = "re"),
           family = tw(link = "log"),
           data = data[!is.na(data$MPs_concentration_weight),],
           method = "REML")

summary(fit)

#Obtain Estimated Marginal Means for 'Tissue'
emm_tissue <- emmeans(fit, ~ Tissue)

#Conduct Tukey post-hoc pairwise comparisons
pairwise_comparisons <- pairs(emm_tissue, adjust = "tukey")

# Convert the pairwise comparisons to a data frame
comparison_df <- as.data.frame(pairwise_comparisons)

# Calculate 95% Confidence Intervals and Determine Significance Levels
comparison_df <- comparison_df %>%
  mutate(
    CI_lower = estimate - qt(0.975, df) * SE,
    CI_upper = estimate + qt(0.975, df) * SE,
    # Define significance based on p-value thresholds
    Significance = case_when(
      p.value < 0.001 ~ "p < 0.001",
      p.value < 0.01  ~ "p < 0.01",
      p.value < 0.05  ~ "p < 0.05",
      TRUE            ~ "ns"  # 'ns' stands for not significant
    )
  )


# Generate the figure
b <-
  ggplot(comparison_df, aes(x = reorder(contrast, estimate), y = estimate, color = Significance)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_point(size = 0.6) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), linewidth = 0.2, width = 0.3) +
  scale_color_manual(values = c("p < 0.001" = "darkred", 
                                "p < 0.01"  = "red", 
                                "p < 0.05"  = "orange", 
                                "ns"        = "black"),
                     name = "Significane") +
  ggtitle("b") +
  ylab ("Difference in log(MP concentration)") +
  coord_flip() +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(linewidth = 0.2),
        axis.line.x = element_line(linewidth = 0.2),
        axis.ticks.x = element_line(linewidth = 0.2),
        axis.ticks.y = element_line(linewidth = 0.2),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size=6, family = "sans", face = "bold"),
        axis.text.y = element_text(size=5, family = "sans"),
        axis.text.x  = element_text(size=5, family = "sans", colour = "black"),
        strip.text.x = element_text(size=9, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 10, family = "sans", face = "bold"),
        legend.title = element_text(size=6, family = "sans", face = "bold"),
        legend.text = element_text(size=6, family = "sans", face = "bold"),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(0.2, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.position = "none",
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm"))





#---------------------------------------------------------------------
# Figure 3c - Tukey's post-hoc tests on sex for the spectroscopy data
#---------------------------------------------------------------------
#Subset the sex specific tissues and those with only one sex
data2 <- data[!(data$Tissue %in% c("Penis","Prostate","Testis","Carotid artery","Myocardium","Ileum")),]

#Test for differences in MP concentrations between tissue types
fit <- gam(MPs_concentration_particles_g ~ Tissue*Sex + s(Article_number, bs = "re") + s(Method_analysis, bs = "re"),
           family = tw(link = "log"),
           data = data2[!is.na(data2$MPs_concentration_particles_g),],
           method = "REML")

summary(fit)

#Obtain Estimated Marginal Means for 'Tissue'
emm_sex_within_tissue <- emmeans(fit, ~ Sex | Tissue)

#Conduct Tukey post-hoc pairwise comparisons
pairwise_comparisons <- pairs(emm_sex_within_tissue, adjust = "tukey")

# Convert the pairwise comparisons to a data frame
comparison_df <- as.data.frame(pairwise_comparisons)

# Calculate 95% Confidence Intervals and Determine Significance Levels
comparison_df <- comparison_df %>%
  mutate(
    CI_lower = estimate - qt(0.975, df) * SE,
    CI_upper = estimate + qt(0.975, df) * SE,
    # Define significance based on p-value thresholds
    Significance = case_when(
      p.value < 0.001 ~ "p < 0.001",
      p.value < 0.01  ~ "p < 0.01",
      p.value < 0.05  ~ "p < 0.05",
      TRUE            ~ "ns"  # 'ns' stands for not significant
    )
  )


# Generate the figure
c <-
  ggplot(comparison_df, aes(x = reorder(Tissue, estimate), y = estimate, color = Significance)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_point(size = 0.4) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), linewidth = 0.2, width = 0.3) +
  scale_color_manual(values = c("p < 0.001" = "darkred", 
                                "p < 0.01"  = "red", 
                                "p < 0.05"  = "orange", 
                                "ns"        = "black"),
                     name = "Significane") +
  ggtitle("c") +
  ylab ("Difference in log(MP concentration)") +
  coord_flip() +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(linewidth = 0.2),
        axis.line.x = element_line(linewidth = 0.2),
        axis.ticks.x = element_line(linewidth = 0.2),
        axis.ticks.y = element_line(linewidth = 0.2),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size=6, family = "sans", face = "bold"),
        axis.text.y = element_text(size=5, family = "sans"),
        axis.text.x  = element_text(size=5, family = "sans", colour = "black"),
        strip.text.x = element_text(size=9, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 10, family = "sans", face = "bold"),
        legend.title = element_text(size=6, family = "sans", face = "bold"),
        legend.text = element_text(size=6, family = "sans", face = "bold"),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(0.2, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.position = "none",
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm"))



#---------------------------------------------------------------------
# Figure 3d - Tukey's post-hoc tests on sex for the spectroscopy data
#---------------------------------------------------------------------
aggregate(MPs_concentration_weight ~ Tissue +Sex, data = data2[!is.na(data2$MPs_concentration_weight),], FUN = "mean")

#Test for differences in MP concentrations between tissue types
fit <- gam(MPs_concentration_weight ~ Tissue*Sex+ s(Article_number, bs = "re"),
           family = tw(link = "log"),
           data = data2[!is.na(data2$MPs_concentration_weight),],
           method = "REML")

summary(fit)

#Obtain Estimated Marginal Means for 'Tissue'
emm_sex_within_tissue <- emmeans(fit, ~ Sex | Tissue)

#Conduct Tukey post-hoc pairwise comparisons
pairwise_comparisons <- pairs(emm_sex_within_tissue, adjust = "tukey")

# Convert the pairwise comparisons to a data frame
comparison_df <- as.data.frame(pairwise_comparisons)

# Calculate 95% Confidence Intervals and Determine Significance Levels
comparison_df <- comparison_df %>%
  mutate(
    CI_lower = estimate - qt(0.975, df) * SE,
    CI_upper = estimate + qt(0.975, df) * SE,
    # Define significance based on p-value thresholds
    Significance = case_when(
      p.value < 0.001 ~ "p < 0.001",
      p.value < 0.01  ~ "p < 0.01",
      p.value < 0.05  ~ "p < 0.05",
      TRUE            ~ "ns"  # 'ns' stands for not significant
    )
  )


# Generate the figure
d <-
  ggplot(comparison_df, aes(x = reorder(Tissue, estimate), y = estimate, color = Significance)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_point(size = 0.4) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), linewidth = 0.2, width = 0.3) +
  scale_color_manual(values = c("p < 0.001" = "darkred", 
                                "p < 0.01"  = "red", 
                                "p < 0.05"  = "orange", 
                                "ns"        = "black"),
                     name = "Significane") +
  ggtitle("d") +
  ylab ("Difference in log(MP concentration)") +
  coord_flip() +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(linewidth = 0.2),
        axis.line.x = element_line(linewidth = 0.2),
        axis.ticks.x = element_line(linewidth = 0.2),
        axis.ticks.y = element_line(linewidth = 0.2),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size=6, family = "sans", face = "bold"),
        axis.text.y = element_text(size=5, family = "sans"),
        axis.text.x  = element_text(size=5, family = "sans", colour = "black"),
        strip.text.x = element_text(size=9, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 10, family = "sans", face = "bold"),
        legend.title = element_text(size=6, family = "sans", face = "bold"),
        legend.text = element_text(size=6, family = "sans", face = "bold"),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(0.2, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.position = "none",
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm"))



#---------------------------------------------------------------------
# Combine and save
#---------------------------------------------------------------------


TOP <-
  grid.arrange(a,b,
               ncol=2,
               nrow=1)



BOT <-
  grid.arrange(c,d,
               ncol=2,
               nrow=1) 

fig3 <-
  grid.arrange(TOP,
               BOT,
               ncol=1,
               nrow=2,
               heights=c(3,1.5))


#Save the figures
ggsave(fig3,
       width = 6.86, height = 7, units = "in",
       dpi = 600,
       bg = "transparent",
       file="figures/Figure_3.png")

