#This script generates the descriptive analyses of the dataset
# used throughout the main text. 

# Last updated: November 14 2024

library(mgcv)
library(emmeans)

# Load data
data <- read.csv("data/MPs_human_tissues.csv", header = TRUE)
data[which(data$Tissue == "Olfactory_bulb"),"Tissue"] <- "Olfactory bulb"
data$Tissue <- as.factor(data$Tissue)
data$Article_number <- as.factor(data$Article_number)
data$Method_analysis <- as.factor(data$Method_analysis)


#---------------------------------------------------------------------
# Descriptive statistics of the study cohort 
#---------------------------------------------------------------------

#Which countries are the data from
table(data$Country)
round(prop.table(table(data$Country)),3)*100

#Which tissues are the data from
table(data$Tissue)
round(prop.table(table(data$Tissue)),3)*100

#Sex ratio of the data
table(data$Sex)
SEX <- aggregate(Species ~ Sex + Tissue, data = data, FUN = "length")
SEX <- reshape(SEX, idvar = "Sex", timevar = "Tissue", direction = "wide")
SEX <- SEX[,!(names(SEX) %in% c("Species.Penis","Species.Prostate","Species.Testis"))]
SEX[is.na(SEX)] <- 0
mean(unlist(SEX[1,-1])/unlist(SEX[2,-1]))
sd(unlist(SEX[1,-1])/unlist(SEX[2,-1]))

#Age distribution
mean(data$Age)
sd(data$Age)
range(data$Age)

#---------------------------------------------------------------------
# Microplastic analysis methods
#---------------------------------------------------------------------

#Distribution of analytical methods
table(data$Method_analysis)

#Detection rates of the different methods
Py_GC <- data[data$Method_analysis == "Py-GC","MPs_concentration_weight"]
which(Py_GC == 0)
sum(Py_GC > 0)/length(Py_GC)

Raman <- data[data$Method_analysis == "Raman","MPs_concentration_particles_g"]
round(sum(Raman > 0)/length(Raman),3)*100

LDIR <- data[data$Method_analysis == "LDIR","MPs_concentration_particles_g"]
round(sum(LDIR > 0)/length(LDIR),3)*100

uFTIR <- data[data$Method_analysis == "uFTIR","MPs_concentration_particles_g"]
round(sum(uFTIR > 0)/length(uFTIR),3)*100


#Concentration by analysis method
aggregate(MPs_concentration_particles_g ~ Method_analysis,
          data = data[!is.na(data$MPs_concentration_particles_g),],
          FUN = "median")

FIT <- gam(MPs_concentration_particles_g ~  Method_analysis +
             s(Article_number, bs = "re") +
             s(Tissue, bs = "re"),
           family = tw(),
           method = "REML",
           data = data[!is.na(data$MPs_concentration_particles_g),])

summary(FIT)


#---------------------------------------------------------------------
# Microplastic polymer abundance between different tissues
#---------------------------------------------------------------------

#Distribution of most prevalent polymers
table(data$Most_prevalent_MP)

aggregate(Most_prevalent_MP ~ Tissue, data = data, FUN = "length")

#---------------------------------------------------------------------
# Microplastics concentrations in different tissues
#---------------------------------------------------------------------

#Concentration by tissue type
aggregate(MPs_concentration_particles_g ~ Tissue,
          data = data[!is.na(data$MPs_concentration_particles_g),],
          FUN = "mean")

aggregate(MPs_concentration_particles_g ~ Tissue,
          data = data[!is.na(data$MPs_concentration_particles_g),],
          FUN = "sd")

#Test for differences in MP concentrations between tissue types in the spectroscopy data
fit <- gam(MPs_concentration_particles_g ~ Tissue +
             s(Article_number, bs = "re") +
             s(Method_analysis, bs = "re"),
           family = tw(link = "log"),
           data = data[!is.na(data$MPs_concentration_particles_g),],
           method = "REML")

summary(fit)

#Conduct Tukey post-hoc pairwise comparisons
pairs(emmeans(fit, ~ Tissue), adjust = "tukey")


#Test for differences in MP concentrations between tissue types in the pyrolysis-GC data
fit <- gam(MPs_concentration_weight ~ Tissue +
             s(Article_number, bs = "re"),
           family = tw(link = "log"),
           data = data[!is.na(data$MPs_concentration_weight),],
           method = "REML")

summary(fit)

#Conduct Tukey post-hoc pairwise comparisons
pairs(emmeans(fit, ~ Tissue), adjust = "tukey")


#Subset the sex specific tissues and those with only one sex
data2 <- data[!(data$Tissue %in% c("Penis","Prostate","Testis","Carotid artery","Myocardium","Ileum")),]

#Test for differences in MP concentrations between sexes in the different tissue types in the spectroscopy data
fit <- gam(MPs_concentration_particles_g ~ Tissue*Sex +
             s(Article_number, bs = "re") +
             s(Method_analysis, bs = "re"),
           family = tw(link = "log"),
           data = data2[!is.na(data2$MPs_concentration_particles_g),],
           method = "REML")

summary(fit)

#Conduct Tukey post-hoc pairwise comparisons
pairs(emmeans(fit, ~ Sex | Tissue), adjust = "tukey")


#Test for differences in MP concentrations between sexes in the different tissue types in the pyrolysis-GC data
fit <- gam(MPs_concentration_weight ~ Tissue*Sex +
             s(Article_number, bs = "re"),
           family = tw(link = "log"),
           data = data2[!is.na(data2$MPs_concentration_weight),],
           method = "REML")

summary(fit)

#Conduct Tukey post-hoc pairwise comparisons
pairs(emmeans(fit, ~ Sex | Tissue), adjust = "tukey")


#---------------------------------------------------------------------
# Microplastics concentrations and age
#---------------------------------------------------------------------

#Median concentration across all samples for the spectroscopy data
MPs_particles <- data[!is.na(data$MPs_concentration_particles_g),"MPs_concentration_particles_g"]
median(MPs_particles)

#CIs on the median
n <- length(MPs_particles)
sort(MPs_particles)[round((n/2)*(1 + (1.96)/sqrt(n)),0)]
sort(MPs_particles)[round((n/2)*(1 - (1.96)/sqrt(n)),0)]


#Median concentration across all samples for the pyrolysis-GC data
MPs_particles <- data[!is.na(data$MPs_concentration_weight),"MPs_concentration_weight"]
median(MPs_particles)

#CIs on the median
n <- length(MPs_particles)
sort(MPs_particles)[round((n/2)*(1 + (1.96)/sqrt(n)),0)]
sort(MPs_particles)[round((n/2)*(1 - (1.96)/sqrt(n)),0)]

#Relationship with age for the spectroscopy data
FIT <- gam(MPs_concentration_particles_g ~  Age + s(Age, Tissue, bs = "re") +
             s(Article_number, bs = "re") +
             s(Tissue, bs = "re") +
             s(Method_analysis, bs = "re"),
           family = tw(),
           method = "REML",
           data = data[!is.na(data$MPs_concentration_particles_g),])

summary(FIT)

#Relationship with age for the pyrolysis-GC data
FIT <- gam(MPs_concentration_weight ~  Age + s(Age, Tissue, bs = "re") +
             s(Article_number, bs = "re") +
             s(Tissue, bs = "re"),
           family = tw(),
           method = "REML",
           data = data[!is.na(data$MPs_concentration_weight),])

summary(FIT)


#Estimated concentrations at different ages

#New dataset for predicting
NEW_DATA <- data.frame(Age = c(10,60),
                       Tissue = factor("new"),
                       Article_number = factor("new"),
                       Method_analysis = factor("new"))

#Make predictions for increase over 50 years
mp_hat <- predict(FIT, newdata = NEW_DATA, type = "response", se.fit = TRUE, exclude = c("s(Article_number)"))
mp_hat$fit[2] - mp_hat$fit[1]
1.96*sqrt(sum(mp_hat$se.fit^2))

(mp_hat$fit[2] - mp_hat$fit[1])/sd(na.omit(data$MPs_concentration_weight))*100






#-------------------------------------
#Effects by sex


mps_males <- data[!is.na(data$MPs_concentration_particles_g),]
mps_males <- mps_males[mps_males$Sex == "Male",]

#Relationship with age for males from the spectroscopy data
FIT <- gam(MPs_concentration_particles_g ~  Age + s(Age, bs = "re") +
             s(Article_number, bs = "re") +
             s(Tissue, bs = "re") +
             s(Method_analysis, bs = "re"),
           family = tw(),
           method = "REML",
           data = mps_males)

summary(FIT)

mps_females <- data[!is.na(data$MPs_concentration_particles_g),]
mps_females <- mps_females[mps_females$Sex == "Female",]

#Relationship with age for females from the spectroscopy data
FIT <- gam(MPs_concentration_particles_g ~  Age + s(Age, bs = "re") +
             s(Article_number, bs = "re") +
             s(Tissue, bs = "re") +
             s(Method_analysis, bs = "re"),
           family = tw(),
           method = "REML",
           data = mps_females)

summary(FIT)


mps_males <- data[!is.na(data$MPs_concentration_weight),]
mps_males <- mps_males[mps_males$Sex == "Male",]

#Relationship with age for the spectroscopy data
FIT <- gam(MPs_concentration_weight ~  Age + s(Age, bs = "re") +
             s(Article_number, bs = "re") +
             s(Tissue, bs = "re"),
           family = tw(),
           method = "REML",
           data = mps_males)

summary(FIT)

mps_females <- data[!is.na(data$MPs_concentration_weight),]
mps_females <- mps_females[mps_females$Sex == "Female",]

#Relationship with age for the spectroscopy data
FIT <- gam(MPs_concentration_weight ~  Age + s(Age, bs = "re") +
             s(Article_number, bs = "re") +
             s(Tissue, bs = "re"),
           family = tw(),
           method = "REML",
           data = mps_females)

summary(FIT)


#-------------------------------------
#Relationship with age in male-specific tissues

males <- data[(data$Tissue %in% c("Penis","Prostate","Testis")),]


#Relationship with age for the spectroscopy data
FIT <- gam(MPs_concentration_particles_g ~  Age + s(Age, Tissue, bs = "re") +
             s(Article_number, bs = "re") +
             s(Tissue, bs = "re"),
           family = tw(),
           method = "REML",
           data = males[!is.na(males$MPs_concentration_particles_g),])

summary(FIT)

#Relationship with age for the pyrolysis-GC data
FIT <- gam(MPs_concentration_weight ~  Age,
           family = tw(),
           method = "REML",
           data = males[!is.na(males$MPs_concentration_weight),])

summary(FIT)

