# Test file

# Load packages
library(plyr); library(abind); library(permute); library(lme4); 
library(R2WinBUGS); library(coda);library(vegan); library(arm); 
library(lattice); library(sciplot); library(reshape); library(cluster);
library(nortest); library(ggplot2);library(nlme); library(lme4);
library(scales); library(labdsv); library(stats); library(glmmML); 
library(MuMIn); library(effects); library(car); library(MASS); 
library(tidyverse); library(rsq); library(arm); library(plotrix); library(base);
library(sciplot); library(simr); library(readxl); library(emmeans);
library(multcomp);library(ggpubr); library(DHARMa); library(indicspecies);
library(readxl); library(glmmTMB); library(AICcmodavg);
library(devtools); library(grid); library(ggbiplot); library(forcats);
library(geosphere); library(mgcv); library(maps); library(ade4);
library(bipartite); library(metafor); library(SciViews); library(grDevices); 
library(DataCombine); library(cowplot); library(performance);


# Load in data 
data <- read_excel("Data/PecenkaPollen_DataForR.xlsx")
str(data)
attach(data)

## Run models -----
### Data checks
# select response variables
testvar = data$Tree
testvar = data$Crop
testvar = data$Clover
testvar = (data$HerbaceousFlower + data$Weed + data$Goldenrod)
# Exclude testvar = data$Grass
testvar = data$SR

# Check outliers of response
dotchart(testvar)
boxplot(testvar)
histogram(testvar)

# Models
# depending on the data, use poisson (SR), gaussian (H) or tweedie (in case of many zeros)
# testing a quadratic effect of time has no effect for any of the plant types
mod1 <- glmmTMB(testvar ~ scale(Time3) + Treatment + scale(Crop_3km) + (1|Site_ID/Year),
                family = gaussian, data = data)
mod1i <- glmmTMB(testvar ~ scale(Time3) + Treatment*scale(Crop_3km) + (1|Site_ID/Year),
                family = gaussian, data = data)
mod2 <- glmmTMB(testvar ~ scale(Time3) + Treatment + scale(Crop_1km) + (1|Site_ID/Year),
                family = gaussian, data = data)
mod2i <- glmmTMB(testvar ~ scale(Time3) + Treatment*scale(Crop_1km) + (1|Site_ID/Year),
                family = gaussian, data = data)
summary(mod1i)
summary(mod2i)

#No interaction, so drop the modxi
check_collinearity(mod1)
check_collinearity(mod2)
AICc(mod1)
AICc(mod2)
summary(mod1)
summary(mod2)
r.squaredGLMM(mod1)
r.squaredGLMM(mod2)
plot(allEffects(mod1))

# Model validation
res <- simulateResiduals(fittedModel = mod1, re.form = NULL)
# plot diagnostics 
plot(res)
testOutliers(res)
testDispersion(res)
# zero inflation test 
testZeroInflation(res)



## Plot ----
# Shannonn and crop cover
mod = glmmTMB(H ~ Time3 + Treatment + Crop_3km + (1|Site_ID) + (1|Year),
               family = gaussian, data = data)
summary(mod)
Predicttestvar<-as.data.frame(effect("Crop_3km",mod,xlevels=50))
Shannon = ggplot(data=data, aes(x=Crop_3km))+
  geom_ribbon(data=Predicttestvar, aes(ymin=lower, ymax=upper), fill="red",alpha=0.1)+
  geom_line(data=Predicttestvar, aes(x=Crop_3km, y=fit),color="red", alpha = 1)+
  geom_jitter(aes(y=H), show.legend = FALSE, shape = 1, size = 0.75, alpha = 0.5, 
              height = 0.1, width = 1) +
  xlab("Crop cover 3 km (%)")+
  ylab("Shannon diversity pollen")+
  theme_classic()
Shannon

# Trees and time
mod1 = glmmTMB(Tree ~ Time3 + Treatment + Crop_3km + (1|Site_ID/Year),
              family = tweedie, data = data)
summary(mod1)
Predicttestvar1<-as.data.frame(effect("Time3",mod1,xlevels=50))
Tree = ggplot(data=data, aes(x=Time3))+
  geom_ribbon(data=Predicttestvar1, aes(ymin=lower, ymax=upper), fill="red",alpha=0.1)+
  geom_line(data=Predicttestvar1, aes(x=Time3, y=fit),color="red", alpha = 1)+
  geom_jitter(aes(y=Tree), show.legend = FALSE, shape = 1, size = 0.75, alpha = 0.5, 
              height = 0.1, width = 1) +
  xlab("Julian day")+
  scale_y_continuous(limits = c(-0, 100))+
  ylab("Tree and shrub pollen (%)")+
  theme_classic()
Tree

# Clovers and time
mod2 = glmmTMB(Clover ~ Time3 + Treatment + Crop_3km + (1|Site_ID/Year),
               family = tweedie, data = data)
summary(mod2)
Predicttestvar2<-as.data.frame(effect("Time3",mod2,xlevels=50))
Clover = ggplot(data=data, aes(x=Time3))+
  geom_ribbon(data=Predicttestvar2, aes(ymin=lower, ymax=upper), fill="red",alpha=0.1)+
  geom_line(data=Predicttestvar2, aes(x=Time3, y=fit),color="red", alpha = 1)+
  geom_jitter(aes(y=Clover), show.legend = FALSE, shape = 1, size = 0.75, alpha = 0.5, 
              height = 0.1, width = 1) +
  xlab("Julian day")+
  scale_y_continuous(limits = c(-0, 100))+
  ylab("Clover pollen (%)")+
  theme_classic()
Clover

# Herbs and time
testvar = (data$HerbaceousFlower + data$Weed + data$Goldenrod)
mod3 = glmmTMB(testvar ~ Time3 + Treatment + Crop_3km + (1|Site_ID/Year),
               family = tweedie, data = data)
summary(mod3)
Predicttestvar3<-as.data.frame(effect("Time3",mod3,xlevels=50))
Herb = ggplot(data=data, aes(x=Time3))+
  geom_ribbon(data=Predicttestvar3, aes(ymin=lower, ymax=upper), fill="red",alpha=0.1)+
  geom_line(data=Predicttestvar3, aes(x=Time3, y=fit),color="red", alpha = 1)+
  geom_jitter(aes(y=testvar), show.legend = FALSE, shape = 1, size = 0.75, alpha = 0.5, 
              height = 0.1, width = 1) +
  xlab("Julian day")+
  scale_y_continuous(limits = c(-0, 100))+
  ylab("Herbaceous pollen (%)")+
  theme_classic()
Herb

# Plot figures together
plot_grid(Shannon, Tree, Clover, Herb, 
           labels = c('A', 'B', 'C', 'D'),  
           ncol = 2)
ggsave("figures/Fig1.png", width = 6, height = 6, dpi = 150)



## Stacked barplot 
str(data)

data2 = data %>%
  mutate(Herbaceous = (HerbaceousFlower+Weed+Goldenrod))

data3 = data2 %>%
  select(Site_id, Year, Site, Treatment, Time2, Tree, Crop, Clover, Herbaceous)

data4 <- data3 %>%  
  select(Site_id, Year, Site, Treatment, Time2, Tree, Crop, Clover, Herbaceous) %>% 
  group_by(Time2) %>% 
  summarise(across(Tree:Herbaceous, .fns = mean))

data5 <- data4 %>%
  pivot_longer(
  cols = c(Tree, Crop, Clover, Herbaceous),
  names_to = "Plant_Category",
  values_to = "Proportion_of_Sample"
)

# Define the desired order of studies
desired_order1 <- c("Late June", "Early July", "Late July", "Early August", 
                    "Late August", "Early September")  # Customize this order
# Convert taxa to a factor with the specified order
data5$time3 <- factor(data5$Time2, levels = desired_order1)

# Define the desired order of taxa
desired_order2 <- c("Tree", "Crop", "Clover", "Herbaceous")  # Customize this order
# Convert taxa to a factor with the specified order
data5$plant_category <- factor(data5$Plant_Category, levels = desired_order2)

# Define custom colors for each taxa
taxa_colors <- c("Tree" = "darkgreen", "Crop" = "yellow", 
                 "Clover" = "#7570b3", "Herbaceous" = "brown")

# # Define custom names for the legend
# taxa_labels <- c("HB" = "Honingbij", "BB" = "Hommels", "SB" = "Solitaire bijen", 
#                  "HF" = "Zweefvliegen", "BF" = "Vlinders", "WA" = "Wespen")

# Create stacked bar plot
ggplot(data5, aes(x = time3, y = Proportion_of_Sample, fill = plant_category)) +
  geom_bar(stat = "identity") +
  labs(x = "Time period", y = "Plant type in pollen (%)", fill = "Plant type") +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = taxa_colors)  # Manually assign colors
ggsave("figures/Fig2.png", width = 6, height = 4, dpi = 150)
