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

### Data checks
# select response variables
testvar = data$Tree
testvar = data$Crop
testvar = data$Clover
testvar = (data$HerbaceousFlower + data$Weed + data$Goldenrod)
# Exclude testvar = data$Grass

# Check outliers of response
dotchart(testvar)
boxplot(testvar)
histogram(testvar)

# Models
# depending on the data, use gaussian or tweedie (in case of many zeros)
# testing a quadratic effect of time has no effect for any of the plant types
mod1 <- glmmTMB(testvar ~ scale(Time3) + Treatment + scale(Crop_3km) + (1|Site_ID/Year),
                family = tweedie, data = data)
mod1i <- glmmTMB(testvar ~ scale(Time3) + Treatment*scale(Crop_3km) + (1|Site_ID/Year),
                family = tweedie, data = data)
mod2 <- glmmTMB(testvar ~ scale(Time3) + Treatment + scale(Crop_1km) + (1|Site_ID/Year),
                family = tweedie, data = data)
mod2i <- glmmTMB(testvar ~ scale(Time3) + Treatment*scale(Crop_1km) + (1|Site_ID/Year),
                family = tweedie, data = data)
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


### Plot
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
ggsave("Fig1.png", width = 3.5, height = 3, dpi = 300)
