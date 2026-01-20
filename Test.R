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

