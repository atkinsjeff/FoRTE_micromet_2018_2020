# analysis of variance scripts
require(tidyverse)

# import data
df <- read.csv("analysis_ready.csv")

df


#### 
summary(aov(temp ~ disturbance_severity * treatment, data = df))
summary(aov(temp.sd ~ disturbance_severity * treatment, data = df))
summary(aov(vwc ~ disturbance_severity * treatment, data = df))
summary(aov(vwc.sd ~ disturbance_severity * treatment, data = df))
summary(aov(par ~ disturbance_severity * treatment, data = df))
