# compositional change


# dependencies
require(fortedata)
require(vegan)
require(tidyverse)
require(ggplot2)
require(randomForest)
require(rfUtilities)
require(caret)
require(car)
require(Hmisc)
# disturbance 
# bring in metadata
df <- fortedata::fd_plot_metadata()

df <- data.frame(df)
df$subplot_id <- paste(df$replicate, 0, df$plot, df$subplot, sep = "")
df$subplot_id <- as.factor(df$subplot_id)

df %>%
    select(subplot_id, disturbance_severity, treatment) %>%
    distinct() %>%
    data.frame() -> dis.meta.data


dis.meta.data <- dis.meta.data[c(1:32), ]

## bring in the inventory and mortalitiy data sets
mor <- fd_mortality()
inv <- fd_inventory()

# make factors
mor$subplot_id <- as.factor(mor$subplot_id)

# 

# making the PFT for mor
# 
# 
# 
# 10 - Temperate broadleaf, mid-successional
# FAGR, ACRU, ACPE, ACSA
# 11 - Temperate broadleaf, late successional
# QURU, TCSA
# 
# 6 - Northern North American pines
# PIST, PIRE
# 9 - Temperate broadleaf, early successional
# AMEL, BEAL, POTR, POGR, BEPA


mor %>%
    mutate(pft = case_when(species == "FAGR" ~ 'TBF',
                           species == "ACRU" ~ 'TBF',
                           species == "ACPE" ~ 'TBF',
                           species == "ACSA" ~ 'TBF',
                           species == "QURU" ~ 'TBF',
                           species == "TCSA" ~ 'TBF',
                           species == "PIST" ~ 'NAP',
                           species == "PIRE" ~ "NAP",
                           species == 'AMEL' ~ 'TBF',
                           species == "BEAL" ~ 'TBF',
                           species == 'POTR' ~ 'TBF',
                           species == 'POGR' ~ 'TBF',
                           species == 'BEPA' ~ 'TBF')) %>%
    mutate(stage = case_when(species == "FAGR" ~ 'MID',
                             species == "ACRU" ~ 'MID',
                             species == "ACPE" ~ 'MID',
                             species == "ACSA" ~ 'MID',
                             species == "QURU" ~ 'LATE',
                             species == "TCSA" ~ 'LATE',
                             species == "PIST" ~ 'LATE',
                             species == "PIRE" ~ 'EARLY',
                             species == 'AMEL' ~ 'EARLY',
                             species == "BEAL" ~ 'EARL',
                             species == 'POTR' ~ 'EARLY',
                             species == 'POGR' ~ 'EARLY',
                             species == 'BEPA' ~ 'EARLY')) %>%
    data.frame() -> mor


# make basal area
mor$basal_area <- pi * (mor$dbh_cm / 2)^2

mor %>%
    group_by(subplot_id, fate) %>%
    dplyr::summarize(ba_plot = sum(basal_area)) -> subplot.ba

mor %>%
    group_by(subplot_id, species, stage, fate) %>%
    dplyr::summarize(ba_species = sum(basal_area)) -> species.ba

mor %>%
    group_by(subplot_id, fate) %>%
    dplyr::summarize(no_stems = n()) -> subplot.stems.2020

# gets the count of stems per species
inv %>%
    filter(health_status == "L") %>%
    group_by(subplot_id, species) %>%
    dplyr::summarize(no_stems_species = n()) %>%
    data.frame() -> species.stems.2018

mor %>%
    filter(health_status == "L") %>%
    filter(fate == "live") %>%
    group_by(subplot_id, species) %>%
    dplyr::summarize(no_stems_species = n()) %>%
    data.frame() -> species.stems.2020

# this make the species richness values for 2020, i.e. the ones that are still alive
mor %>% 
    filter(health_status == "L") %>%
    filter(fate == "live") %>%
    group_by(subplot_id) %>%
    dplyr::summarize(richness = n_distinct(species)) -> subplot.rich.2020

subplot.rich.2020$year <- 2020

inv %>%
    filter(health_status == "L") %>%
    group_by(subplot_id) %>%
    dplyr::summarize(richness = n_distinct(species)) -> subplot.rich.2018

subplot.rich.2018$year <- 2018

# now to compbine
species.richness <- rbind(subplot.rich.2018, subplot.rich.2020)





### species diversity with Shannon-Weinter
# first the community matrix for 2018
species.matrix.18 <- spread(species.stems.2018, species, no_stems_species)
species.matrix.18[is.na(species.matrix.18)] <- 0

# second the community matrix for 2020
species.matrix.20 <- spread(species.stems.2020, species, no_stems_species)
species.matrix.20[is.na(species.matrix.20)] <- 0

# make the values
shannon.18 <- diversity(species.matrix.18[, c(2:12)])
shannon.20 <- diversity(species.matrix.20[, c(2:12)])

subplots <- species.matrix.18$subplot_id

# reorganizing
shannon.2018 <- data.frame(subplots, shannon.18, year = 2018)
colnames(shannon.2018) <- c("subplot_id", "sw.index", "year")

shannon.2020 <- data.frame(subplots, shannon.20, year = 2020)
colnames(shannon.2020) <- c("subplot_id", "sw.index", "year")

shannon <- rbind(shannon.2018, shannon.2020)




# species evennness and bringing it all together
biodiversity <- merge(species.richness, shannon)

# replacing the shannon weiner 0 with 0.000001
biodiversity$evenness <- biodiversity$sw.index / log(as.numeric(biodiversity$richness))
biodiversity[biodiversity == 0] <- 0.001

# replace with 0
biodiversity[is.na(biodiversity)] <- 0
#biodiversity[sapply(biodiversity, is.infinite)] <- 0

# bio wide
biodiversity %>%
    pivot_wider(names_from = "year",
                values_from = c("richness", "sw.index", "evenness")) %>%
    data.frame() -> bio.2020

#### writing the data
#write.csv(biodiversity, "./data/forte_biodiversity_2018_2020.csv")

# bringin meta
bio.2020 <- merge(bio.2020, dis.meta.data)
bio.2020$txt <- as.factor(paste(bio.2020$disturbance_severity, bio.2020$treatment, sep = ""))

bio.2020$evenness_2020[bio.2020$evenness_2020 == 0] <- 0.0001


###### calculate difference from rep level control
bio.2020$replicate <- substr(bio.2020$subplot_id, 0, 1)

# find control values
bio.2020 %>%
    filter(disturbance_severity == 0) -> controls


# 
summary(aov(richness_2020 ~ replicate * disturbance_severity * treatment, data = bio.2020))

x11()
ggplot(bio.2020, aes(x = replicate, y = richness_2020))+
    geom_boxplot()+
    facet_grid(disturbance_severity ~ treatment)



# make the deltas
bio.2020$drichness <- log (bio.2020$richness_2020 / bio.2020$richness_2018)
bio.2020$devenness <- log(bio.2020$evenness_2020 / bio.2020$evenness_2018)
bio.2020$dsw.index <- log(bio.2020$sw.index_2020 / bio.2020$sw.index_2018)




bio.2020 %>%
    select(subplot_id, disturbance_severity, treatment, txt, drichness, dsw.index, devenness) %>%
    data.frame() -> bio.delta

# test for multicollinearity 
multi.collinear(bio.delta[, c(5:7)])

#### this tests each severity level independtivly to detect differences among treatments
#### 0%
bio.delta %>% 
    filter(disturbance_severity == 0) %>%
    data.frame() -> b.0

multi.collinear(b.0[, c(5:7)])

# testing biodiversity model by YEAR

rf.bio.imp.0 <- rf.modelSel(b.0[,c(5:6)], as.factor(b.0$treatment), 
                            seed = 666, 
                            imp.scale="se",
                            final.model = TRUE) 

x11()
plot(rf.bio.imp.0, imp = "sel")
rf.bio.0 <- randomForest(as.factor(treatment) ~ drichness + dsw.index, 
                         data = b.0, 
                         ntree = 2000,
                         importance = TRUE)

# test for signfincance
rf.significance(rf.bio.0, b.0[, c(5:6)], ntree = 2000, nperm = 999)


#### 45%
bio.delta %>% 
    filter(disturbance_severity == 45) %>%
    data.frame() -> b.45

multi.collinear(b.45[, c(5:7)])

# testing biodiversity model by YEAR

rf.bio.imp.45 <- rf.modelSel(b.45[,c(5:7)], as.factor(b.45$treatment), 
                             seed = 666, 
                             imp.scale="se",
                             final.model = TRUE) 

x11()
plot(rf.bio.imp.45, imp = "sel")
rf.bio.45 <- randomForest(as.factor(treatment) ~ devenness + dsw.index, 
                          data = b.45, 
                          ntree = 2000,
                          importance = TRUE)

# test for signfincance
rf.significance(rf.bio.45, b.45[, c(6:7)], ntree = 2000, nperm = 999)


#### 65%
bio.delta %>% 
    filter(disturbance_severity == 65) %>%
    data.frame() -> b.65

multi.collinear(b.65[, c(5:7)])

# testing biodiversity model by YEAR

rf.bio.imp.65 <- rf.modelSel(b.65[,c(5:7)], as.factor(b.65$treatment), 
                             seed = 666, 
                             imp.scale="se",
                             final.model = TRUE) 

x11()
plot(rf.bio.imp.65, imp = "sel")
rf.bio.65 <- randomForest(as.factor(treatment) ~ drichness + dsw.index, 
                          data = b.65, 
                          ntree = 2000,
                          importance = TRUE)

# test for signfincance
rf.significance(rf.bio.65, b.65[, c(5:6)], ntree = 2000, nperm = 999)



#### 85%
bio.delta %>% 
    filter(disturbance_severity == 85) %>%
    data.frame() -> b.85

multi.collinear(b.85[, c(5:7)])

# testing biodiversity model by YEAR

rf.bio.imp.85 <- rf.modelSel(b.85[,c(5:7)], as.factor(b.85$treatment), 
                             seed = 666, 
                             imp.scale="se",
                             final.model = TRUE) 

x11()
plot(rf.bio.imp.85, imp = "sel")
rf.bio.85 <- randomForest(as.factor(treatment) ~ drichness + dsw.index, 
                          data = b.85, 
                          ntree = 2000,
                          importance = TRUE)

# find significance
rf.bio.85
rf.significance(rf.bio.85, b.85[, c(5:6)], ntree = 2000, nperm = 999)





###### now doing by severity level with in treatment
# BOTTOM
bio.delta %>% 
    filter(treatment == "B") %>%
    data.frame() -> bottom

multi.collinear(bottom[, c(5:7)])

# testing biodiversity model by YEAR

rf.bio.imp.bottom <- rf.modelSel(bottom[,c(5:7)], as.factor(bottom$disturbance_severity), 
                                 seed = 666, 
                                 imp.scale="se",
                                 final.model = TRUE) 

x11()
plot(rf.bio.imp.bottom, imp = "sel")
rf.bio.bottom <- randomForest(as.factor(disturbance_severity) ~ drichness + dsw.index, 
                              data = bottom, 
                              ntree = 2000,
                              importance = TRUE)

# find significance
rf.bio.bottom
rf.significance(rf.bio.bottom, bottom[, c(5:6)], ntree = 2000, nperm = 999)




# TOP
bio.delta %>% 
    filter(treatment == "T") %>%
    data.frame() -> top

multi.collinear(top[, c(5:7)])

# testing biodiversity model by YEAR

rf.bio.imp.top <- rf.modelSel(top[,c(5:7)], as.factor(top$disturbance_severity), 
                              seed = 666, 
                              imp.scale="se",
                              final.model = TRUE) 

x11()
plot(rf.bio.imp.top, imp = "sel")
rf.bio.top <- randomForest(as.factor(disturbance_severity) ~ drichness + dsw.index, 
                           data = top, 
                           ntree = 2000,
                           importance = TRUE)

# find significance
rf.bio.top
rf.significance(rf.bio.top, top[, c(5:6)], ntree = 2000, nperm = 999)

### 
x11()
ggplot(b.85, aes(x = treatment, y = drichness))+
    geom_boxplot()

x11()
ggplot(b.85, aes(x = treatment, y = dsw.index))+
    geom_boxplot()







#### ADDING IN DISTURBANCE
biodiversity <- merge(biodiversity, dis.meta.data)
biodiversity$txt <- as.factor(paste(biodiversity$disturbance_severity, biodiversity$treatment, sep = ""))

# filter to 2020
biodiversity %>%
    filter(year == 2020) -> bio.2020

rf.bio.imp.SEVERITY <- rf.modelSel(bio.2020[,3:5], as.factor(bio.2020$disturbance_severity),
                                   seed = 666, 
                                   imp.scale="se",
                                   final.model = TRUE) 

rf.bio.imp.TREAT <- rf.modelSel(bio.2020[,3:5], as.factor(bio.2020$treatment),
                                seed = 666, 
                                imp.scale="mir",
                                final.model = TRUE)

rf.bio.imp.TXT <- rf.modelSel(bio.2020[,3:5], as.factor(bio.2020$txt),
                              seed = 666, 
                              imp.scale="mir",
                              final.model = TRUE)

rf.TXT <- randomForest(txt ~  sw.index + richness,
                       data = bio.2020,
                       ntree = 2000,
                       mtry = 2,
                       importance = TRUE)

rf.significance(rf.TXT, biodiversity[, c(3:5)], q = 0.99, p = 0.05, nperm = 999)

x11()
plot(rf.bio.imp.SEVERITY, imp = "sel")

# short modelrich
rf.bio.SEV <- randomForest(as.factor(disturbance_severity) ~  evenness + sw.index + richness,
                           data = bio.2020,
                           ntree = 2000,
                           mtry = 3,
                           importance = TRUE)

rf.bio.SEV
caret::varImp(rf.bio.SEV)




