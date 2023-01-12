# FoRTE Structural Analysis 2021 
# Jeff W. Atkins (jeffrey.atkins@usda.gov)
# 
# 01_data_import_and_clenaing_script.R

# dependencies
require(fortedata)
require(tidyverse)
require(randomForest)
require(ggridges)
require(ggplot2)
require(vegan)
require(tidyverse)
require(randomForest)
require(rfUtilities)
require(caret)
require(car)
require(Hmisc)
require(reldist)
require(corrplot)
require(viridis)



####### Disturbance metadata
# bring in metadata
df <- fortedata::fd_plot_metadata()

df <- data.frame(df)
df$subplot_id <- paste(df$replicate, 0, df$plot, df$subplot, sep = "")
df$subplot_id <- as.factor(df$subplot_id)

# sorts disturbance metadata
df %>%
    dplyr::select(subplot_id, disturbance_severity, treatment) %>%
    distinct() %>%
    data.frame() %>%
    na.omit() -> dis.meta.data



# ###### Lidar derived forest canopy structural data
# 
# # forest structure from lidar
csc <- fd_canopy_structure()

# makes columsn into factors for analysis
csc$replicate <- as.factor(csc$replicate)
csc$year <- as.factor(csc$year)
csc$subplot_id <- as.factor(csc$subplot_id)


# fix CSC to plot level

csc %>%
    dplyr::group_by(year, subplot_id) %>%
    dplyr::summarise_at(vars(mean_height_mean:p90),~round(mean(., na.rm = TRUE), 2)) %>%
    data.frame() -> csc

# merges csc w/ disturbance metadata
csc <- merge(csc, dis.meta.data)


csc %>%
    dplyr::select(subplot_id, disturbance_severity, treatment, year, rugosity, top_rugosity, enl, fhd, porosity,
                  clumping_index, mean_height_mean, moch, max_ht, p10, p25, p50, p75, p90, cover_fraction, vai_mean) %>%
    dplyr::filter(year != 2019) %>%
    dplyr::group_by(subplot_id, disturbance_severity, treatment, year) %>%
    dplyr::summarise_at(vars(rugosity:vai_mean),~round(mean(., na.rm = TRUE), 2)) %>%
    pivot_longer( cols = rugosity:vai_mean) %>%
    pivot_wider(names_from = year, names_glue = "mean{year}") %>%
    data.frame() -> big.boi.means


# csc %>%
#     dplyr::select(disturbance_severity, year, rugosity, top_rugosity, enl, fhd, porosity,
#                   clumping_index) %>%
#     dplyr::filter(year != 2019) %>%
#     dplyr::group_by(disturbance_severity, year) %>%
#     dplyr::summarise_at(vars(rugosity:clumping_index),~round(sd(., na.rm = TRUE), 2)) %>%
#     pivot_wider(names_from = year, values_from = c(rugosity:clumping_index)) %>%

# make all the csc

# csc treatments
big.boi.means$diffMean <- (big.boi.means$mean2020 - big.boi.means$mean2018) / big.boi.means$mean2018
big.boi.means$logRR <- log(big.boi.means$mean2018 / big.boi.means$mean2020)
#df.csc <- merge(big.boi.means, big.boi.sd)
df.csc <- big.boi.means
df.csc$s = sqrt( ((df.csc$mean2018^2) + (df.csc$mean2020^2))/2 )
df.csc$cohens.d = (df.csc$mean2018 - df.csc$mean2020) / df.csc$s


# #### code effect size
# df.csc$effect_Size <- if (df.csc$cohens.d < 0.19){
#     effect_Size = "Very Small"
# } else if (df.csc$cohens.d >= 0.2 && df.csc$cohens.d < 0.5){
#     effect_Size = "Small"
# } else if (df.csc$cohens.d >= 0.5 && df.csc$cohens.d < 0.8){
#     effect_Size = "Medium"
# } else {effect_size = "crazy"}









##### STAND STRUCTURE Data

## bring in the inventory and mortalitiy data sets
mor <- fd_mortality()
#inv <- fd_inventory()

# make factors
mor$subplot_id <- as.factor(mor$subplot_id)

# making the PFT for mor
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

# make basal area per plot by kill and live
mor %>%
    dplyr::group_by(subplot_id, fate) %>%
    dplyr::filter(health_status == "L") %>%
    dplyr::summarize(ba_plot = sum(basal_area)) -> subplot.ba

# basal area per species per kill,live per plot
mor %>%
    group_by(subplot_id, species, stage, fate) %>%
    dplyr::filter(health_status == "L") %>%
    dplyr::summarize(ba_species = sum(basal_area)) -> species.ba

# counts no. of stems killed 
mor %>%
    group_by(subplot_id, fate) %>%
    dplyr::filter(health_status == "L") %>%
    dplyr::summarize(no_stems = n()) -> subplot.stems.2020


# make dbh _class
mor$dbh_class <- round(mor$dbh_cm/5) * 5
mor$ba <- ((mor$dbh_cm / 2) ^2) * pi

# gets the count of stems per species
mor %>%
    filter(health_status == "L") %>%
    filter(fate == "live") %>%
    group_by(subplot_id) %>%
    dplyr::summarize(dbh.sd = sd(dbh_cm, na.rm = TRUE),
                     ba = sum(ba, na.rm = TRUE)) %>%
    data.frame() -> dbh.live

mor %>%
    filter(health_status == "L") %>%
    group_by(subplot_id) %>%
    dplyr::summarize(dbh.sd = sd(dbh_cm, na.rm = TRUE),
                     ba = sum(ba, na.rm = TRUE)) %>%
    data.frame() -> dbh.total

####### count how many of each
mor %>%
    filter(health_status == "L") %>%
    filter(fate == "live") -> live


# Make tree size diversity index 
dbh.classes <- data.frame(table(live$dbh_class, live$subplot_id))
a <- spread(dbh.classes, Var1, Freq)

b <- diversity(a[, c(2:12)], index = "simpson")

##### dbh diversity
dbh.div <- data.frame(a$Var2)
dbh.div$dbh_shan <- b

colnames(dbh.div)[1] <- "subplot_id"


##### Gini inefficiency 

giniplots <- aggregate(dbh_cm ~ subplot_id,
                       data = live,
                       FUN = "gini")

names(giniplots) <- c("subplot_id", "gini")

df <- merge(dbh.live, giniplots)
df <- merge(df, dbh.div)

df.post <- merge(df, dis.meta.data)   


###### now for pre 
dbh.classes <- data.frame(table(mor$dbh_class, mor$subplot_id))
a <- spread(dbh.classes, Var1, Freq)

b <- diversity(a[, c(2:12)], index = "simpson")

##### dbh diversity
dbh.div <- data.frame(a$Var2)
dbh.div$dbh_shan <- b

colnames(dbh.div)[1] <- "subplot_id"


#####

giniplots <- aggregate(dbh_cm ~ subplot_id,
                       data = mor,
                       FUN = "gini")

names(giniplots) <- c("subplot_id", "gini")

df <- merge(dbh.total, giniplots)
df <- merge(df, dbh.div)

df.pre <- merge(df, dis.meta.data)   

df.pre$year <- as.factor(2018)
df.post$year <- as.factor(2020)

df <- rbind(df.pre, df.post)










#### MAKE WIDE AND CALCULATE SCALAR DIFFERENCES
df %>%
    #dplyr::select(disturbance_severity, treatment, year, mean_height_mean, moch, max_ht, p10, p25, p50, p75, p90, cover_fraction, vai_mean) %>%
    #dplyr::filter(year != 2019) %>%
    dplyr::group_by(disturbance_severity, treatment, year) %>%
    #dplyr::summarise_at(vars(mean_height_mean:vai_mean),~round(mean(., na.rm = TRUE), 2)) %>%
    pivot_longer( cols = dbh.sd:dbh_shan    ) %>%
    pivot_wider(names_from = year, names_glue = "mean{year}") %>%
    data.frame() -> stand.struct.means



stand.struct.means$diffMean <- (stand.struct.means$mean2020 - stand.struct.means$mean2018) / stand.struct.means$mean2018
stand.struct.means$logRR <- log(stand.struct.means$mean2018 / stand.struct.means$mean2020)
df.ss <- stand.struct.means
df.ss$s = sqrt( ((df.ss$mean2018^2) + (df.ss$mean2020^2))/2 )

df.ss$cohens.d = (df.ss$mean2018 - df.ss$mean2020) / df.ss$s

x11()
ggplot(df.ss, aes(x = diffMean, y = logRR))+
    geom_point()
#### code effect size
# df.ss$effect_size <- if (df.ss$cohens.d < 0.2) {
#     effect_size = "Very Small"
# } else if (df.ss$cohens.d >= 0.2 && df.ss$cohens.d < 0.5){
#     effect_size = "Small"
# } else if (df.ss$cohens.d >= 0.5 && df.ss$cohens.d < 0.8){
#     effect_size = "Medium"
# } else {effect_size = "crazy"}


#write.csv(df.csc, "annual_differences_stats.csv")














#######################
# Community Analyses


# gets the count of stems per species
mor %>%
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

mor %>%
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
biodiversity <- merge(biodiversity, dis.meta.data)   


# bio wide
biodiversity %>%
    dplyr::group_by(disturbance_severity, treatment, year) %>%
    #dplyr::summarise_at(vars(mean_height_mean:vai_mean),~round(mean(., na.rm = TRUE), 2)) %>%
    pivot_longer( cols = richness:evenness) %>%
    pivot_wider(names_from = year, names_glue = "mean{year}") %>%
    data.frame() -> df.bio






df.bio$diffMean <- (df.bio$mean2020 - df.bio$mean2018) / df.bio$mean2018
df.bio$logRR <- log(df.bio$mean2018 / df.bio$mean2020)

df.bio$s = sqrt( ((df.bio$mean2018^2) + (df.bio$mean2020^2))/2 )

df.bio$cohens.d = (df.bio$mean2018 - df.bio$mean2020) / df.bio$s

df.bio$txt <- NULL

# bind
df.all <- rbind(df.bio, df.ss, df.csc)


#### writing the data
# #write.csv(biodiversity, "./data/forte_biodiversity_2018_2020.csv")
# 
# # bringin meta
# bio.2020 <- merge(bio.2020, dis.meta.data)
# bio.2020$txt <- as.factor(paste(bio.2020$disturbance_severity, bio.2020$treatment, sep = ""))
# 
# bio.2020$evenness_2020[bio.2020$evenness_2020 == 0] <- 0.0001
# 
# 
# ###### calculate difference from rep level control
# bio.2020$replicate <- substr(bio.2020$subplot_id, 0, 1)
# 
# # find control values
# bio.2020 %>%
#     filter(disturbance_severity == 0) -> controls
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # modeling
# summary(lm(dbh.sd ~ disturbance_severity, data = dbh))
# summary(lm(gini ~ disturbance_severity, data = df.post))
# summary(lm(dbh_shan ~ disturbance_severity, data = df.post))
# 
# summary(lm(dbh.sd ~ treatment, data = dbh))
# summary(lm(gini ~ treatment, data = df.post))
# summary(lm(dbh_shan ~ treatment, data = df.post))
# 
# #####
# mor %>%
#     filter(fate == "live") %>%
#     dplyr::group_by(subplot_id, pft) %>%
#     dplyr::summarize(ba = sum(basal_area, na.rm = TRUE)) %>%
#     data.frame() -> pft
# 
# # merge
# pft <- merge(pft, dis.meta.data)
# 
# summary(lm(ba ~ disturbance_severity, data = pft))
# 
# #
# mor %>%
#     count()
# mtcars %>%
#     count(am, gear) %>%
#     group_by(am) %>%
#     mutate(freq = n / sum(n))
# 
# x11()
# ggplot(pft, aes(x = as.factor(disturbance_severity), y = ba))+
#     geom_boxplot()