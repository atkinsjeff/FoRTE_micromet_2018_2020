# comparing pre and post disturbance
require(fortedata)
require(lubridate)
require(tidyverse)
require(Boruta)
require(Metrics)
require(carat)
require(corrplot)
require(MuMIn)
require(relaimpo)
require(gridExtra)
require(cowplot)
require(vegan)
require(rfUtilities)

# bring in metadata
df <- fortedata::fd_plot_metadata()

df <- data.frame(df)
df$subplot_id <- paste(df$replicate, 0, df$plot, df$subplot, sep = "")
df$subplot_id <- as.factor(df$subplot_id)

df %>%
    dplyr::select(subplot_id, disturbance_severity, treatment) %>%
    distinct() %>%
    data.frame() %>%
    na.omit() -> dis.meta.data


##### IMPORT DATA FRAME FROM 01_data_import
df.18 <- df.all[, c("subplot_id","disturbance_severity", "treatment", "name", "mean2018")]

# pivot wider
df.18 %>%
    pivot_wider(names_from = "name",
                values_from = c("mean2018")) %>%
    data.frame() -> df.18

# remove NA
df.18 <- na.omit(df.18)

# bring in met data
met <- data.frame(fortedata::fd_soil_respiration())
met$year <- as.factor(format(as.Date(met$date, format = "%Y-%m-%d"), "%Y"))


# organize 
met %>%
    # select(subplot_id, date, soil_temp, vwc) %>%
    group_by(subplot_id) %>%
    dplyr::filter(between(date, as.Date("2018-06-01"), as.Date("2018-08-31"))) %>%
    dplyr::summarize(temp = mean(soil_temp, na.rm = TRUE),
                     vwc = mean(vwc, na.rm = TRUE)) %>%
    data.frame() -> x

met %>%
    # select(subplot_id, date, soil_temp, vwc) %>%
    group_by(subplot_id) %>%
    dplyr::filter(between(date, as.Date("2018-06-01"), as.Date("2018-08-31"))) %>%
    dplyr::summarize(temp.sd = sd(soil_temp, na.rm = TRUE),
                     vwc.sd = sd(vwc, na.rm = TRUE)) %>%
    data.frame() -> y
# bring em together
x <- merge(x, y)
df.18 <- merge(df.18, x)


# light data
par <- data.frame(fortedata::fd_ceptometer())


# quick par test
par$TEST <- ifelse(is.na(par$notes), 1, 0)

x11()
ggplot(par, aes(x = fapar, color = as.factor(TEST)))+
    geom_boxplot()

par$year <- year(par$timestamp)
####
par %>%
    #filter(between(timestamp, as.Date("2020-06-01"), as.Date("2020-08-31"))) %>%
    filter(year == 2020) %>%
    group_by(subplot_id) %>%
    dplyr::summarize(par = mean(fapar, na.rm = TRUE)) %>%
    data.frame() -> z


#### all together now
df.18 <- merge(df.18, z)



### SORT TO TOP AND BOTTOM

# filter to top and bottom
df.18 %>%
    dplyr::filter(treatment == "B") -> df.bottom

df.18 %>%
    dplyr::filter(treatment == "T") -> df.top



### CORRELATION
require(corrplot)
m <- cor(df.bottom[4:20], method = "pearson", use = "complete.obs")
x11()
corrplot(m, method = "number", type = "lower")
m.bot <- rcorr(as.matrix(df.bottom[4:20]))
m.top <- rcorr(as.matrix(df.top[4:20]))

x11()
corrplot(m.bot[[1]], 
         method = "number",
         #order = 'AOE', 
         addCoef.col = 'black', 
         cl.pos = 'n',
         #type = "lower",    # Correlation plot style (also "upper" and "lower")
         diag = FALSE,
         p.mat = m.bot[[3]], 
         sig.level = 0.05, # displays only corr. coeff. for p < 0.05
         insig = "blank",  # else leave the cell blank
         tl.srt = 0,       # control the orintation of text labels
         tl.offset = 1) 

x11()
corrplot(m.top[[1]], 
         method = "number",
         #order = 'AOE', 
         addCoef.col = 'black', 
         cl.pos = 'n',
         #type = "upper",    # Correlation plot style (also "upper" and "lower")
         diag = FALSE,
         p.mat = m.top[[3]], 
         sig.level = 0.05, # displays only corr. coeff. for p < 0.05
         insig = "blank",  # else leave the cell blank
         tl.srt = 0,       # control the orintation of text labels
         tl.offset = 1) 


###### RANDOM FOREST MODEL
col.names <- list(colnames(df.18))
print(col.names)

c("richness","sw.index","evenness","dbh.sd", "ba", "gini","dbh_shan",
  "mean_height_mean","moch","max_ht", "p10", "p25", "p50", "p75", "p90",                 
  "cover_fraction", "vai_mean")
##### TEMP
# bottom-up
# multi colinearity test
multi.collinear(df.bottom[ ,4:20], p=0.05)

set.seed(666)
rf.bot.temp <- rf.modelSel(df.bottom[ , c("richness","sw.index","evenness","dbh.sd", "ba", 
                                          "mean_height_mean","moch","max_ht", "p10", "p25", "p90")], 
                           df.bottom[,"temp"], seed = 666, imp.scale = "se", final.model = TRUE,  ntree = 2000, parsimony = 0.05 ) 
rf.bot.temp$rf.final

# correlation
cor(rf.bot.temp$rf.final$y, rf.bot.temp$rf.final$predicted)

rf.bot.temp$selvars
rf.regression.fit(rf.bot.temp$rf.final)
rf.significance(rf.bot.temp$rf.final, df.bottom[ ,4:13], p = 0.05, nperm = 999, seed = 666)


### This shows how much our model accuracy decreases if we leave out that variable. 
### Mean Decrease Gini (IncNodePurity) - This is a measure of variable importance based on the Gini 
### impurity index used for the calculating the splits in trees.
randomForest::importance(rf.bot.temp$rf.final, scale = FALSE)
randomForest::importance(rf.bot.temp$rf.final)
# rf.bot.temp$sel.importance
# rf.bot.temp$importance

#
rf.effectSize(rf.bot.temp$rf.final, y = df.bottom$drug, pred.data = df.bottom[ ,4:13], x.var = drug )
rf.effectSize(rf.bot.temp$rf.final, y = df.bottom$dcc, pred.data = df.bottom[ ,4:13], x.var = dcc )
rf.effectSize(rf.bot.temp$rf.final, y = df.bottom$ddbh_shan , pred.data = df.bottom[ ,4:13], x.var = ddbh_shan  )


rf.bot.temp$importance
rf.bot.temp$rf.final$importance

p <- as.matrix(rf.bot.temp$rf.final$importance)   

ord <- rev(order(p[,1], decreasing=TRUE)[1:dim(p)[1]]) 

x11()
dotchart(p[ord,1], main="Scaled Variable Importance", pch=19)

### CROSS VALIDATION
x11()
par(mfrow=c(2,2))
for(i in rf.bot.temp$selvars){
    rf.partial.ci(m = rf.bot.temp$rf.final, x = df.bottom[ ,4:13], yname = "", xname = i, delta = TRUE, lci = 0.25, uci = 0.75)
}

x11()
plot(rf.bot.temp)

rf.cv <- rf.crossValidation(rf.bot.temp$rf.final, df.bottom[ ,4:13], 
                            p=0.10, n=99, ntree=501) 

x11()
par(mfrow=c(2,2))
plot(rf.cv)  
plot(rf.cv, stat = "mse")
plot(rf.cv, stat = "var.exp")
plot(rf.cv, stat = "mae")

# Convergence
x11()
plot(rf.bot.temp$rf.final)


# variable importance
p <- as.matrix(rf.bot.temp$rf.final$importance[,3])   

ord <- rev(order(p[,1], decreasing=TRUE)[1:dim(p)[1]]) 

dotchart(p[ord,1], main="Scaled Variable Importance", pch=19)  









# top
rf.top.temp <- rf.modelSel(df.top[ ,4:12], df.top[,"temp"], seed = 666, imp.scale = "se", final.model = TRUE,  ntree = 2000, parsimony = 0.01 ) 

# correlation
cor(rf.top.temp$rf.final$y, rf.top.temp$rf.final$predicted)

rf.top.temp$selvars
rf.regression.fit(rf.top.temp$rf.final)
rf.significance(rf.top.temp$rf.final, df.top[ ,4:13], p = 0.05, nperm = 999, seed = 666)

# selected variables
rf.top.temp$selvars

rf.effectSize(rf.top.temp$rf.final, y = df.top$dmoch, pred.data = df.top[ ,4:13], x.var = dmoch )
rf.effectSize(rf.top.temp$rf.final, y = df.top$dcc, pred.data = df.top[ ,4:13], x.var = dcc )





##### TEMP.SD
# bottom-up
rf.bot.temp.sd <- rf.modelSel(df.bottom[ ,4:13], df.bottom[,"temp.sd"], seed = 666, imp.scale = "se", final.model = TRUE,  ntree = 2000, parsimony = 0.1 ) 

# correlation
cor(rf.bot.temp.sd$rf.final$y, rf.bot.temp.sd$rf.final$predicted)

# print vars
rf.bot.temp.sd$selvars
rf.regression.fit(rf.bot.temp.sd$rf.final)
rf.significance(rf.bot.temp.sd$rf.final, df.bottom[ ,4:12], p = 0.05, nperm = 999, seed = 666)


# top
rf.top.temp.sd <- rf.modelSel(df.top[ ,4:12], df.top[,"temp.sd"], seed = 666, imp.scale = "se", final.model = TRUE,  ntree = 2000, parsimony = 0.01 ) 

cor(rf.top.temp.sd$rf.final$y, rf.top.temp.sd$rf.final$predicted)

rf.top.temp.sd$selvars
rf.regression.fit(rf.top.temp.sd$rf.final)
rf.significance(rf.top.temp.sd$rf.final, df.top[ ,4:12], p = 0.05, nperm = 999, seed = 666)





##### VWC
# bottom-up
rf.bot.vwc <- rf.modelSel(df.bottom[ ,4:12], df.bottom[,"vwc"], seed = 666, imp.scale = "se", final.model = TRUE,  ntree = 2000 ) 

cor(rf.bot.vwc$rf.final$y, rf.bot.vwc$rf.final$predicted)

rf.bot.vwc$selvars
rf.regression.fit(rf.bot.vwc$rf.final)
rf.significance(rf.bot.vwc$rf.final, df.bottom[ ,4:13], p = 0.05, nperm = 999, seed = 666)



# top
rf.top.vwc <- rf.modelSel(df.top[ ,4:12], df.top[,"vwc"], seed = 666, imp.scale = "se", final.model = TRUE, ntree = 2000 ) 

cor(rf.top.vwc$rf.final$y, rf.top.vwc$rf.final$predicted)

rf.top.vwc$selvars
rf.regression.fit(rf.top.vwc$rf.final)
rf.significance(rf.top.vwc$rf.final, df.top[ ,4:13], p = 0.05, nperm = 999, seed = 666)


##### VWC
# bottom-up
rf.bot.vwc.sd <- rf.modelSel(df.bottom[ ,4:12], df.bottom[,"vwc.sd"], seed = 666, imp.scale = "se", final.model = TRUE,  ntree = 2000 ) 
cor(rf.bot.vwc.sd$rf.final$y, rf.bot.vwc.sd$rf.final$predicted)

rf.bot.vwc.sd$selvars
rf.regression.fit(rf.bot.vwc.sd$rf.final)
rf.significance(rf.bot.vwc.sd$rf.final, df.bottom[ ,4:13], p = 0.05, nperm = 999, seed = 666)



# top
rf.top.vwc.sd <- rf.modelSel(df.top[ ,4:12], df.top[,"vwc.sd"], seed = 666, imp.scale = "se", final.model = TRUE,  ntree = 2000 ) 
cor(rf.top.vwc.sd$rf.final$y, rf.top.vwc.sd$rf.final$predicted)

rf.top.vwc.sd$selvars
rf.regression.fit(rf.top.vwc.sd$rf.final)
rf.significance(rf.top.vwc.sd$rf.final, df.top[ ,4:13], p = 0.05, nperm = 999, seed = 666)




##### VWC
# bottom-up
rf.bot.par<- rf.modelSel(df.bottom[ ,4:12], df.bottom[,"par"], seed = 666, imp.scale = "se", final.model = TRUE,  ntree = 2000 ) 
cor(rf.bot.par$rf.final$y, rf.bot.par$rf.final$predicted)

rf.bot.par$selvars
rf.regression.fit(rf.bot.par$rf.final)
rf.significance(rf.bot.par$rf.final, df.bottom[ ,4:13], p = 0.05, nperm = 999, seed = 666)



# top
rf.top.par<- rf.modelSel(df.top[ ,4:12], df.top[,"par"], seed = 666, imp.scale = "se", final.model = TRUE,  ntree = 2000 ) 
cor(rf.top.par$rf.final$y, rf.top.par$rf.final$predicted)

rf.top.par$selvars
rf.regression.fit(rf.top.par$rf.final)
rf.significance(rf.top.par$rf.final, df.top[ ,4:13], p = 0.05, nperm = 999, seed = 666)









# Classification on iris data
require(randomForest)
data(iris)
iris$Species <- as.factor(iris$Species)
( rf.class <- rf.modelSel(iris[,1:4], iris[,"Species"], seed=1234, imp.scale="mir") )
( rf.class <- rf.modelSel(iris[,1:4], iris[,"Species"], seed=1234, imp.scale="mir", 
                          parsimony=0.03) )

plot(rf.class)              # plot importance for selected variables
plot(rf.class, imp = "all") # plot importance for all variables 

vars <- rf.class$selvars
( rf.fit <- randomForest(x=iris[,vars], y=iris[,"Species"]) )




