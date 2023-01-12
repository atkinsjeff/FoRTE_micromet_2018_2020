# condensed test script for jeff evans
# require(fortedata)
# require(lubridate)
require(tidyverse)
# require(Boruta)
# require(Metrics)
# require(carat)
# require(corrplot)
# require(MuMIn)
# require(relaimpo)
# require(gridExtra)
# require(cowplot)
# require(vegan)
require(rfUtilities)
require(randomForest)


# read in data
df <- read.csv("atkins_test_data.csv")

# filter to top and bottom treatments
df %>%
    dplyr::filter(treatment == "B") -> df.bottom

df %>%
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
# full variable list
c("richness","sw.index","evenness","dbh.sd","ba","gini","dbh_shan","rugosity","top_rugosity",
  "enl","fhd","porosity","clumping_index","mean_height_mean","moch","max_ht", "p10", "p25",
  "p50", "p75","p90","cover_fraction")

##### TEMP
# bottom-up
# multi colinearity test
df.bottom <- na.omit(df.bottom)
multi.collinear(df.bottom[ , c("richness","sw.index","evenness","dbh.sd","ba","gini","dbh_shan","rugosity","top_rugosity",
                               "enl","fhd","porosity","clumping_index","mean_height_mean","moch","max_ht", "p10", "p25",
                               "p50", "p75","p90","cover_fraction")], p=0.05,  leave.out = TRUE, n = 999)
# 
# multi.collinear(df.bottom[ ,c("richness","sw.index","evenness","dbh.sd","ba","gini","dbh_shan","rugosity","top_rugosity",
#                               "enl","fhd","porosity","clumping_index","mean_height_mean","p90","cover_fraction")], p=0.05, perm = TRUE,  leave.out = TRUE, n = 999)


#To improve model fit we will test for multicolinearity using the “multi.collinear” 
#function in the rfUtilities library. If any variables are identified as multicolinear we will remove them from the model.
df.bottom.test <- df.bottom[ ,c("richness","sw.index","evenness","dbh.sd","ba","gini","dbh_shan","rugosity","top_rugosity",
                                "enl","fhd","porosity","clumping_index","mean_height_mean","moch","max_ht", "p10", "p25",
                                "p50", "p75","p90","cover_fraction")]
cl <- multi.collinear(df.bottom.test, p=0.05)



#We observe that four variables: "mat30", "sar", "spost27", and "spost9" are multicollinear. There are however, variables that can act as "hingepin" variables and cause other variables to appear multicolnear whent there are in fact, not redundant. Here we can perform a "leave one out" test to evaluate if any of the varialbes are forcing other varibles to appear multicollinear. 



for(l in cl) {
    
    cl.test <-df.bottom.test[,-which(names(df.bottom.test)==l)]
    
    print(paste("Remove variable", l, sep=": "))
    
    multi.collinear(cl.test, p=0.05) 
    
}

set.seed(666)
rf.bot.temp <- rf.modelSel(df.bottom[ , c("richness","sw.index","evenness","dbh.sd","ba","gini","dbh_shan","rugosity","top_rugosity",
                                          "enl","fhd","porosity","clumping_index","mean_height_mean", "cover_fraction")], 
                           df.bottom[,"temp"], seed = 666, r = c(0.25, 0.5, 0.75), imp.scale = "mir", final.model = TRUE, ntree = 5000, parsimony = 0.1, mtry = 4) 
rf.bot.temp <- rf.modelSel(df.bottom[ ,c( "rugosity", "top_rugosity", "clumping_index", "cover_fraction" )], 
                           df.bottom[,"temp"], seed = 666, r = c(0.25, 0.5, 0.75), imp.scale = "mir", final.model = TRUE, ntree = 5000, parsimony = 0.1 ) 
x11()
plot(rf.bot.temp) # plot importance for selected variables

rf.bot.temp$rf.final

# correlation
cor(rf.bot.temp$rf.final$y, rf.bot.temp$rf.final$predicted)

rf.bot.temp$selvars
rf.regression.fit(rf.bot.temp$rf.final)
rf.significance(rf.bot.temp$rf.final, df.bottom[ , c("rugosity", "top_rugosity", "cover_fraction")], p = 0.05, nperm = 999, seed = 666)


### This shows how much our model accuracy decreases if we leave out that variable. 
### Mean Decrease Gini (IncNodePurity) - This is a measure of variable importance based on the Gini 
### impurity index used for the calculating the splits in trees.
randomForest::importance(rf.bot.temp$rf.final, scale = FALSE)
randomForest::importance(rf.bot.temp$rf.final)
# rf.bot.temp$sel.importance
# rf.bot.temp$importance

#
rf.effectSize(rf.bot.temp$rf.final, y = df.bottom$rugosity, pred.data = df.bottom[ ,c("rugosity", "top_rugosity","cover_fraction")], x.var = rugosity )
rf.effectSize(rf.bot.temp$rf.final, y = df.bottom$cover_fraction, pred.data = df.bottom[ ,c("rugosity", "top_rugosity","cover_fraction")], x.var = cover_fraction )
rf.effectSize(rf.bot.temp$rf.final, y = df.bottom$top_rugosity , pred.data = df.bottom[ ,c("rugosity", "top_rugosity","cover_fraction")], x.var = top_rugosity  )


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
    rf.partial.ci(m = rf.bot.temp$rf.final, x = df.bottom[ , c("rugosity", "top_rugosity","cover_fraction")], yname = "", xname = i, delta = TRUE, lci = 0.025, uci = 0.975)
}

x11()
plot(rf.bot.temp)

rf.cv <- rf.crossValidation(rf.bot.temp$rf.final, df.bottom[ ,c("rugosity", "top_rugosity","cover_fraction")], 
                            p=0.05, n=99, ntree=501) 

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
p <- as.matrix(rf.bot.temp$rf.final$importance[,2])   

ord <- rev(order(p[,1], decreasing=TRUE)[1:dim(p)[1]]) 

dotchart(p[ord,1], main="Scaled Variable Importance", pch=19)  









# top
multi.collinear(df.top[ , c("richness","sw.index","evenness","dbh.sd","ba","gini","dbh_shan","rugosity","top_rugosity",
                            "enl","fhd","porosity","clumping_index","mean_height_mean","moch","max_ht", "p10", "p25",
                            "p50", "p75","p90","cover_fraction")], p=0.05,  leave.out = TRUE, n = 999)

rf.top.temp <- rf.modelSel(df.top[ , c("richness","sw.index","evenness","dbh.sd","ba","gini","dbh_shan","rugosity","top_rugosity",
                                       "enl","fhd","porosity","clumping_index","mean_height_mean","moch","cover_fraction")], 
                           df.top[,"temp"], seed = 666, imp.scale = "se", final.model = TRUE,  ntree = 2000, parsimony = 0.01 ) 

x11()
plot(rf.top.temp) # plot importance for selected variables

rf.bot.temp$rf.final



# correlation
cor(rf.top.temp$rf.final$y, rf.top.temp$rf.final$predicted)

rf.top.temp$selvars
rf.regression.fit(rf.top.temp$rf.final)
rf.significance(rf.top.temp$rf.final, df.top[ ,4:13], p = 0.05, nperm = 999, seed = 666)

# selected variables
rf.top.temp$selvars

rf.effectSize(rf.top.temp$rf.final, y = df.top$dmoch, pred.data = df.top[ ,4:13], x.var = dmoch )
rf.effectSize(rf.top.temp$rf.final, y = df.top$dcc, pred.data = df.top[ ,4:13], x.var = dcc )
