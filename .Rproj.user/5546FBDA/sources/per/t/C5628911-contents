# packages
invisible(lapply(c("randomForest", "rfUtilities",    
                   "ranger", "pdp", "spatialEco"),require, 
                 character.only=TRUE))

# set source directory

source("./R/rf.ImpScale.R")

pfun <- function(object, newdata, i = 2) {
    predict(object, data = newdata)$predictions[,i]
}
rmse <- function(y,x) { sqrt(mean((y - x)^2)) }

set.seed(666)

#***********************************************
# read and filter top/bottom treatments
dat <- read.csv("atkins_test_data.csv")
# Evaluate colinearity and multicolinearity
# full variable list
all.vars <- c("richness","sw.index","evenness","dbh.sd","ba","gini","dbh_shan",
              "rugosity","top_rugosity", "enl","fhd","porosity","clumping_index",
              "mean_height_mean","moch","max_ht", "p10", "p25",
              "p50", "p75","p90","cover_fraction")

# remove collinearity
cl.vars <- spatialEco::collinear(dat[,all.vars], p = 0.85, 
                                 nonlinear = FALSE, p.value = 0.05)
if(length(cl.vars) > 0) {
    cat("Dropping colinear variables; ", cl.vars, "\n")  
    dat2 <- dat[,-which(names(dat) %in% cl.vars)]
    vars <- all.vars[-which(all.vars %in% cl.vars)]
} else {
    cat("No collinear variable found", "\n")
}

# multicolinearity
( mc <- multi.collinear(dat2[,vars], p=0.005,  perm = TRUE, n = 999) )
mc.vars <- mc[which( mc$frequency > 5 ),]$variables 
if(length(mc.vars) > 0) { 
    dat2 <- dat2[,-which(names(dat2) %in% mc.vars)]
} else {
    cat("No multicollinear variable found", "\n")
}


###
df.bottom <- na.omit(dat[dat$treatment == "B",])
df.top <- na.omit(dat[dat$treatment == "T",])



#***********************************************
# # Evaluate colinearity and multicolinearity
# 
# # remove collinearity
# cl.vars <- spatialEco::collinear(df.bottom[,all.vars], p = 0.85, 
#                                  nonlinear = FALSE, p.value = 0.001)
# if(length(cl.vars) > 0) {
#     cat("Dropping colinear variables; ", cl.vars, "\n")  
#     df.bottom <- df.bottom[,-which(names(df.bottom) %in% cl.vars)]
#     vars <- all.vars[-which(all.vars %in% cl.vars)]
# } else {
#     cat("No collinear variable found", "\n")
# }
# 
# # multicolinearity
# ( mc <- multi.collinear(df.bottom[,vars], p=0.005,  perm = TRUE, n = 999) )
# mc.vars <- mc[which( mc$frequency > 5 ),]$variables 
# if(length(mc.vars) > 0) { 
#     df.bottom <- df.bottom[,-which(names(df.bottom) %in% mc.vars)]
# } else {
#     cat("No multicollinear variable found", "\n")
# }

#***********************************************
# BOTTOM DOWN
# Random Forests Model

b = 501 # number of bootstraps

# Evaluate importance significance
( imp.reg <- ranger(x=df.bottom[,vars], y=df.bottom[,"temp"],
                    importance = "permutation", 
                    num.trees=b) )
( imp <- na.omit(rf.ImpScale(imp.reg, scaling="p", n=59)) )
( imp.vars <- imp[imp$importance > 0 & imp$pvalue <= 0.1,]$parameter )
# ( imp.vars <- imp[imp$importance > 0,]$parameter )

# Fit model
( bot.fit <- ranger(x=df.bottom[,imp.vars], y=df.bottom[,"temp"],
                    importance = "permutation", num.trees=b,
                    keep.inbag=TRUE) )

cat("R-square of model:", bot.fit$r.squared, "\n")
cat("Prediction RMSE:", rmse(df.bottom[,"temp"], bot.fit$predictions), "\n")
cat("Prediction correlation:", cor(df.bottom[,"temp"], bot.fit$predictions), "\n")

# effect size
rf.bot.temp <- randomForest(temp ~ rugosity + cover_fraction, data = df.bottom,
                              mtry = 2, ntree = 2000)
rf.regression.fit(rf.bot.temp)

# Plot importance
p <- imp[which(imp$parameter %in% imp.vars),]   
p <- p[order(p$importance),] 
x11()
dotchart(p$importance, labels = p$parameter, 
         main="Scaled Variable Importance", pch=19)  

#*****************************************************
# Jackknife validation
x <- df.bottom[,imp.vars]
y <- df.bottom[,"temp"] 
midx <- 1:length(y)     
error <- vector()
serror <- vector()
r.sqr <- vector() 
j=0
for(i in midx) { 
    j=j+1
    cat("Running obs ", j, " of ", length(midx), "\n")
    rfv <- ranger(x=x[-i,], y=y[-i],
                  num.trees = b, 
                  importance="permutation")
    error <- append(error, rfv$prediction.error)
    serror <- append(serror, rmse(y[-i], rfv$predictions))
    r.sqr <- append(r.sqr, rfv$r.squared)	
}
jackknife <- data.frame(y=y, error=error, 
                        rmse=serror, r.sqr=r.sqr)

cat("Median Jackknife R-square:", median(r.sqr), "\n")
cat("Median Jackknife RMSE:", median(serror), "\n")
cat("Median Jackknife Prediction Error:", median(error), "\n")

# Plot Jackknife prediction and RMSE error
x11()
plot(density(error), main="Jackknife prediction error")
abline(v=bot.fit$prediction.error, col="red")
x11()
plot(density(serror), main="Jackknife RMSE")
abline(v=rmse(df.bottom[,"temp"], bot.fit$predictions), col="red")

# Calculate and plot confidence region
se <- suppressWarnings(predict(bot.fit, df.bottom[,imp.vars],
                               type = "se", se.method = "infjack"))
scalar = (se$predictions * 1.96) / max((se$predictions * 1.96)) 
ci95 <- data.frame(
    y = df.bottom[,"temp"],
    pred = bot.fit$predictions, 
    std.err=se$se,
    lower.ci = bot.fit$predictions - scalar, 
    upper.ci = bot.fit$predictions + scalar)

# plot confidence 95% region 
sort.idx <- sort(ci95$pred, index.return = TRUE)$ix
x11()
plot(ci95$pred[sort.idx], type="l", xaxt="n",
     ylab="estimates", xlab="observations (n=14)", 
     main="Estimate with 95% confidence region",  
     ylim=c(min(ci95[,4:5]), max(ci95[,4:5])))
lines(ci95$upper.ci[sort.idx], col="grey", lty=2)
lines(ci95$lower.ci[sort.idx], col="grey", lty=2)

# Partial dependency and individual conditional expectation (ICE)
pd=list()
x11()
for(i in 1:length(imp.vars)){
    pd[[i]] <- partial(bot.fit, pred.var=imp.vars[i], 
                       train=df.bottom[,imp.vars],
                       main=paste0("Partial dependency of ", imp.vars[i]),
                       plot = TRUE, ice=FALSE)
}
print(pd[[1]], split = c(1, 1, 2, 2), more = TRUE)
print(pd[[2]], split = c(2, 1, 2, 2), more = TRUE)
print(pd[[3]], split = c(1, 2, 2, 2), more = TRUE)  

# save.image("C:/evans/USFS/SRS/bottom_model.RData")







#### BOTTOM UP TEMP.SD
b = 501 # number of bootstraps

# Evaluate importance significance
( imp.reg <- ranger(x=df.bottom[,vars], y=df.bottom[,"temp.sd"],
                    importance = "permutation", 
                    num.trees=b) )
( imp <- na.omit(rf.ImpScale(imp.reg, scaling="p", n=59)) )
( imp.vars <- imp[imp$importance > 0 & imp$pvalue <= 0.1,]$parameter )
# ( imp.vars <- imp[imp$importance > 0,]$parameter )

# Fit model
( bot.fit <- ranger(x=df.bottom[,imp.vars], y=df.bottom[,"temp.sd"],
                    importance = "permutation", num.trees=b,
                    keep.inbag=TRUE) )

cat("R-square of model:", bot.fit$r.squared, "\n")
cat("Prediction RMSE:", rmse(df.bottom[,"temp.sd"], bot.fit$predictions), "\n")
cat("Prediction correlation:", cor(df.bottom[,"temp.sd"], bot.fit$predictions), "\n")

# # effect size
# rf.bot.temp.sd <- randomForest(temp ~ rugosity + cover_fraction, data = df.bottom,
#                             mtry = 2, ntree = 2000)
# rf.regression.fit(rf.bot.temp)

# rugosity found to be important, but not signficant
#summary(lm(temp.sd ~ rugosity, data = df.bottom))
# Plot importance
p <- imp[which(imp$parameter %in% imp.vars),]   
p <- p[order(p$importance),] 
x11()
dotchart(p$importance, labels = p$parameter, 
         main="Scaled Variable Importance", pch=19)  

#*****************************************************
# Jackknife validation
x <- df.bottom[,imp.vars]
y <- df.bottom[,"temp.sd"] 
midx <- 1:length(y)     
error <- vector()
serror <- vector()
r.sqr <- vector() 
j=0
for(i in midx) { 
    j=j+1
    cat("Running obs ", j, " of ", length(midx), "\n")
    rfv <- ranger(x=x[-i,], y=y[-i],
                  num.trees = b, 
                  importance="permutation")
    error <- append(error, rfv$prediction.error)
    serror <- append(serror, rmse(y[-i], rfv$predictions))
    r.sqr <- append(r.sqr, rfv$r.squared)	
}
jackknife <- data.frame(y=y, error=error, 
                        rmse=serror, r.sqr=r.sqr)

cat("Median Jackknife R-square:", median(r.sqr), "\n")
cat("Median Jackknife RMSE:", median(serror), "\n")
cat("Median Jackknife Prediction Error:", median(error), "\n")

# Plot Jackknife prediction and RMSE error
x11()
plot(density(error), main="Jackknife prediction error")
abline(v=bot.fit$prediction.error, col="red")
x11()
plot(density(serror), main="Jackknife RMSE")
abline(v=rmse(df.bottom[,"temp"], bot.fit$predictions), col="red")

# Calculate and plot confidence region
se <- suppressWarnings(predict(bot.fit, df.bottom[,imp.vars],
                               type = "se", se.method = "infjack"))
scalar = (se$predictions * 1.96) / max((se$predictions * 1.96)) 
ci95 <- data.frame(
    y = df.bottom[,"temp"],
    pred = bot.fit$predictions, 
    std.err=se$se,
    lower.ci = bot.fit$predictions - scalar, 
    upper.ci = bot.fit$predictions + scalar)

# plot confidence 95% region 
sort.idx <- sort(ci95$pred, index.return = TRUE)$ix
x11()
plot(ci95$pred[sort.idx], type="l", xaxt="n",
     ylab="estimates", xlab="observations (n=14)", 
     main="Estimate with 95% confidence region",  
     ylim=c(min(ci95[,4:5]), max(ci95[,4:5])))
lines(ci95$upper.ci[sort.idx], col="grey", lty=2)
lines(ci95$lower.ci[sort.idx], col="grey", lty=2)

# Partial dependency and individual conditional expectation (ICE)
pd=list()
x11()
for(i in 1:length(imp.vars)){
    pd[[i]] <- partial(bot.fit, pred.var=imp.vars[i], 
                       train=df.bottom[,imp.vars],
                       main=paste0("Partial dependency of ", imp.vars[i]),
                       plot = TRUE, ice=FALSE)
}
print(pd[[1]], split = c(1, 1, 2, 2), more = TRUE)
print(pd[[2]], split = c(2, 1, 2, 2), more = TRUE)
#print(pd[[3]], split = c(1, 2, 2, 2), more = TRUE)  





#### BOTTOM UP VWC
#### BOTTOM UP TEMP.SD
b = 501 # number of bootstraps

# Evaluate importance significance
( imp.reg <- ranger(x=df.bottom[,vars], y=df.bottom[,"vwc"],
                    importance = "permutation", 
                    num.trees=b) )
( imp <- na.omit(rf.ImpScale(imp.reg, scaling="p", n=59)) )
( imp.vars <- imp[imp$importance > 0 & imp$pvalue <= 0.1,]$parameter )
# ( imp.vars <- imp[imp$importance > 0,]$parameter )

# Fit model
( bot.fit <- ranger(x=df.bottom[,imp.vars], y=df.bottom[,"vwc"],
                    importance = "permutation", num.trees=b,
                    keep.inbag=TRUE) )

cat("R-square of model:", bot.fit$r.squared, "\n")
cat("Prediction RMSE:", rmse(df.bottom[,"vwc"], bot.fit$predictions), "\n")
cat("Prediction correlation:", cor(df.bottom[,"vwc"], bot.fit$predictions), "\n")

# Plot importance
p <- imp[which(imp$parameter %in% imp.vars),]   
p <- p[order(p$importance),] 
dotchart(p$importance, labels = p$parameter, 
         main="Scaled Variable Importance", pch=19)  

#*****************************************************
# Jackknife validation
x <- df.bottom[,imp.vars]
y <- df.bottom[,"temp"] 
midx <- 1:length(y)     
error <- vector()
serror <- vector()
r.sqr <- vector() 
j=0
for(i in midx) { 
    j=j+1
    cat("Running obs ", j, " of ", length(midx), "\n")
    rfv <- ranger(x=x[-i,], y=y[-i],
                  num.trees = b, 
                  importance="permutation")
    error <- append(error, rfv$prediction.error)
    serror <- append(serror, rmse(y[-i], rfv$predictions))
    r.sqr <- append(r.sqr, rfv$r.squared)	
}
jackknife <- data.frame(y=y, error=error, 
                        rmse=serror, r.sqr=r.sqr)

cat("Median Jackknife R-square:", median(serror), "\n")
cat("Median Jackknife RMSE:", median(serror), "\n")
cat("Median Jackknife Prediction Error:", median(error), "\n")

# Plot Jackknife prediction and RMSE error
x11()
plot(density(error), main="Jackknife prediction error")
abline(v=bot.fit$prediction.error, col="red")
x11()
plot(density(serror), main="Jackknife RMSE")
abline(v=rmse(df.bottom[,"temp"], bot.fit$predictions), col="red")

# Calculate and plot confidence region
se <- suppressWarnings(predict(bot.fit, df.bottom[,imp.vars],
                               type = "se", se.method = "infjack"))
scalar = (se$predictions * 1.96) / max((se$predictions * 1.96)) 
ci95 <- data.frame(
    y = df.bottom[,"temp"],
    pred = bot.fit$predictions, 
    std.err=se$se,
    lower.ci = bot.fit$predictions - scalar, 
    upper.ci = bot.fit$predictions + scalar)

# plot confidence 95% region 
sort.idx <- sort(ci95$pred, index.return = TRUE)$ix
x11()
plot(ci95$pred[sort.idx], type="l", xaxt="n",
     ylab="estimates", xlab="observations (n=14)", 
     main="Estimate with 95% confidence region",  
     ylim=c(min(ci95[,4:5]), max(ci95[,4:5])))
lines(ci95$upper.ci[sort.idx], col="grey", lty=2)
lines(ci95$lower.ci[sort.idx], col="grey", lty=2)

# Partial dependency and individual conditional expectation (ICE)
pd=list()
for(i in 1:length(imp.vars)){
    pd[[i]] <- partial(bot.fit, pred.var=imp.vars[i], 
                       train=df.bottom[,imp.vars],
                       main=paste0("Partial dependency of ", imp.vars[i]),
                       plot = TRUE, ice=FALSE)
}
print(pd[[1]], split = c(1, 1, 2, 2), more = TRUE)
print(pd[[2]], split = c(2, 1, 2, 2), more = TRUE)
print(pd[[3]], split = c(1, 2, 2, 2), more = TRUE)  










###############################
# BOTTOM UP VWC.sd
#### BOTTOM UP TEMP.SD
b = 501 # number of bootstraps

# Evaluate importance significance
( imp.reg <- ranger(x=df.bottom[,vars], y=df.bottom[,"vwc.sd"],
                    importance = "permutation", 
                    num.trees=b) )
( imp <- na.omit(rf.ImpScale(imp.reg, scaling="p", n=100)) )
( imp.vars <- imp[imp$importance > 0 & imp$pvalue <= 0.1,]$parameter )
# ( imp.vars <- imp[imp$importance > 0,]$parameter )

# Fit model
( bot.fit <- ranger(x=df.bottom[,c("top_rugosity", "p75")], y=df.bottom[,"vwc.sd"],
                    importance = "permutation", num.trees=b,
                    keep.inbag=TRUE) )

cat("R-square of model:", bot.fit$r.squared, "\n")
cat("Prediction RMSE:", rmse(df.bottom[,"vwc.sd"], bot.fit$predictions), "\n")
cat("Prediction correlation:", cor(df.bottom[,"vwc.sd"], bot.fit$predictions), "\n")

# effect size
rf.bot.vwc.sd <- randomForest(vwc.sd ~ top_rugosity + p75 + moch , data = df.bottom,
                              mtry = 2, ntree = 2000)
rf.regression.fit(rf.bot.vwc.sd)
# Plot importance
p <- imp[which(imp$parameter %in% c("top_rugosity", "p75")),]   
p <- p[order(p$importance),] 
x11()
dotchart(p$importance, labels = p$parameter, 
         main="Scaled Variable Importance", pch=19)  

#*****************************************************
# Jackknife validation
x <- df.bottom[,c("top_rugosity", "p75")]
y <- df.bottom[,"vwc.sd"] 
midx <- 1:length(y)     
error <- vector()
serror <- vector()
r.sqr <- vector() 
j=0
for(i in midx) { 
    j=j+1
    cat("Running obs ", j, " of ", length(midx), "\n")
    rfv <- ranger(x=x[-i,], y=y[-i],
                  num.trees = b, 
                  importance="permutation")
    error <- append(error, rfv$prediction.error)
    serror <- append(serror, rmse(y[-i], rfv$predictions))
    r.sqr <- append(r.sqr, rfv$r.squared)	
}
jackknife <- data.frame(y=y, error=error, 
                        rmse=serror, r.sqr=r.sqr)

cat("Median Jackknife R-square:", median(r.sqr), "\n")
cat("Median Jackknife RMSE:", median(serror), "\n")
cat("Median Jackknife Prediction Error:", median(error), "\n")

# Plot Jackknife prediction and RMSE error
x11()
plot(density(error), main="Jackknife prediction error")
abline(v=bot.fit$prediction.error, col="red")
x11()
plot(density(serror), main="Jackknife RMSE")
abline(v=rmse(df.bottom[,"vwc"], bot.fit$predictions), col="red")

# Calculate and plot confidence region
se <- suppressWarnings(predict(bot.fit, df.bottom[,imp.vars],
                               type = "se", se.method = "infjack"))
scalar = (se$predictions * 1.96) / max((se$predictions * 1.96)) 
ci95 <- data.frame(
    y = df.bottom[,"vwc"],
    pred = bot.fit$predictions, 
    std.err=se$se,
    lower.ci = bot.fit$predictions - scalar, 
    upper.ci = bot.fit$predictions + scalar)

# plot confidence 95% region 
sort.idx <- sort(ci95$pred, index.return = TRUE)$ix
x11()
plot(ci95$pred[sort.idx], type="l", xaxt="n",
     ylab="estimates", xlab="observations (n=14)", 
     main="Estimate with 95% confidence region",  
     ylim=c(min(ci95[,4:5]), max(ci95[,4:5])))
lines(ci95$upper.ci[sort.idx], col="grey", lty=2)
lines(ci95$lower.ci[sort.idx], col="grey", lty=2)

# Partial dependency and individual conditional expectation (ICE)
pd=list()
for(i in 1:length(imp.vars)){
    pd[[i]] <- partial(bot.fit, pred.var=imp.vars[i], 
                       train=df.bottom[,imp.vars],
                       main=paste0("Partial dependency of ", imp.vars[i]),
                       plot = TRUE, ice=FALSE)
}
x11()
print(pd[[1]], split = c(1, 1, 2, 2), more = TRUE)
print(pd[[2]], split = c(2, 1, 2, 2), more = TRUE)
print(pd[[3]], split = c(1, 2, 2, 2), more = TRUE)  





##########
# Bottom-up par
#### BOTTOM UP TEMP.SD
b = 501 # number of bootstraps

# Evaluate importance significance
( imp.reg <- ranger(x=df.bottom[,vars], y=df.bottom[,"par"],
                    importance = "permutation", 
                    num.trees=b) )
( imp <- na.omit(rf.ImpScale(imp.reg, scaling="p", n=59)) )
( imp.vars <- imp[imp$importance > 0 & imp$pvalue <= 0.1,]$parameter )
# ( imp.vars <- imp[imp$importance > 0,]$parameter )

# Fit model
( bot.fit <- ranger(x=df.bottom[,imp.vars], y=df.bottom[,"par"],
                    importance = "permutation", num.trees=b,
                    keep.inbag=TRUE) )

cat("R-square of model:", bot.fit$r.squared, "\n")
cat("Prediction RMSE:", rmse(df.bottom[,"par"], bot.fit$predictions), "\n")
cat("Prediction correlation:", cor(df.bottom[,"par"], bot.fit$predictions), "\n")

# Plot importance
p <- imp[which(imp$parameter %in% imp.vars),]   
p <- p[order(p$importance),] 
dotchart(p$importance, labels = p$parameter, 
         main="Scaled Variable Importance", pch=19)  

#*****************************************************
# Jackknife validation
x <- df.bottom[,imp.vars]
y <- df.bottom[,"par"] 
midx <- 1:length(y)     
error <- vector()
serror <- vector()
r.sqr <- vector() 
j=0
for(i in midx) { 
    j=j+1
    cat("Running obs ", j, " of ", length(midx), "\n")
    rfv <- ranger(x=x[-i,], y=y[-i],
                  num.trees = b, 
                  importance="permutation")
    error <- append(error, rfv$prediction.error)
    serror <- append(serror, rmse(y[-i], rfv$predictions))
    r.sqr <- append(r.sqr, rfv$r.squared)	
}
jackknife <- data.frame(y=y, error=error, 
                        rmse=serror, r.sqr=r.sqr)

cat("Median Jackknife R-square:", median(serror), "\n")
cat("Median Jackknife RMSE:", median(serror), "\n")
cat("Median Jackknife Prediction Error:", median(error), "\n")

# Plot Jackknife prediction and RMSE error
x11()
plot(density(error), main="Jackknife prediction error")
abline(v=bot.fit$prediction.error, col="red")
x11()
plot(density(serror), main="Jackknife RMSE")
abline(v=rmse(df.bottom[,"par"], bot.fit$predictions), col="red")

# Calculate and plot confidence region
se <- suppressWarnings(predict(bot.fit, df.bottom[,imp.vars],
                               type = "se", se.method = "infjack"))
scalar = (se$predictions * 1.96) / max((se$predictions * 1.96)) 
ci95 <- data.frame(
    y = df.bottom[,"par"],
    pred = bot.fit$predictions, 
    std.err=se$se,
    lower.ci = bot.fit$predictions - scalar, 
    upper.ci = bot.fit$predictions + scalar)

# plot confidence 95% region 
sort.idx <- sort(ci95$pred, index.return = TRUE)$ix
x11()
plot(ci95$pred[sort.idx], type="l", xaxt="n",
     ylab="estimates", xlab="observations (n=14)", 
     main="Estimate with 95% confidence region",  
     ylim=c(min(ci95[,4:5]), max(ci95[,4:5])))
lines(ci95$upper.ci[sort.idx], col="grey", lty=2)
lines(ci95$lower.ci[sort.idx], col="grey", lty=2)

# Partial dependency and individual conditional expectation (ICE)
pd=list()
for(i in 1:length(imp.vars)){
    pd[[i]] <- partial(bot.fit, pred.var=imp.vars[i], 
                       train=df.bottom[,imp.vars],
                       main=paste0("Partial dependency of ", imp.vars[i]),
                       plot = TRUE, ice=FALSE)
}
print(pd[[1]], split = c(1, 1, 2, 2), more = TRUE)
print(pd[[2]], split = c(2, 1, 2, 2), more = TRUE)
print(pd[[3]], split = c(1, 2, 2, 2), more = TRUE)  








################ TOP
#***********************************************


#***********************************************
# top DOWN
# Random Forests Model

b = 501 # number of bootstraps

# Evaluate importance significance
( imp.reg <- ranger(x=df.top[,vars], y=df.top[,"temp"],
                    importance = "permutation", 
                    num.trees=b) )
( imp <- na.omit(rf.ImpScale(imp.reg, scaling="p", n=1000)) )
( imp.vars <- imp[imp$importance > 0 & imp$pvalue <= 0.1,]$parameter )
# ( imp.vars <- imp[imp$importance > 0,]$parameter )

# Fit model
( bot.fit <- ranger(x=df.top[,imp.vars], y=df.top[,"temp"],
                    importance = "permutation", num.trees=b,
                    keep.inbag=TRUE) )

cat("R-square of model:", bot.fit$r.squared, "\n")
cat("Prediction RMSE:", rmse(df.top[,"temp"], bot.fit$predictions), "\n")
cat("Prediction correlation:", cor(df.top[,"temp"], bot.fit$predictions), "\n")

imp.vars
# effect size
rf.top.temp <- randomForest(temp ~ top_rugosity + fhd + cover_fraction , data = df.top,
                              mtry = 2, ntree = 2000)
rf.regression.fit(rf.top.temp)

# Plot importance
p <- imp[which(imp$parameter %in% imp.vars),]   
p <- p[order(p$importance),] 
x11()
dotchart(p$importance, labels = p$parameter, 
         main="Scaled Variable Importance", pch=19)  

#*****************************************************
# Jackknife validation
x <- df.top[,imp.vars]
y <- df.top[,"temp"] 
midx <- 1:length(y)     
error <- vector()
serror <- vector()
r.sqr <- vector() 
j=0
for(i in midx) { 
    j=j+1
    cat("Running obs ", j, " of ", length(midx), "\n")
    rfv <- ranger(x=x[-i,], y=y[-i],
                  num.trees = b, 
                  importance="permutation")
    error <- append(error, rfv$prediction.error)
    serror <- append(serror, rmse(y[-i], rfv$predictions))
    r.sqr <- append(r.sqr, rfv$r.squared)	
}
jackknife <- data.frame(y=y, error=error, 
                        rmse=serror, r.sqr=r.sqr)

cat("Median Jackknife R-square:", median(r.sqr), "\n")
cat("Median Jackknife RMSE:", median(serror), "\n")
cat("Median Jackknife Prediction Error:", median(error), "\n")

# Plot Jackknife prediction and RMSE error
x11()
plot(density(error), main="Jackknife prediction error")
abline(v=bot.fit$prediction.error, col="red")
x11()
plot(density(serror), main="Jackknife RMSE")
abline(v=rmse(df.top[,"temp"], bot.fit$predictions), col="red")

# Calculate and plot confidence region
se <- suppressWarnings(predict(bot.fit, df.top[,imp.vars],
                               type = "se", se.method = "infjack"))
scalar = (se$predictions * 1.96) / max((se$predictions * 1.96)) 
ci95 <- data.frame(
    y = df.top[,"temp"],
    pred = bot.fit$predictions, 
    std.err=se$se,
    lower.ci = bot.fit$predictions - scalar, 
    upper.ci = bot.fit$predictions + scalar)

# plot confidence 95% region 
sort.idx <- sort(ci95$pred, index.return = TRUE)$ix
x11()
plot(ci95$pred[sort.idx], type="l", xaxt="n",
     ylab="estimates", xlab="observations (n=14)", 
     main="Estimate with 95% confidence region",  
     ylim=c(min(ci95[,4:5]), max(ci95[,4:5])))
lines(ci95$upper.ci[sort.idx], col="grey", lty=2)
lines(ci95$lower.ci[sort.idx], col="grey", lty=2)

# Partial dependency and individual conditional expectation (ICE)
pd=list()
for(i in 1:length(imp.vars)){
    pd[[i]] <- partial(bot.fit, pred.var=imp.vars[i], 
                       train=df.top[,imp.vars],
                       main=paste0("Partial dependency of ", imp.vars[i]),
                       plot = TRUE, ice=FALSE)
}
x11()
print(pd[[1]], split = c(1, 1, 2, 2), more = TRUE)
print(pd[[2]], split = c(2, 1, 2, 2), more = TRUE)
print(pd[[3]], split = c(1, 2, 2, 2), more = TRUE)  

# save.image("C:/evans/USFS/SRS/top_model.RData")







#### top UP TEMP.SD
b = 501 # number of bootstraps

# Evaluate importance significance
( imp.reg <- ranger(x=df.top[,vars], y=df.top[,"temp.sd"],
                    importance = "permutation", 
                    num.trees=b) )
( imp <- na.omit(rf.ImpScale(imp.reg, scaling="p", n=59)) )
( imp.vars <- imp[imp$importance > 0 & imp$pvalue <= 0.1,]$parameter )
# ( imp.vars <- imp[imp$importance > 0,]$parameter )

# Fit model
( bot.fit <- ranger(x=df.top[,imp.vars], y=df.top[,"temp.sd"],
                    importance = "permutation", num.trees=b,
                    keep.inbag=TRUE) )

cat("R-square of model:", bot.fit$r.squared, "\n")
cat("Prediction RMSE:", rmse(df.top[,"temp.sd"], bot.fit$predictions), "\n")
cat("Prediction correlation:", cor(df.top[,"temp.sd"], bot.fit$predictions), "\n")
# 
# lm.top.temp.sd <- lm(temp.sd ~ max_ht, data = df.top)
# summary(lm.top.temp.sd)
# Plot importance
p <- imp[which(imp$parameter %in% imp.vars),]   
p <- p[order(p$importance),] 
dotchart(p$importance, labels = p$parameter, 
         main="Scaled Variable Importance", pch=19)  

#*****************************************************
# Jackknife validation
x <- df.top[,imp.vars]
y <- df.top[,"temp"] 
midx <- 1:length(y)     
error <- vector()
serror <- vector()
r.sqr <- vector() 
j=0
for(i in midx) { 
    j=j+1
    cat("Running obs ", j, " of ", length(midx), "\n")
    rfv <- ranger(x=x[-i,], y=y[-i],
                  num.trees = b, 
                  importance="permutation")
    error <- append(error, rfv$prediction.error)
    serror <- append(serror, rmse(y[-i], rfv$predictions))
    r.sqr <- append(r.sqr, rfv$r.squared)	
}
jackknife <- data.frame(y=y, error=error, 
                        rmse=serror, r.sqr=r.sqr)

cat("Median Jackknife R-square:", median(serror), "\n")
cat("Median Jackknife RMSE:", median(serror), "\n")
cat("Median Jackknife Prediction Error:", median(error), "\n")

# Plot Jackknife prediction and RMSE error
x11()
plot(density(error), main="Jackknife prediction error")
abline(v=bot.fit$prediction.error, col="red")
x11()
plot(density(serror), main="Jackknife RMSE")
abline(v=rmse(df.top[,"temp"], bot.fit$predictions), col="red")

# Calculate and plot confidence region
se <- suppressWarnings(predict(bot.fit, df.top[,imp.vars],
                               type = "se", se.method = "infjack"))
scalar = (se$predictions * 1.96) / max((se$predictions * 1.96)) 
ci95 <- data.frame(
    y = df.top[,"temp"],
    pred = bot.fit$predictions, 
    std.err=se$se,
    lower.ci = bot.fit$predictions - scalar, 
    upper.ci = bot.fit$predictions + scalar)

# plot confidence 95% region 
sort.idx <- sort(ci95$pred, index.return = TRUE)$ix
x11()
plot(ci95$pred[sort.idx], type="l", xaxt="n",
     ylab="estimates", xlab="observations (n=14)", 
     main="Estimate with 95% confidence region",  
     ylim=c(min(ci95[,4:5]), max(ci95[,4:5])))
lines(ci95$upper.ci[sort.idx], col="grey", lty=2)
lines(ci95$lower.ci[sort.idx], col="grey", lty=2)

# Partial dependency and individual conditional expectation (ICE)
pd=list()
for(i in 1:length(imp.vars)){
    pd[[i]] <- partial(bot.fit, pred.var=imp.vars[i], 
                       train=df.top[,imp.vars],
                       main=paste0("Partial dependency of ", imp.vars[i]),
                       plot = TRUE, ice=FALSE)
}
print(pd[[1]], split = c(1, 1, 2, 2), more = TRUE)
print(pd[[2]], split = c(2, 1, 2, 2), more = TRUE)
print(pd[[3]], split = c(1, 2, 2, 2), more = TRUE)  





#### top UP VWC
#### top UP TEMP.SD
b = 501 # number of bootstraps

# Evaluate importance significance
( imp.reg <- ranger(x=df.top[,vars], y=df.top[,"vwc"],
                    importance = "permutation", 
                    num.trees=b) )
( imp <- na.omit(rf.ImpScale(imp.reg, scaling="p", n=1000)) )
( imp.vars <- imp[imp$importance > 0 & imp$pvalue <= 0.1,]$parameter )
# ( imp.vars <- imp[imp$importance > 0,]$parameter )

# Fit model
( bot.fit <- ranger(x=df.top[,imp.vars], y=df.top[,"vwc"],
                    importance = "permutation", num.trees=b,
                    keep.inbag=TRUE) )

cat("R-square of model:", bot.fit$r.squared, "\n")
cat("Prediction RMSE:", rmse(df.top[,"vwc"], bot.fit$predictions), "\n")
cat("Prediction correlation:", cor(df.top[,"vwc"], bot.fit$predictions), "\n")

# LINEAR MODEL
lm.top.vwc <- lm(vwc ~ cover_fraction, data = df.top)
summary(lm.top.vwc)
#rmse
sqrt(mean(lm.top.vwc$residuals^2))
# r corr
cor(df.top$vwc, df.top$cover_fraction)

# cohen's f2 from linear model
cohen_f2 <- function(fit, fit2){
    R2 <- summary(fit)$r.squared
    if(missing(fit2)) {
        R2/(1 - R2)
    } else {
        R2B <- summary(fit2)$r.squared
        (R2B - R2)/(1 - R2B)
    }
}

cohen_f2(lm.top.vwc)

# plot of soil moisture
x11(width = 4, height = 4)
ggplot(df.top, aes(x = cover_fraction, y = vwc))+
    geom_point()+
    theme_classic()+
    xlab("Change in Cover Fraction")+
    ylab("VWC [%]")+
    stat_smooth(method = "lm", se = FALSE)
# # Plot importance
# p <- imp[which(imp$parameter %in% imp.vars),]   
# p <- p[order(p$importance),] 
# x11()
# dotchart(p$importance, labels = p$parameter, 
#          main="Scaled Variable Importance", pch=19)  
# 
# #*****************************************************
# # Jackknife validation
# x <- df.top[,imp.vars]
# y <- df.top[,"vwc"] 
# midx <- 1:length(y)     
# error <- vector()
# serror <- vector()
# r.sqr <- vector() 
# j=0
# for(i in midx) { 
#     j=j+1
#     cat("Running obs ", j, " of ", length(midx), "\n")
#     rfv <- ranger(x=x[-i,], y=y[-i],
#                   num.trees = b, 
#                   importance="permutation")
#     error <- append(error, rfv$prediction.error)
#     serror <- append(serror, rmse(y[-i], rfv$predictions))
#     r.sqr <- append(r.sqr, rfv$r.squared)	
# }
# jackknife <- data.frame(y=y, error=error, 
#                         rmse=serror, r.sqr=r.sqr)
# 
# cat("Median Jackknife R-square:", median(r.sqr), "\n")
# cat("Median Jackknife RMSE:", median(serror), "\n")
# cat("Median Jackknife Prediction Error:", median(error), "\n")
# 
# # Plot Jackknife prediction and RMSE error
# x11()
# plot(density(error), main="Jackknife prediction error")
# abline(v=bot.fit$prediction.error, col="red")
# x11()
# plot(density(serror), main="Jackknife RMSE")
# abline(v=rmse(df.top[,"vwc"], bot.fit$predictions), col="red")
# 
# # Calculate and plot confidence region
# se <- suppressWarnings(predict(bot.fit, df.top[,imp.vars],
#                                type = "se", se.method = "infjack"))
# scalar = (se$predictions * 1.96) / max((se$predictions * 1.96)) 
# ci95 <- data.frame(
#     y = df.top[,"vwc"],
#     pred = bot.fit$predictions, 
#     std.err=se$se,
#     lower.ci = bot.fit$predictions - scalar, 
#     upper.ci = bot.fit$predictions + scalar)
# 
# # plot confidence 95% region 
# sort.idx <- sort(ci95$pred, index.return = TRUE)$ix
# x11()
# plot(ci95$pred[sort.idx], type="l", xaxt="n",
#      ylab="estimates", xlab="observations (n=14)", 
#      main="Estimate with 95% confidence region",  
#      ylim=c(min(ci95[,4:5]), max(ci95[,4:5])))
# lines(ci95$upper.ci[sort.idx], col="grey", lty=2)
# lines(ci95$lower.ci[sort.idx], col="grey", lty=2)
# 
# # Partial dependency and individual conditional expectation (ICE)
# pd=list()
# for(i in 1:length(imp.vars)){
#     pd[[i]] <- partial(bot.fit, pred.var=imp.vars[i], 
#                        train=df.top[,imp.vars],
#                        main=paste0("Partial dependency of ", imp.vars[i]),
#                        plot = TRUE, ice=FALSE)
# }
# x11()
# print(pd[[1]], split = c(1, 1, 2, 2), more = TRUE)
# print(pd[[2]], split = c(2, 1, 2, 2), more = TRUE)
# print(pd[[3]], split = c(1, 2, 2, 2), more = TRUE)  










###############################
# top UP VWC.sd
#### top UP TEMP.SD
b = 501 # number of bootstraps

# Evaluate importance significance
( imp.reg <- ranger(x=df.top[,vars], y=df.top[,"vwc.sd"],
                    importance = "permutation", 
                    num.trees=b) )
( imp <- na.omit(rf.ImpScale(imp.reg, scaling="p", n=100)) )
( imp.vars <- imp[imp$importance > 0 & imp$pvalue <= 0.1,]$parameter )
# ( imp.vars <- imp[imp$importance > 0,]$parameter )

# Fit model
( bot.fit <- ranger(x=df.top[,imp.vars], y=df.top[,"vwc.sd"],
                    importance = "permutation", num.trees=b,
                    keep.inbag=TRUE) )
imp.vars
cat("R-square of model:", bot.fit$r.squared, "\n")
cat("Prediction RMSE:", rmse(df.top[,"vwc.sd"], bot.fit$predictions), "\n")
cat("Prediction correlation:", cor(df.top[,"vwc.sd"], bot.fit$predictions), "\n")

# effect size
rf.top.vwc.sd <- randomForest(vwc.sd ~ richness + ba + gini + dbh_shan, data = df.top,
                            mtry = 2, ntree = 2000)
rf.regression.fit(rf.top.vwc.sd)

# Plot importance
p <- imp[which(imp$parameter %in% imp.vars),]   
p <- p[order(p$importance),] 
x11()
dotchart(p$importance, labels = p$parameter, 
         main="Scaled Variable Importance", pch=19)  

#*****************************************************
# Jackknife validation
x <- df.top[,imp.vars]
y <- df.top[,"vwc.sd"] 
midx <- 1:length(y)     
error <- vector()
serror <- vector()
r.sqr <- vector() 
j=0
for(i in midx) { 
    j=j+1
    cat("Running obs ", j, " of ", length(midx), "\n")
    rfv <- ranger(x=x[-i,], y=y[-i],
                  num.trees = b, 
                  importance="permutation")
    error <- append(error, rfv$prediction.error)
    serror <- append(serror, rmse(y[-i], rfv$predictions))
    r.sqr <- append(r.sqr, rfv$r.squared)	
}
jackknife <- data.frame(y=y, error=error, 
                        rmse=serror, r.sqr=r.sqr)

cat("Median Jackknife R-square:", median(r.sqr), "\n")
cat("Median Jackknife RMSE:", median(serror), "\n")
cat("Median Jackknife Prediction Error:", median(error), "\n")

# Plot Jackknife prediction and RMSE error
x11()
plot(density(error), main="Jackknife prediction error")
abline(v=bot.fit$prediction.error, col="red")
x11()
plot(density(serror), main="Jackknife RMSE")
abline(v=rmse(df.top[,"vwc.sd"], bot.fit$predictions), col="red")

# Calculate and plot confidence region
se <- suppressWarnings(predict(bot.fit, df.top[,imp.vars],
                               type = "se", se.method = "infjack"))
scalar = (se$predictions * 1.96) / max((se$predictions * 1.96)) 
ci95 <- data.frame(
    y = df.top[,"vwc.sd"],
    pred = bot.fit$predictions, 
    std.err=se$se,
    lower.ci = bot.fit$predictions - scalar, 
    upper.ci = bot.fit$predictions + scalar)

# plot confidence 95% region 
sort.idx <- sort(ci95$pred, index.return = TRUE)$ix
x11()
plot(ci95$pred[sort.idx], type="l", xaxt="n",
     ylab="estimates", xlab="observations (n=14)", 
     main="Estimate with 95% confidence region",  
     ylim=c(min(ci95[,4:5]), max(ci95[,4:5])))
lines(ci95$upper.ci[sort.idx], col="grey", lty=2)
lines(ci95$lower.ci[sort.idx], col="grey", lty=2)

# Partial dependency and individual conditional expectation (ICE)
pd=list()
for(i in 1:length(imp.vars)){
    pd[[i]] <- partial(bot.fit, pred.var=imp.vars[i], 
                       train=df.top[,imp.vars],
                       main=paste0("Partial dependency of ", imp.vars[i]),
                       plot = TRUE, ice=FALSE)
}
x11()
print(pd[[1]], split = c(1, 1, 2, 2), more = TRUE)
print(pd[[2]], split = c(2, 1, 2, 2), more = TRUE)
print(pd[[3]], split = c(1, 2, 2, 2), more = TRUE)  
print(pd[[4]], split = c(2, 2, 2, 2), more = TRUE)  





##########
# top-down par
b = 501 # number of bootstraps

# Evaluate importance significance
( imp.reg <- ranger(x=df.top[,vars], y=df.top[,"par"],
                    importance = "permutation", 
                    num.trees=b) )
( imp <- na.omit(rf.ImpScale(imp.reg, scaling="p", n=1000)) )
( imp.vars <- imp[imp$importance > 0 & imp$pvalue <= 0.1,]$parameter )
# ( imp.vars <- imp[imp$importance > 0,]$parameter )

# Fit model
( bot.fit <- ranger(x=df.top[,imp.vars], y=df.top[,"par"],
                    importance = "permutation", num.trees=b,
                    keep.inbag=TRUE) )

cat("R-square of model:", bot.fit$r.squared, "\n")
cat("Prediction RMSE:", rmse(df.top[,"par"], bot.fit$predictions), "\n")
cat("Prediction correlation:", cor(df.top[,"par"], bot.fit$predictions), "\n")

# LINEAR MODEL
lm.top.par <- lm(par ~ cover_fraction, data = df.top)
summary(lm.top.par)
#rmse
sqrt(mean(lm.top.par$residuals^2))
# r corr
cor(df.top$par, df.top$cover_fraction)

# cohen's f2 from linear model
cohen_f2 <- function(fit, fit2){
    R2 <- summary(fit)$r.squared
    if(missing(fit2)) {
        R2/(1 - R2)
    } else {
        R2B <- summary(fit2)$r.squared
        (R2B - R2)/(1 - R2B)
    }
}

cohen_f2(lm.top.par)

# plot of soil moisture
x11(width = 4, height = 4)
ggplot(df.top, aes(x = cover_fraction, y = par))+
    geom_point()+
    theme_classic()+
    xlab("Change in Cover Fraction")+
    ylab("faPAR")+
    stat_smooth(method = "lm", se = FALSE)

# Plot importance
p <- imp[which(imp$parameter %in% imp.vars),]   
p <- p[order(p$importance),] 
x11()
dotchart(p$importance, labels = p$parameter, 
         main="Scaled Variable Importance", pch=19)  

#*****************************************************
# Jackknife validation
x <- df.top[,imp.vars]
y <- df.top[,"par"] 
midx <- 1:length(y)     
error <- vector()
serror <- vector()
r.sqr <- vector() 
j=0
for(i in midx) { 
    j=j+1
    cat("Running obs ", j, " of ", length(midx), "\n")
    rfv <- ranger(x=x[-i,], y=y[-i],
                  num.trees = b, 
                  importance="permutation")
    error <- append(error, rfv$prediction.error)
    serror <- append(serror, rmse(y[-i], rfv$predictions))
    r.sqr <- append(r.sqr, rfv$r.squared)	
}
jackknife <- data.frame(y=y, error=error, 
                        rmse=serror, r.sqr=r.sqr)

cat("Median Jackknife R-square:", median(r.sqr), "\n")
cat("Median Jackknife RMSE:", median(serror), "\n")
cat("Median Jackknife Prediction Error:", median(error), "\n")

# Plot Jackknife prediction and RMSE error
x11()
plot(density(error), main="Jackknife prediction error")
abline(v=bot.fit$prediction.error, col="red")
x11()
plot(density(serror), main="Jackknife RMSE")
abline(v=rmse(df.top[,"par"], bot.fit$predictions), col="red")

# Calculate and plot confidence region
se <- suppressWarnings(predict(bot.fit, df.top[,imp.vars],
                               type = "se", se.method = "infjack"))
scalar = (se$predictions * 1.96) / max((se$predictions * 1.96)) 
ci95 <- data.frame(
    y = df.top[,"par"],
    pred = bot.fit$predictions, 
    std.err=se$se,
    lower.ci = bot.fit$predictions - scalar, 
    upper.ci = bot.fit$predictions + scalar)

# plot confidence 95% region 
sort.idx <- sort(ci95$pred, index.return = TRUE)$ix
x11()
plot(ci95$pred[sort.idx], type="l", xaxt="n",
     ylab="estimates", xlab="observations (n=14)", 
     main="Estimate with 95% confidence region",  
     ylim=c(min(ci95[,4:5]), max(ci95[,4:5])))
lines(ci95$upper.ci[sort.idx], col="grey", lty=2)
lines(ci95$lower.ci[sort.idx], col="grey", lty=2)

# Partial dependency and individual conditional expectation (ICE)
pd=list()
for(i in 1:length(imp.vars)){
    pd[[i]] <- partial(bot.fit, pred.var=imp.vars[i], 
                       train=df.top[,imp.vars],
                       main=paste0("Partial dependency of ", imp.vars[i]),
                       plot = TRUE, ice=FALSE)
}
x11()
print(pd[[1]], split = c(1, 1, 2, 2), more = TRUE)
print(pd[[2]], split = c(2, 1, 2, 2), more = TRUE)

