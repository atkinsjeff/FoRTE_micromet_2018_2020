invisible(lapply(c("randomForest", "rfUtilities",    
                   "ranger", "pdp", "spatialEco"),require, 
                 character.only=TRUE))

# set working directory

source("./R/rf.ImpScale.R")

pfun <- function(object, newdata, i = 2) {
    predict(object, data = newdata)$predictions[,i]
}
rmse <- function(y,x) { sqrt(mean((y - x)^2)) }

set.seed(666)

#***********************************************
# read and filter top/bottom treatments
dat <- read.csv("atkins_test_data.csv")
df.bottom <- na.omit(dat[dat$treatment == "B",])
df.top <- na.omit(dat[dat$treatment == "T",])

# full variable list
all.vars <- c("richness","sw.index","evenness","dbh.sd","ba","gini","dbh_shan",
              "rugosity","top_rugosity", "enl","fhd","porosity","clumping_index",
              "mean_height_mean","moch","max_ht", "p10", "p25",
              "p50", "p75","p90","cover_fraction")

#***********************************************
# Evaluate colinearity and multicolinearity

# remove collinearity
cl.vars <- spatialEco::collinear(df.bottom[,all.vars], p = 0.85, 
                                 nonlinear = FALSE, p.value = 0.001)
if(length(cl.vars) > 0) {
    cat("Dropping colinear variables; ", cl.vars, "\n")  
    df.bottom <- df.bottom[,-which(names(df.bottom) %in% cl.vars)]
    vars <- all.vars[-which(all.vars %in% cl.vars)]
} else {
    cat("No collinear variable found", "\n")
}

# multicolinearity
( mc <- multi.collinear(df.bottom[,vars], p=0.005,  perm = TRUE, n = 999) )
mc.vars <- mc[which( mc$frequency > 5 ),]$variables 
if(length(mc.vars) > 0) { 
    df.bottom <- df.bottom[,-which(names(df.bottom) %in% mc.vars)]
} else {
    cat("No multicollinear variable found", "\n")
}

#***********************************************
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

save.image("C:/evans/USFS/SRS/bottom_model.RData")
