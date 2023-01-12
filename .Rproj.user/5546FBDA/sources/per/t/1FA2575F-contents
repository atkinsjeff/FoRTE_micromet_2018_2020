#### plots for publication

library(ggplot2)

# Create data
data <- data.frame(
    x=LETTERS[1:26],
    y=abs(rnorm(26))
)

# Change baseline
ggplot(data, aes(x=x, y=y)) +
    geom_segment( aes(x=x, xend=x, y=1, yend=y), color="grey") +
    geom_point( color="orange", size=4) +
    theme_light() +
    theme(
        panel.grid.major.x = element_blank(),
        panel.border = element_blank(),
        axis.ticks.x = element_blank()
    ) +
    xlab("") +
    ylab("Value of Y")

# import the data:   not these data have been log transofmred and must be adjusted......see ~lines 30-40
comm <- read.csv("community_response_ratio.csv")
strut <- read.csv("structure_response_ratio.csv")
forest <- read.csv("forest_structure_rr.csv")

a <- comm[, c("subplot_id","disturbance_severity", "treatment", "drich")]
b <- strut[, c("subplot_id", "drug", "denl", "dp10", "dmoch", "dcc")]
c <- forest[, c(2, 11:16)]

df <- merge(a, b)
df <- merge(df, c)

# make regular, i.e. unlogtransform them
df$drich <- exp(df$drich)
df$drug <- exp(df$drug)
df$denl <- exp(df$denl)
df$dp10 <- exp(df$dp10)
df$dmoch <- exp(df$dmoch)
df$dcc <- exp(df$dcc)
df$dsd.dbh <- exp(df$dsd.dbh)
df$dba <- exp(df$dba)
df$dgini <- exp(df$dgini)
df$ddbh_shan <- exp(df$ddbh_shan)

# make radar graphs
df %>%
    group_by(disturbance_severity, treatment) %>%
    summarise_at(vars(drich:ddbh_shan), mean, na.rm = TRUE) %>%
    pivot_longer(cols = c(drich:ddbh_shan), names_to = "Variable", values_to = "Scalar") %>%
    data.frame() -> df.sums

# make color
df.sums$trend = ifelse(df.sums$Scalar >= 1, "Increase", "Decrease")
# Change baseline
df.sums %>%
    filter(treatment == "B") %>%
    data.frame -> bottom


# make a factor
bottom$disturbance_severity <- as.factor(bottom$disturbance_severity)
x11(width = 3, height = 6)
ggplot(bottom, aes(x = disturbance_severity, y = Scalar, color = trend)) +
    #geom_vline(xintercept = levels(as.factor(df.sums$Variable)), color = "grey")+
    geom_segment( aes(x = disturbance_severity, xend = disturbance_severity, y = 1, yend = Scalar), color = "black") +
    geom_point( size=4) +
    theme_minimal() +
    theme(
        panel.grid.major.x = element_blank(),
        panel.border = element_blank(),
        axis.ticks.x = element_blank())+
    theme(legend.position = "none")+
    ylim(c(0, 1.5))+
    xlab("") +
    ylab("Value of Y")+
    geom_hline(yintercept = 1, color = "black")+
    #theme(axis.text.x = element_text(angle = 90))+
    facet_wrap(.~ Variable, ncol = 2)




# Change baseline
df.sums %>%
    filter(treatment == "T") %>%
    data.frame -> top


# make a factor
top$disturbance_severity <- as.factor(top$disturbance_severity)
x11(width = 3, height = 6)
ggplot(top, aes(x = disturbance_severity, y = Scalar, color = trend)) +
    #geom_vline(xintercept = levels(as.factor(df.sums$Variable)), color = "grey")+
    geom_segment( aes(x = disturbance_severity, xend = disturbance_severity, y = 1, yend = Scalar), color = "black") +
    geom_point( size=4) +
    theme_minimal() +
    theme(
        panel.grid.major.x = element_blank(),
        panel.border = element_blank(),
        axis.ticks.x = element_blank())+
    theme(legend.position = "none")+
    ylim(c(0, 1.5))+
    xlab("") +
    ylab("Value of Y")+
    geom_hline(yintercept = 1, color = "black")+
    #theme(axis.text.x = element_text(angle = 90))+
    facet_wrap(.~ Variable, ncol = 2)



require(corrplot)
require(Hmisc)
df %>%
    filter(treatment == "B") %>%
    na.omit() %>%
    data.frame() -> bot.wide

m <- cor(bot.wide[4:13])
m2 <- rcorr(as.matrix(bot.wide[4:13]))
write.csv(m2, "bottom_up_full_correlation.csv")
x11()
corrplot(m, 
         order = 'AOE', addCoef.col = 'black', tl.pos = 'l',
         cl.pos = 'n',# Correlation matrix
         method = "color", # Correlation plot method
         type = "lower",    # Correlation plot style (also "upper" and "lower")
         diag = FALSE,      # If TRUE (default), adds the diagonal
         tl.col = "black", # Labels color
         bg = "white",     # Background color
         title = "",       # Main title
         col = NULL, row = NULL)    

df %>%
    filter(treatment == "T") %>%
    na.omit() %>%
    data.frame() -> top.wide

m2 <- cor(top.wide[4:13])


x11()
corrplot(m2, 
         order = 'AOE', addCoef.col = 'black', tl.pos = 'l',
         cl.pos = 'n',# Correlation matrix
         method = "color", # Correlation plot method
         type = "upper",    # Correlation plot style (also "upper" and "lower")
         diag = FALSE,      # If TRUE (default), adds the diagonal
         tl.col = "black", # Labels color
         bg = "white",     # Background color
         title = "",       # Main title
         col = NULL, row = NULL)    