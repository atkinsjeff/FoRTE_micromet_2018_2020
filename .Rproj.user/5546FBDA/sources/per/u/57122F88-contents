# extra codes

# # organize 
# met %>%
#     # select(subplot_id, date, soil_temp, vwc) %>%
#     group_by(subplot_id, year) %>%
#     dplyr::filter(between(date, as.Date("2018-06-01"), as.Date("2018-08-31"))) %>%
#     dplyr::summarize(temp = mean(soil_temp, na.rm = TRUE),
#                      vwc = mean(vwc, na.rm = TRUE)) %>%
#     data.frame() -> x2018
# 
# met %>%
#     # select(subplot_id, date, soil_temp, vwc) %>%
#     group_by(subplot_id, year) %>%
#     dplyr::filter(between(date, as.Date("2018-06-01"), as.Date("2018-08-31"))) %>%
#     dplyr::summarize(temp.sd = sd(soil_temp, na.rm = TRUE),
#                      vwc.sd = sd(vwc, na.rm = TRUE)) %>%
#     data.frame() -> y2018
# 
# met %>%
#     # select(subplot_id, date, soil_temp, vwc) %>%
#     group_by(subplot_id, year) %>%
#     dplyr::filter(between(date, as.Date("2020-06-01"), as.Date("2020-08-31"))) %>%
#     dplyr::summarize(temp = mean(soil_temp, na.rm = TRUE),
#                      vwc = mean(vwc, na.rm = TRUE)) %>%
#     data.frame() -> x2020
# 
# met %>%
#     # select(subplot_id, date, soil_temp, vwc) %>%
#     group_by(subplot_id, year) %>%
#     dplyr::filter(between(date, as.Date("2020-06-01"), as.Date("2020-08-31"))) %>%
#     dplyr::summarize(temp.sd = sd(soil_temp, na.rm = TRUE),
#                      vwc.sd = sd(vwc, na.rm = TRUE)) %>%
#     data.frame() -> y2020
# 
# ####
# par %>%
#     #filter(between(timestamp, as.Date("2020-06-01"), as.Date("2020-08-31"))) %>%
#     filter(year != 2019) %>%
#     group_by(subplot_id, year) %>%
#     dplyr::summarize(par = mean(fapar, na.rm = TRUE)) %>%
#     data.frame() -> z2
# 
# 
# 
# 
# met.sums2018 <- merge(x2018, y2018)
# met.sums2020 <- merge(x2020, y2020)
# 
# met.sums <- rbind(met.sums2018, met.sums2020)
# met.sums <- merge(met.sums, z2)
# met.sums <- merge(met.sums, dis.meta.data)
# 
# 

require(ggplot2)


# forte colors
######
forte.pal <- forte_colors()
treat.lab <- c("Bottom-Up", "Top-Down")
names(treat.lab) <- c("B", "T")


x11(height = 8, width = 8)
ggplot(met.sums, aes(x = as.factor(disturbance_severity), y = temp, fill = as.factor(disturbance_severity)))+
    geom_boxplot(alpha = 0.5)+
    geom_point( size = 3, shape = 21, alpha = 1, 
                position = position_jitterdodge())+
    scale_fill_manual(values = forte.pal)+
    #scale_fill_viridis(discrete = TRUE)+
    theme_bw()+
    theme(axis.text.x = element_text(color = "grey20", size = 14, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          axis.text.y = element_text(color = "grey20", size = 14, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
          axis.title.x = element_text(color = "grey20", size = 20, angle = 0, hjust = .5, vjust = 0, face = "plain"),
          axis.title.y = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"))+
    ylab(expression("Growing Season"~ ~T[Soil]~( degree~C)))+
    xlab("Disturbance Severity (% Targeted Leaf Area Reduction)")+
    theme(legend.position = "none")+
    facet_grid(treatment ~ .,
               labeller = labeller(treatment = treat.lab))+
    theme(strip.text.y = element_text(size = 16))


x11(height = 8, width = 8)
ggplot(met.sums, aes(x = as.factor(disturbance_severity), y = vwc, fill = as.factor(disturbance_severity)))+
    geom_boxplot(alpha = 0.5)+
    geom_point( size = 3, shape = 21, alpha = 1, 
                position = position_jitterdodge())+
    scale_fill_manual(values = forte.pal)+
    #scale_fill_viridis(discrete = TRUE)+
    theme_bw()+
    theme(axis.text.x = element_text(color = "grey20", size = 14, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          axis.text.y = element_text(color = "grey20", size = 14, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
          axis.title.x = element_text(color = "grey20", size = 20, angle = 0, hjust = .5, vjust = 0, face = "plain"),
          axis.title.y = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"))+
    ylab("Growing Season  VWC (%)")+
    xlab("Disturbance Severity (% Targeted Leaf Area Reduction)")+
    theme(legend.position = "none")+
    facet_grid(treatment ~ .,
               labeller = labeller(treatment = treat.lab))+
    theme(strip.text.y = element_text(size = 16))


x11(height = 8, width = 8)
ggplot(met.sums, aes(x = as.factor(disturbance_severity), y = par, fill = as.factor(disturbance_severity)))+
    geom_boxplot(alpha = 0.5)+
    geom_point( size = 3, shape = 21, alpha = 1, 
                position = position_jitterdodge())+
    scale_fill_manual(values = forte.pal)+
    #scale_fill_viridis(discrete = TRUE)+
    theme_bw()+    
    theme(axis.text.x = element_text(color = "grey20", size = 14, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          axis.text.y = element_text(color = "grey20", size = 14, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
          axis.title.x = element_text(color = "grey20", size = 20, angle = 0, hjust = .5, vjust = 0, face = "plain"),
          axis.title.y = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"))+
    ylab(expression(paste("faPAR (", mu, m^2, sec^-1,")", sep="")))+
    xlab("Disturbance Severity (% Targeted Leaf Area Reduction)")+
    theme(legend.position = "none")+
    facet_grid(treatment ~ .,
               labeller = labeller(treatment = treat.lab))+
    theme(strip.text.y = element_text(size = 16))

ylab(expression(paste("faPAR (", mu, m^2, sec^-1,")", sep="")))


pred <- lm(par ~ dcc * treatment, data = df)

summary(lm(par ~ dcc, data = df.bottom))
summary(lm(par ~ dcc, data = df.top))


df$pred <- predict(pred, type = "response")

x11()
ggplot(df, aes(x = dcc, y = par, fill = as.factor(disturbance_severity)))+
    # geom_boxplot(alpha = 0.5)+
    geom_point( size = 3, shape = 21, alpha = 1)+
    scale_fill_manual(values = forte.pal)+
    #scale_fill_viridis(discrete = TRUE)+
    theme_bw()+
    geom_smooth(method = "lm", se = FALSE, data = df, aes(x = dcc, y = pred), color = "black")+
    ylab(expression(paste("faPAR (", mu, m^2, sec^-1,")", sep="")))+
    xlab("Disturbance Severity [% Targeted Leaf Area Reduction]")+
    theme(legend.position = "none")+
    facet_grid(treatment ~ .,
               labeller = labeller(treatment = treat.lab))




#### changes in structure
x11(height = 8, width = 8)
ggplot(df, aes(x = as.factor(disturbance_severity), y = drug, fill = as.factor(disturbance_severity)))+
    geom_boxplot(alpha = 0.5)+
    geom_point( size = 3, shape = 21, alpha = 1, 
                position = position_jitterdodge())+
    scale_fill_manual(values = forte.pal)+
    #scale_fill_viridis(discrete = TRUE)+
    theme_bw()+    
    theme(axis.text.x = element_text(color = "grey20", size = 14, angle = 0, hjust = .5, vjust = .5, face = "plain"),
          axis.text.y = element_text(color = "grey20", size = 14, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
          axis.title.x = element_text(color = "grey20", size = 20, angle = 0, hjust = .5, vjust = 0, face = "plain"),
          axis.title.y = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"))+
    ylab(expression(paste(widehat(R_c))))+
    xlab("Disturbance Severity (% Targeted Leaf Area Reduction)")+
    theme(legend.position = "none")+
    facet_grid(treatment ~ .,
               labeller = labeller(treatment = treat.lab))+
    theme(strip.text.y = element_text(size = 16))
















# 
# 
# 
# #### AMERIFLUX CHANGE PLOTS
# 
# x11(height = 8, width = 8)
# ggplot(csc, aes(x = as.factor(disturbance_severity), y = rugosity, fill = as.factor(disturbance_severity)))+
#     geom_boxplot(alpha = 0.5)+
#     geom_point( size = 3, shape = 21, alpha = 1, 
#                 position = position_jitterdodge())+
#     scale_fill_manual(values = forte.pal)+
#     #scale_fill_viridis(discrete = TRUE)+
#     theme_bw()+
#     theme(axis.text.x = element_text(color = "grey20", size = 14, angle = 0, hjust = .5, vjust = .5, face = "plain"),
#           axis.text.y = element_text(color = "grey20", size = 14, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
#           axis.title.x = element_text(color = "grey20", size = 20, angle = 0, hjust = .5, vjust = 0, face = "plain"),
#           axis.title.y = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"))+
#     #ylab(expression("Growing Season"~ ~T[Soil]~( degree~C)))+
#     ylab("Canopy Rugosity")+
#     xlab("Disturbance Severity (% Targeted Leaf Area Reduction)")+
#     theme(legend.position = "none")+
#     facet_grid(treatment ~ year,
#                labeller = labeller(treatment = treat.lab))+
#     theme(strip.text.y = element_text(size = 16))
# 
# 
# 
# x11(height = 8, width = 8)
# ggplot(csc, aes(x = as.factor(year), y = vai_mean, fill = as.factor(disturbance_severity)))+
#     geom_boxplot(alpha = 0.5)+
#     geom_point( size = 3, shape = 21, alpha = 1, 
#                 position = position_jitterdodge())+
#     scale_fill_manual(values = forte.pal)+
#     #scale_fill_viridis(discrete = TRUE)+
#     theme_bw()+
#     theme(axis.text.x = element_text(color = "grey20", size = 14, angle = 0, hjust = .5, vjust = .5, face = "plain"),
#           axis.text.y = element_text(color = "grey20", size = 14, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
#           axis.title.x = element_text(color = "grey20", size = 20, angle = 0, hjust = .5, vjust = 0, face = "plain"),
#           axis.title.y = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"))+
#     #ylab(expression("Growing Season"~ ~T[Soil]~( degree~C)))+
#     ylab("Leaf Area")+
#     xlab("Disturbance Severity (% Targeted Leaf Area Reduction)")+
#     theme(legend.position = "none")+
#     facet_grid(treatment ~ disturbance_severity,
#                labeller = labeller(treatment = treat.lab))+
#     theme(strip.text.y = element_text(size = 16))
# 
# 
# ####
# csc <- tibble::rownames_to_column(csc, "VALUE")
# csc %>%
#     dplyr::select(subplot_id, rugosity, year, disturbance_severity, treatment) %>%
#     group_by(subplot_id, year) %>%
#     summarise(mean.rug = mean(rugosity, na.rm = TRUE)) %>%
#     pivot_wider(names_from = year, values_from = mean.rug) %>%
#     data.frame() -> csc.wide
# 
# csc.wide <- merge(csc.wide, dis.meta.data)
# 
# x11(height = 8, width = 5)
# ggplot(csc.wide, aes(x = X2018, y = X2020, fill = as.factor(disturbance_severity), color = as.factor(disturbance_severity)))+
#     geom_point( size = 3, shape = 21, alpha = 1)+
#     scale_fill_manual(values = forte.pal)+
#     scale_color_manual(values = forte.pal)+
#     #scale_fill_viridis(discrete = TRUE)+
#     theme_bw()+
#     theme(axis.text.x = element_text(color = "grey20", size = 14, angle = 0, hjust = .5, vjust = .5, face = "plain"),
#           axis.text.y = element_text(color = "grey20", size = 14, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
#           axis.title.x = element_text(color = "grey20", size = 20, angle = 0, hjust = .5, vjust = 0, face = "plain"),
#           axis.title.y = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"))+
#     ylab("Canopy Complexity (2020)")+
#     xlab("Canopy Complexity (2018)")+
#     stat_smooth(method = "lm", se = FALSE)+
#     geom_abline(slope = 1, color = "#404040", size = 0.5,  alpha)+
#     xlim(c(0, 50))+
#     ylim(c(0, 50))+
#     theme(legend.position = "none")+
#     facet_grid(treatment ~ .,
#                labeller = labeller(treatment = treat.lab))+
#     theme(strip.text.y = element_text(size = 16))
# 
# 
# 
# 
# ####
# csc <- tibble::rownames_to_column(csc, "VALUE")
# csc %>%
#     dplyr::select(subplot_id, vai_mean, year, disturbance_severity, treatment) %>%
#     group_by(subplot_id, year) %>%
#     summarise(mean.rug = mean(vai_mean, na.rm = TRUE)) %>%
#     pivot_wider(names_from = year, values_from = mean.rug) %>%
#     data.frame() -> csc.vai
# 
# csc.vai <- merge(csc.vai, dis.meta.data)
# 
# x11(height = 8, width = 5)
# ggplot(csc.vai, aes(x = X2018, y = X2020, fill = as.factor(disturbance_severity), color = as.factor(disturbance_severity)))+
#     geom_point( size = 5, shape = 21, alpha = 0.5, color = "black")+
#     scale_fill_manual(values = forte.pal)+
#     scale_color_manual(values = forte.pal)+
#     #scale_fill_viridis(discrete = TRUE)+
#     theme_bw()+
#     theme(axis.text.x = element_text(color = "grey20", size = 14, angle = 0, hjust = .5, vjust = .5, face = "plain"),
#           axis.text.y = element_text(color = "grey20", size = 14, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
#           axis.title.x = element_text(color = "grey20", size = 20, angle = 0, hjust = .5, vjust = 0, face = "plain"),
#           axis.title.y = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"))+
#     ylab("Leaf Area (2020)")+
#     xlab("Leaf Area (2018)")+
#     stat_smooth(method = "lm", se = FALSE, size = 2.5)+
#     geom_abline(slope = 1, color = "#404040", size = 0.5)+
#     xlim(c(3, 8))+
#     ylim(c(3, 8))+
#     theme(legend.position = "none")+
#     facet_grid(treatment ~ .,
#                labeller = labeller(treatment = treat.lab))+
#     theme(strip.text.y = element_text(size = 16))

