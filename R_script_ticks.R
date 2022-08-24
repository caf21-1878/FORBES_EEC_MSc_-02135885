# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#  Project R script - Craig Forbes MSc EEC                                    #
#  25/08/2022                                                                 #
#  Research Question: How will asymmetric warming affect host-seeking         #
#                      behaviour in Ixodes ricinus larvae                     #
#                                                                             #
#  Cited packages: rTPC (Padfield & O'Sullivan 2021)                          #
#                  MuMIn (Nakagawa & Schielzeth 2013)                         #
#                  lme4 (Bates et al. 2015)                                   #
#                  pROC (Robin et al. 2011)                                   #
#                  multcomp (Hothorn et al. 2008)                             #
#                  nlme (Pinheiro et al. 2017)                                #
#                                                                             #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


#### Set Up####

#clear environment
rm(list = ls())

#Set working directory
setwd("/Users/Craig Forbes/Desktop/PG/Project")

#Install packages#
install.packages("lme4")
install.packages("usdm")
install.packages("lmerTest")
install.packages("psych")
install.packages("dplyr")
install.packages("tidyr")
install.packages("ggplot2")
install.packages("nls.multstart")
remotes::install_github("padpadpadpad/rTPC")
install.packages("mgcv")
install.packages("ggiraph")
install.packages("ggiraphExtra")
install.packages("plyr")
install.packages("purrr")
install.packages("multcomp")
install.packages("broom")
install.packages("MuMIn")
install.packages("blmeco")
install.packages("rnaturnalearth")
install.packages("sf")
install.packages("pROC")
install.packages("ggeffects")
install.packages("sjmisc")
install.packages("ggpubr")
install.packages("sjPlot")
install.packages("sjlabelled")
install.packages("nlme")


#Load packages 
library(lme4)
library(usdm)
library(lmerTest)
library(psych)
library(dplyr)
library(tidyr)
library(ggplot2)
library(nls.multstart)
library(rTPC)
library(mgcv)
library(ggiraph)
library(ggiraphExtra)
library(purrr)
library(plyr)
library(purrr)
library(multcomp)
library(broom)
library(MuMIn)
library(blmeco)
library(rnaturalearth)
library(sf)
library(pROC)
library(ggeffects)
library(sjmisc)
library(ggpubr)
library(sjPlot)
library(sjlabelled)
library(nlme)
library(grid)
library(gridExtra)
library(lattice)




#### Study site map ####

#Load in a map of the UK
UK <- map_data(map = "world", region = "UK")

#Creat a custom plot of the UK, highlighting Berkshire
  big_map <- ggplot(data = UK, aes(x = long, y = lat, group = group)) + 
    geom_polygon(colour = "black", fill =NA, size = 1.3) +
    coord_map() +
    geom_rect(aes(xmin = -1, xmax = -0.46, ymin = 51.3, ymax = 51.59), color = 
                "red", fill = "grey", size = 1.2) + 
    labs(x = "Latitude", y = "Longitude", size = 2) + 
    theme(plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm")) +
    theme(axis.title.y = element_text(vjust=4, size = 27)) + 
    theme(axis.title.x = element_text(vjust=-4, size = 27)) + 
    theme(axis.text = element_text(size = 27, colour = "black")) + 
    theme(legend.position="none")+
    ylim(50, 60) +
    theme(panel.background = element_blank()) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 3, size = 35)) +
    theme(axis.line = element_line(colour = "black"),
          panel.background = element_blank()) + 
    theme(axis.line.x.bottom=element_line(size=1)) + 
    theme(axis.line.y.left =element_line(size=1))
  
#Check the map 
  big_map
  
#load in data on transect locations
  coords <- read.csv("coords.csv", header = T)
  
#Change the names to match the UK data
  colnames(coords)[1] <- "long"
  coords$group <- seq(1,1,98)
  
#Create a zoomed in map of Berkshire, with the transect points
  small_map <- ggplot(data = UK, aes(x = long, y = lat, group = group)) + 
    geom_polygon(colour = "black", fill =NA, size = 1.3) +
    coord_map() +
    geom_rect(aes(xmin = -1, xmax = -0.46, ymin = 51.3, ymax = 51.59), color = 
                "red", fill = "grey", size = 1.5) + 
    xlim(-1, -0.46) + 
    ylim(51.3, 51.59) + 
    geom_point(data = coords, aes(x = coords$lat, y = coords$long), size = 4, pch = 4) +
    labs(x = "Latitude", y = "Longitude", size = 2) + 
    theme(plot.margin = unit(c(1,1,1,1), "cm")) +
    theme(axis.title = element_blank()) + 
    theme(axis.text = element_blank()) + 
    theme(legend.position="none")+
    theme(panel.background = element_blank()) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 3, size = 35)) +
    theme(axis.line = element_blank(),
          axis.ticks = element_blank(),
          panel.background = element_blank())
 
#Check the map 
  small_map

#Combine the maps into one figure 
  ggarrange(big_map, small_map, nrow = 1, ncol = 2, widths = c(1.7,1), heights = c(2,1)) 
  
#### Data Wrangling ####

#Load in the height and questing numbers data
  p_data <- read.csv("tick_props_R.csv", header = T)
  h_data <-read.csv("tick_height_R.csv", header = T)

#Change the column names
  colnames(p_data) [1] <- "pop"
  colnames(h_data) [1] <- "pop"

#Check the structure of the dataframes
  str(h_data)
  str(p_data)
  
#Convert the 6th height column to numeric
  h_data$height.6 <- as.numeric(h_data$height.6)

#Edit the height data to have one column for questing heights
  h_data<-pivot_longer(h_data, names_to = "qh", values_to = "q_height", cols=height:height.11, values_drop_na = T)
  h_data <- h_data[-c(9)]

#Check dataframe structures
  str(h_data)
  str(p_data)

#Convert proportion data to binary response#

#Remove the percentage column from the dataframe
  p_data_npq <- subset(p_data, select = -c(perc_q))
  
#Rename the population column
  colnames(p_data_npq) [1] <- "pop"


#Create a function to extract binary responses (0,1) from binomial data (proportion 
# of questing ticks) whilst keeping other columns in the dataframe
  binary_data = pmap_dfr(p_data_npq, 
                      function(pop, min_t, rec_time, site, lab_time, time_since_min
                               , sample_t, mean_nt,numb_q, numb_tick) {
                        data.frame(pop = pop,
                                   min_t = min_t,
                                   rec_time = rec_time,
                                   site = site, 
                                   lab_time = lab_time, 
                                   time_since_min = time_since_min, 
                                   sample_t = sample_t, 
                                   mean_nt = mean_nt, 
                                   q = c( rep(1, numb_q),
                                             rep(0, numb_tick - numb_q) ) ) })
#Check that the function worked
  str(binary_data)

#### Questing Height Analysis ####

#Calculate overall summary statistics for height data 
  mean(h_data$q_height)
  sd(h_data$q_height)
  min(h_data$q_height)
  max(h_data$q_height)

#Calculate summary statistics for each treatment 
  h_summary <- dplyr::group_by(h_data, min_t) %>%
    dplyr::summarise(., sd = sd(q_height), q_height = mean(q_height)) %>%
    ungroup()

#Check that the stats were calculated
  str(h_summary)

#convert minimum temperature to a factor for analysis
  h_data$min_t <- as.factor(h_data$min_t)

#Run a linear mixed model with site, sample temperature and population as random 
# effects 
  lmm1 <- lmer(q_height ~ min_t + (1|site) + (1|sample_t) + (1|pop), data = h_data)

#Check the summary output
  summary(lmm1)

#Do a type II anova to check overall impact of minimum temperature
  car::Anova(lmm1, type = 2)
  
#Create a tabel showing the output for use in the report
  tab_model(lmm1, show.se = T, show.fstat = T, auto.label = F, dv.labels =
            "LMM (Questing Height)", pred.labels = c("5C (reference category)",
                                                            "7.5C", "10C", "12.5C", "15C", "17.5C")
          , string.pred = "Treatment", string.se = "SE")


#Calculate the R squared
  r.squaredGLMM(lmm1)

#Get means and confidence intervals for each group using ggpredict
  h <- as.data.frame(ggpredict(lmm1, "min_t"))
  
#Create a plot for questing height 
  lmm_plot <- ggplot() + geom_point(aes(x = h_summary$min_t, y = h_summary$q_height), size = 2.5) + 
    geom_errorbar(aes(x = h_summary$min_t, y = h_summary$q_height, 
                    ymin = conf.low, ymax = conf.high), h, width = 0.2, size = 1) +
    labs(x = "Minimum Temperature (ºC)", y = "Questing Height (cm)", size = 2) + 
    theme(plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm")) +
    theme(axis.title.y = element_text(vjust=4, size = 27)) + 
    theme(axis.title.x = element_text(vjust=-4, size = 27)) + 
    theme(axis.text = element_text(size = 27, colour = "black")) + 
    theme(legend.position="none")+
    theme(panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill=NA, size=1.5)) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 3, size = 35)) +
    theme(axis.line = element_line(colour = "black"),
        panel.background = element_blank())

#Perform a TUKEY HSD test to compare treatments to each other. 
  h_tukey<- glht(lmm1, mcp(min_t="Tukey"))

#Check the summary
  summary(h_tukey)

#Calculate confidence intervals for differences and create a dataframe
  h_tukey_df <- confint(h_tukey) %>% tidy

#Convert to a dataframe
  h_tukey_df <- as.data.frame(h_tukey_df)
  
#rename the treatment comparisons
  h_tukey_df$contrast<- c("5/7.5", "5/10", "5/12.5", "5/15", "5/17.5", 
                        "7.5/10", "7.5/12.5", "7.5/15", "7.5/17.5", 
                        "10/12.5", "10/15", "10/17.5", "12.5/15", "12.5/17.5", 
                        "15/17.5") 
  
#Add a column with asterixes to represent level of significance for each treatment 
  h_tukey_df$sig <- c("", "*", "***", "***", "*", "*", "***", "***", "*", "", "", "", "", "", "" )

#Create a plot showing the confidence intervals of pairwise comparisons 
  ht_plot <- ggplot()+ geom_point(aes(x = estimate, y = contrast), h_tukey_df) + 
    geom_errorbar(aes(x = estimate, y = contrast, 
                    xmin = conf.low, xmax = conf.high),h_tukey_df, width = 0.2, size = 1) +
    geom_text(aes(x = estimate, y = contrast, label = sig, vjust = -0.3, size = 10), h_tukey_df) +
    labs(x = "Effect Size", y = "Treatment Comparison")+ 
    geom_vline(xintercept = 0, linetype = "dashed", size = 1, alpha = 0.4, colour = "black")+
    xlim(-2.5, 2.5) +
    theme(plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm")) +
    theme(axis.title.y = element_text(vjust=4, size = 27)) + 
    theme(axis.title.x = element_text(vjust=-4, size = 27)) + 
    theme(axis.text.y = element_text(size = 20, colour = "black")) +
    theme(axis.text.x = element_text(size = 27, colour = "black")) +
    theme(legend.position="none")+
    theme(panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill=NA, size=1.5)) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 3, size = 35)) +
    theme(axis.line = element_line(colour = "black"),
        panel.background = element_blank()) +
    theme(panel.grid.major.y = element_line(colour = "grey", size = 0.4))

# Combine the height plot and the tukey plot for use in the report 
    ggarrange(lmm_plot, ht_plot, nrow = 1, ncol = 2, labels = "AUTO", 
          hjust = -1.5, vjust = 3, font.label = list(size = 25, 
              face = "bold.italic", color ="black"), widths = c(1, 1.05))



#### Questing Probabality Analysis ####

#### Calulate summary statistics for each group
  p_summary <- dplyr::group_by(p_data, min_t) %>%
   dplyr::summarise(., sd = sd(perc_q), perc_q = mean(perc_q)) %>%
   ungroup()

#Add mean temperature to the summary stats 
    p_summary$mean_nt <- c(10, 11.5, 13, 14.5, 16, 17.5)
 
#Check summary stats        
  str(p_summary)
  
#Convert minimum temperature to a factor for analysis
  binary_data$min_t <- as.factor(binary_data$min_t)

#Check structure of dataframe 
  str(binary_data)

#Split data into training (75%) and test (25%) data to assess AUC  
sample <- sample.int(n = nrow(binary_data), size = floor(.75*nrow(binary_data)), replace = F)
train <- binary_data[sample, ]
test  <- binary_data[-sample, ]

#Dataset is randomly split every time code is run so output parameters will vary slightly 
#to those presented in the study

#Run a GLMM on the training data to establish questing probabaility for each treatment
  glm1 <- glmer(formula = q ~ min_t + (1|pop) + (1|site) + (1|sample_t), 
              family = binomial, data = train)
  
#Check summary output
  summary(glm1)

#Run a type II anova to check overall impact 
  car::Anova(glm1)

#Calculate  pseudo R2
  r.squaredGLMM(glm1)

#Predict the model over the retained test data
  predicted <- predict(glm1, test, type="response")

#Assess the model's capacity to predict the test data by calculating AUC   
  auc(test$q, predicted)
  
#Generate an ROC curve for presentation in the report  
  roc1 <- roc(test$q, predicted)

#Plot the ROC curve for use in the report
 ggroc(roc1) + 
   geom_segment(aes(x = 1, y = 0, xend = 0, yend = 1), size = 1, lineend = "butt", linetype = "dashed") +
   geom_text(aes(x = 0.74, y = 0.7), label = "AUC = 0.6346", size = 10, vjust = 0.5) +
   labs(x = "False Positive Rate", y = "True Positive Rate") + 
   theme(plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm")) +
   theme(axis.title.y = element_text(vjust=4, size = 27)) + 
   theme(axis.title.x = element_text(vjust=-4, size = 27)) + 
   theme(axis.text = element_text(size = 27, colour = "black")) + 
   theme(legend.position="none")+
   theme(panel.background = element_blank(), 
         panel.border = element_rect(colour = "black", fill=NA, size=1.5)) +
   theme(plot.title = element_text(hjust = 0.5, vjust = 3, size = 35)) +
   theme(axis.line = element_line(colour = "black"),
         panel.background = element_blank())
 
#Calculate the probabaility of questing at each treatment
  g<- as.data.frame(ggpredict(glm1, "min_t", interactive = T))

#Plot model predictions 
  glmm_plot <- ggplot() + geom_point(aes(x = g$x, y = g$predicted), size = 2.5) + 
    geom_errorbar(aes(x = g$x, y = g$predicted, 
                    ymin = g$conf.low, ymax = g$conf.high),g, width = 0.2, size = 1) +
    labs(x = "Nocturnal Minimum Temperature (ºC)", y = "Probability of Questing (%)") + 
    ylim(0.2, 0.7) +
    theme(plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm")) +
    theme(axis.title.y = element_text(vjust=4, size = 27)) + 
    theme(axis.title.x = element_text(vjust=-4, size = 27)) + 
    theme(axis.text = element_text(size = 27, colour = "black")) + 
    theme(legend.position="none")+
    theme(panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill=NA, size=1.5)) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 3, size = 35)) +
    theme(axis.line = element_line(colour = "black"),
      panel.background = element_blank())

#Create a table of model output for use in report
  tab_model(glm1, show.se = T, show.fstat = T, auto.label = F, dv.labels =
            "GLMM (Questing Probabaility)", pred.labels = c("5C (reference category)",
                                                            "7.5C", "10C", "12.5C", "15C", "17.5C")
          , string.pred = "Treatment", string.se = "SE")



##Determine the post-hoc comparisons of interest
p_tukey<- glht(glm1, mcp(min_t="Tukey"))
summary(p_tukey)
 confint(p_tukey)
# Convert output to dataframe for plotting 
 p_tukey_df <- confint(p_tukey) %>% tidy
 p_tukey_df <- as.data.frame(p_tukey_df)
 p_tukey_df$contrast<- c("5/7.5", "5/10", "5/12.5", "5/15", "5/17.5", 
                          "7.5/10", "7.5/12.5", "7.5/15", "7.5/17.5", 
                          "10/12.5", "10/15", "10/17.5", "12.5/15", "12.5/17.5", 
                          "15/17.5")
 p_tukey_df$sig <- c("**", "***", "***", "", "", "***", "***", "", "***", "", "***", "***", "***", "***", "**")
 

 #generate CI plot for pairwise comparisons
 pt_plot <- ggplot()+ geom_point(aes(x = estimate, y = contrast), p_tukey_df) + 
   geom_errorbar(aes(x = estimate, y = contrast, 
                     xmin = conf.low, xmax = conf.high),p_tukey_df, width = 0.2, size = 1) +
   labs(x = "Effect Size", y = "Treatment Comparison")+ 
   geom_text(aes(x = estimate, y = contrast, label = sig, vjust = -0.3, size = 10), p_tukey_df) +
   xlim(-2, 2) + 
   geom_vline(xintercept = 0, linetype = "dashed", size = 1, alpha = 0.4, colour = "black")+
   theme(plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm")) +
   theme(axis.title.y = element_text(vjust=4, size = 27)) + 
   theme(axis.title.x = element_text(vjust=-4, size = 27)) + 
   theme(axis.text.y = element_text(size = 20, colour = "black")) +
   theme(axis.text.x = element_text(size = 27, colour = "black")) +
   theme(legend.position="none")+
   theme(panel.background = element_blank(), 
         panel.border = element_rect(colour = "black", fill=NA, size=1.5)) +
   theme(plot.title = element_text(hjust = 0.5, vjust = 3, size = 35)) +
   theme(axis.line = element_line(colour = "black"),
         panel.background = element_blank()) +
   theme(panel.grid.major.y = element_line(colour = "grey", size = 0.4))
 
#Combine Glm plot and tukey plot 
  ggarrange(glmm_plot, pt_plot, nrow = 1, ncol = 2, labels = "AUTO", 
           hjust = -1.5, vjust = 3, font.label = list(size = 25, 
              face = "bold.italic", color ="black"), widths = c(1, 1.05))
 
 #Asymmetric vs constant####
 
#Create a range of temperatures used by Gilbert et al. (2014)
  temp <- seq(6, 15, 1)
  
#Add the associated questing % values 
  qp <- c(0.125, 0.23, 0.35, 0.48, 0.59, 0.71, 0.77, 0.83, 0.84, 0.89)
  
#Add strandard errors and convert to SD for use in report 
  sd <- c(0.04*sqrt(4), 0.009*sqrt(4), 0.025*sqrt(4), 0.03*sqrt(4), 0.02*sqrt(4), 0.018*sqrt(4), 0.015*sqrt(4), 0.009*sqrt(4), 0.012*sqrt(4), 0.01*sqrt(4))

#Combine into a single dataframe
  Gilb_data <- data.frame(temp, qp, sd) 

#Generate a plot showing Gilbert et al. (2014) values alongside my own results 
comp_plot <- ggplot () + 
  geom_point(aes(mean_nt, perc_q), p_summary, pch = 19, size = 3) +
  geom_point(aes(temp, qp), Gilb_data, pch = 23, bg = "grey", size = 3) + 
  geom_segment(aes(x = 10, y = 0.15, xend = 17.5, yend = 0.15), size = 1, lineend = "butt") + 
  scale_x_continuous(breaks = round(seq(min(Gilb_data$temp), 18, by = 1),1)) +
  geom_errorbar(aes(x = mean_nt, y = perc_q, 
                    ymin = perc_q - sd, ymax = perc_q +sd),p_summary, width = 0.08, size = 0.8, alpha = 0.5) +
  geom_errorbar(aes(x = temp, y = qp, 
                    ymin = qp - sd, ymax = qp +sd),Gilb_data, width = 0.08, size = 0.8, alpha = 0.6) +
  annotate("text", x = 13.6, y = 0.1, label = "Mean Nocturnal Temperature (10-17.5ºC)", size = 8) + 
  labs(x = "Temperature (ºC)", y = "Questing Proportion") + 
  ylim(0, 1) +
  theme(plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm")) +
  theme(axis.title.y = element_text(vjust=4, size = 27)) + 
  theme(axis.title.x = element_text(vjust=-4, size = 27)) + 
  theme(axis.text = element_text(size = 27, colour = "black")) + 
  theme(legend.position="none")+
  theme(panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill=NA, size=1.5)) +
  theme(plot.title = element_text(hjust = 0.5, vjust = 3, size = 35)) +
  theme(axis.line = element_line(colour = "black"),
        panel.background = element_blank())

#Check the plot 
  comp_plot
  
 
#### Thermal Performance ####

#Disclaimer: This section of code was adapted from Dr. Padfields gudie on 
  #using the rTPC package

#Isolate temperature and questing data from the wider dataframe and rename columns
  d <-  data.frame(p_data$min_t,p_data$perc_q)
  colnames(d)[1] <- "temp"
  colnames(d)[2] <- "rate"


#Fit the 6 chosen TPC models to the data 
  d_fits <- nest(d, data = c(temp, rate)) %>%
    dplyr::mutate(boatman = map(data, ~nls_multstart(rate~boatman_2017(temp = temp, rmax, tmin, tmax, a,b),
                                            data = .x,
                                            iter = c(4,4,4,4,4),
                                            start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'boatman_2017') - 10,
                                            start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'boatman_2017') + 10,
                                            lower = get_lower_lims(.x$temp, .x$rate, model_name = 'boatman_2017'),
                                            upper = get_upper_lims(.x$temp, .x$rate, model_name = 'boatman_2017'),
                                            supp_errors = 'Y',
                                            convergence_count = FALSE)),
         gaussian = map(data, ~nls_multstart(rate~gaussian_1987(temp = temp, rmax, topt, a),
                                             data = .x,
                                             iter = c(4,4,4),
                                             start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'gaussian_1987') - 10,
                                             start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'gaussian_1987') + 10,
                                             lower = get_lower_lims(.x$temp, .x$rate, model_name = 'gaussian_1987'),
                                             upper = get_upper_lims(.x$temp, .x$rate, model_name = 'gaussian_1987'),
                                             supp_errors = 'Y',
                                             convergence_count = FALSE)),
         oneill = map(data, ~nls_multstart(rate~oneill_1972(temp = temp, rmax, ctmax, topt, q10),
                                           data = .x,
                                           iter = c(4,4,4,4),
                                           start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'oneill_1972') - 10,
                                           start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'oneill_1972') + 10,
                                           lower = get_lower_lims(.x$temp, .x$rate, model_name = 'oneill_1972'),
                                           upper = get_upper_lims(.x$temp, .x$rate, model_name = 'oneill_1972'),
                                           supp_errors = 'Y',
                                           convergence_count = FALSE)),
         quadratic = map(data, ~nls_multstart(rate~quadratic_2008(temp = temp, a, b, c),
                                              data = .x,
                                              iter = c(4,4,4),
                                              start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'quadratic_2008') - 0.5,
                                              start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'quadratic_2008') + 0.5,
                                              lower = get_lower_lims(.x$temp, .x$rate, model_name = 'quadratic_2008'),
                                              upper = get_upper_lims(.x$temp, .x$rate, model_name = 'quadratic_2008'),
                                              supp_errors = 'Y',
                                              convergence_count = FALSE)),
         ratkowsky = map(data, ~nls_multstart(rate~ratkowsky_1983(temp = temp, tmin, tmax, a, b),
                                    data = d,
                                    iter = c(4,4,4,4),
                                    start_lower = get_start_vals(d$temp, d$rate, model_name = 'ratkowsky_1983') - 10,
                                    start_upper = get_start_vals(d$temp, d$rate, model_name = 'ratkowsky_1983') + 10,
                                    lower = get_lower_lims(d$temp, d$rate, model_name = 'ratkowsky_1983'),
                                    upper = get_upper_lims(d$temp, d$rate, model_name = 'ratkowsky_1983'),
                                    supp_errors = 'Y',
                                    convergence_count = FALSE)),
         sharpeschoolhigh = map(data, ~nls_multstart(rate~sharpeschoolhigh_1981(temp = temp, r_tref,e,eh,th, tref = 15),
                                                     data = .x,
                                                     iter = c(4,4,4,4),
                                                     start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'sharpeschoolhigh_1981') - 10,
                                                     start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'sharpeschoolhigh_1981') + 10,
                                                     lower = get_lower_lims(.x$temp, .x$rate, model_name = 'sharpeschoolhigh_1981'),
                                                     upper = get_upper_lims(.x$temp, .x$rate, model_name = 'sharpeschoolhigh_1981'),
                                                     supp_errors = 'Y',
                                                     convergence_count = FALSE)))


#Stack the models 
  d_stack <- dplyr::select(d_fits, -data) %>%
  pivot_longer(., names_to = 'model_name', values_to = 'fit', boatman:sharpeschoolhigh)

#Generate a new data frame to predict the models over
  newdata <- tibble(temp = seq(-5, 30, length.out = 10000))

#Predict the each model over the new data frame
  d_preds <- d_stack %>%
  dplyr::mutate(., preds = map(fit, broom::augment,  newdata = newdata)) %>%
  dplyr::select(-fit) %>%
  unnest(preds)

#Calculate evaluative parameters for each model 
  d_ic <- d_stack %>%
  mutate(., info = map(fit, broom::glance),
         AICc = map_dbl(fit, MuMIn::AICc)) %>%
  dplyr::select(-fit) %>%
  unnest(info) %>%
  dplyr::select(model_name, sigma, AIC, AICc, BIC, df.residual)

#Select relevant paramters and coerce into a new dataframe
  table.aic <- d_ic %>% dplyr::select(model_name, AIC, df.residual)

#Rename columns 
  colnames(table.aic)[1] <- "Model"
  colnames(table.aic) [3] <- "df"

#Convert into a table for use in report 
  table.aic.p <- ggtexttable(table.aic, rows = NULL, theme = ttheme())

#Check the table 
  table.aic.p

#Calculate thermal parameters for each model
  params <- d_stack %>%
    mutate(., params = map(fit, calc_params)) %>%
    dplyr::select(-fit) %>%
    unnest(params)

#Check parameter estimates 
  params
#Filter the AIC to select the best model 
  best_model = filter(d_ic, AICc == min(AICc)) %>% pull(model_name)
  
#Print the best model
  best_model

#Create a plot showing the best model, compared to the others
multi_model_plot <- ggplot(d_preds, aes(temp, .fitted)) +
  geom_line(aes(group = model_name), col = 'grey50', alpha = 0.5, size = 1) +
  geom_line(data = filter(d_preds, model_name == best_model), col = "blue", size = 1.5) + 
  geom_point(aes(min_t, perc_q), p_summary) +
  geom_errorbar(aes(x = min_t, y = perc_q, 
                    ymin = perc_q - sd, ymax = perc_q + sd),p_summary, width = 0.2, size = 1) +
  ylim(0, 0.8) +
  xlim(-2.5, 25) +
  labs(x = 'Nocturnal Minimum Temperature (ºC)',
       y = 'Questing Proportion') + 
  theme(plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm")) +
  theme(axis.title.y = element_text(vjust=4, size = 27)) + 
  theme(axis.title.x = element_text(vjust=-4, size = 27)) + 
  theme(axis.text = element_text(size = 27, colour = "black")) + 
  theme(legend.position="none")+
  theme(panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill=NA, size=1.5)) +
  theme(plot.title = element_text(hjust = 0.5, vjust = 3, size = 35)) +
  theme(axis.line = element_line(colour = "black"),
        panel.background = element_blank()) 

#Check plot 
  multi_model_plot

# Bootstrapping 
# Fit the Gaussian model again and predict it over the same range (-5 - 30)
  d_fit <- nest(d, data = c(temp, rate)) %>%
   mutate(gaussian = map(data, ~nls_multstart(rate~gaussian_1987(temp = temp, rmax, topt, a),
                                             data = .x,
                                             iter = c(4,4,4),
                                             start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'gaussian_1987') - 10,
                                             start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'gaussian_1987') + 10,
                                             lower = get_lower_lims(.x$temp, .x$rate, model_name = 'gaussian_1987'),
                                             upper = get_upper_lims(.x$temp, .x$rate, model_name = 'gaussian_1987'),
                                             supp_errors = 'Y',
                                             convergence_count = FALSE)),
         #create new dataframe for predictions 
         new_data = map(data, ~tibble(temp = seq(-5, 30, length.out = 10000))),
         #predict the model over the new dataframe 
         preds =  map2(gaussian, new_data, ~augment(.x, newdata = .y)))



# Refit the Gaussian model using nlsLM
  fit_nlsLM <- minpack.lm::nlsLM(rate~gaussian_1987(temp = temp, rmax, topt, a),
                               data = d,
                               start = coef(d_fit$gaussian[[1]]),
                               lower = get_lower_lims(d$temp, d$rate, model_name = 'gaussian_1987'),
                               upper = get_upper_lims(d$temp, d$rate, model_name = 'gaussian_1987'),
                               weights = rep(1, times = nrow(d)))

#Bootstrap the model using case resampling
  boot1 <- car::Boot(fit_nlsLM, method = 'case')

#Check the data
  head(boot1$t)

#Create predictions for each bootstrapped model
  boot1_preds <- boot1$t %>%
    as.data.frame() %>%
    drop_na() %>%
    dplyr::mutate(iter = 1:dplyr::n()) %>%
    group_by_all() %>%
    do(data.frame(temp = seq(-5, 30, length.out = 1000))) %>%
    ungroup() %>%
    mutate(pred = gaussian_1987(temp, rmax, topt, a))


#Generate confidence intervals for the model predictions
  boot1_conf_preds <- group_by(boot1_preds, temp) %>%
  dplyr::mutate(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()

#Create a plot to show the Gaussian model w2ith 95% CI of predictions 
  boot_plot <- ggplot() +
    geom_line(aes(x = temp, y = .fitted), data = filter(d_preds, model_name == best_model), col = "blue", size = 0.8) +
    geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot1_conf_preds, fill = 'red', alpha = 0.3) +
    geom_point(aes(min_t, perc_q), p_summary) +
    geom_errorbar(aes(x = min_t, y = perc_q, 
                    ymin = perc_q - sd, ymax = perc_q +sd),p_summary, width = 0.2, size = 1) +
    ylim(0, 0.8) +
    xlim(-2.5, 25) +
    labs(x = 'Nocturnal Mean Temperature (ºC)') + 
    theme(plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm")) +
    theme(axis.title.y = element_blank()) + 
    theme(axis.title.x = element_text(vjust=-4, size = 27)) + 
    theme(axis.text  = element_text(size = 27, colour = "black")) +
    theme(legend.position="none")+
    theme(panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill=NA, size=1.5)) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 3, size = 35)) +
    theme(axis.line = element_line(colour = "black"),
        panel.background = element_blank())
  
#Check the plot 
  boot_plot
  
#Create a figure with both TPC plots 
  figure <- ggarrange(multi_model_plot + rremove("xlab") + rremove("ylab"), boot_plot + rremove ("xlab"), nrow = 1, ncol = 2, labels = "AUTO", 
          hjust = -0.5, vjust = 3, font.label = list(size = 25, 
          face = "bold.italic", color ="black"), widths = c(1.05,1))

#Annotate the figure to have common axes
  final_TPC_plot <- annotate_figure(figure, left = textGrob("Questing Proportion", rot = 90, gp = gpar(cex = 3)),
                bottom = textGrob("Minimum Nocturnal Temperature (ºC)", vjust = -0.6, gp = gpar(cex = 3)))

# Check the figure 
  final_TPC_plot
#### Population TPCs####

#Add population to the dataframe used for TPC analysis 
  d$pop <- p_data$pop 

#Calculate starting values for the Gaussian model 
  start <- rTPC::get_start_vals(d$temp, d$rate, "gaussian_1987")

#Assign these values to the different parameters 
  start ["rmax"]<- 0.933 
  start["topt"] <- 11.25
  start["a"] <- 12.5

#Run a standard Non-Linear model (no random effects)
nls_mod <- nls(rate ~ rTPC::gaussian_1987(temp,  rmax, topt, a),
               data = d, 
               start = start)

#Establish predictions for the nls model 
  d$pred <- predict(nls_mod)

#Run a non-linear mixed effects model using nlme package
  mod <- nlme(rate ~ rTPC::gaussian_1987(temp,  rmax, topt, a),
            data = d,
            fixed = rmax + topt + a ~ 1,
            random = rmax + topt + a ~ 1,
            group = ~pop,
            start = coef(nls_mod),
            method = "ML",
            control = nlmeControl(pnlsTol = 0.2, msMaxIter =  500))


#Establish a vector with a range of temperatures
  xtemp <- seq(-5, 30, length.out = 10000)

#Establish a vector of matching lengths for each population (1 in this case)
  pop <- seq(1,1, length.out = 10000)

#Create a dataframe with these vectors 
  df1 <- data.frame(xtemp,pop)

#Match the column names with the nlme model 
  colnames(df1)[1] <-"temp"

#Predict the model over this data frame
df1$pred <- predict(mod, df1)

#Repeat the process for each population, altering the population to predict over each time 

#Population 2
pop2 <- seq(2,2, length.out = 10000)
df2 <- data.frame(xtemp,pop2)
colnames(df2)[1] <- "temp"
colnames(df2)[2] <- "pop"
df2$pred <- predict(mod, df2)

#Population 3 
pop3 <- seq(3,3, length.out = 10000)
df3 <- data.frame(xtemp,pop3)
colnames(df3)[1] <- "temp"
colnames(df3)[2] <- "pop"
df3$pred <- predict(mod, df3)

#Population 4 
pop4 <- seq(4,4, length.out = 10000)
df4 <- data.frame(xtemp,pop4)
colnames(df4)[1] <- "temp"
colnames(df4)[2] <- "pop"
df4$pred <- predict(mod, df4)

#Population 5
pop5 <- seq(5,5, length.out = 10000)
df5 <- data.frame(xtemp,pop5)
colnames(df5)[1] <- "temp"
colnames(df5)[2] <- "pop"
df5$pred <- predict(mod, df5)

#Population 6 
pop6 <- seq(6,6, length.out = 10000)
df6 <- data.frame(xtemp,pop6)
colnames(df6)[1] <- "temp"
colnames(df6)[2] <- "pop"
df6$pred <- predict(mod, df6)

#Population 7 
pop7 <- seq(7,7, length.out = 10000)
df7 <- data.frame(xtemp,pop7)
colnames(df7)[1] <- "temp"
colnames(df7)[2] <- "pop"
df7$pred <-predict(mod, df7)

#Population 8
pop8 <- seq(8,8, length.out = 10000)
df8 <- data.frame(xtemp,pop8)
colnames(df8)[1] <- "temp"
colnames(df8)[2] <- "pop"
df8$pred <- predict(mod, df8)

#Population 9
pop9 <- seq(9,9, length.out = 10000)
df9 <- data.frame(xtemp,pop9)
colnames(df9)[1] <- "temp"
colnames(df9)[2] <- "pop"
df9$pred <- predict(mod, df9)

#Population 10 
pop10 <- seq(10,10, length.out = 10000)
df10 <- data.frame(xtemp,pop10)
colnames(df10)[1] <- "temp"
colnames(df10)[2] <- "pop"
df10$pred <- predict(mod, df10)

#Create a plot to compare each individual TPC for populations
pop_plot <- ggplot () + 
  geom_line(aes(temp, pred), df1, col = "cornflowerblue", size = 1) +
  geom_line(aes(temp, pred), df2, col = "darkgreen", size = 1) + 
  geom_line(aes(temp, pred), df3, col = "gold3", size = 1) + 
  geom_line(aes(temp, pred), df4, col = "mediumblue", size = 1) + 
  geom_line(aes(temp, pred), df5, col = "purple", size = 1 ) +
  geom_line(aes(temp, pred), df6, col = "lawngreen", size = 1) + 
  geom_line(aes(temp, pred), df7, col = "maroon4", size = 1) + 
  geom_line(aes(temp, pred), df8, col = "orange", size = 1) + 
  geom_line(aes(temp, pred), df9, col = "red3", size = 1 ) + 
  geom_line(aes(temp, pred), df10, col = "yellow1") + 
  geom_point(aes(min_t, perc_q), p_summary) +
  geom_errorbar(aes(x = min_t, y = perc_q, 
                    ymin = perc_q - sd, ymax = perc_q + sd),p_summary, width = 0.2, size = 1) +
  ylim(0, 0.8) +
  xlim(-2.5, 25) +
  labs(x = 'Nocturnal Minimum Temperature (ºC)',
       y = 'Questing Proportion') + 
  theme(plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm")) +
  theme(axis.title.y = element_text(vjust=4, size = 27)) + 
  theme(axis.title.x = element_text(vjust=-4, size = 27)) + 
  theme(axis.text = element_text(size = 27, colour = "black")) + 
  theme(legend.position="none")+
  theme(panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill=NA, size=1.5)) +
  theme(plot.title = element_text(hjust = 0.5, vjust = 3, size = 35)) +
  theme(axis.line = element_line(colour = "black"),
        panel.background = element_blank()) 

#Check Plot 
  pop_plot

#Calculate the optimal temperature and associated maximum questing proportion for 
  #each population
  df1_opt <- df1 %>% filter(pred == max(pred))
  df2_opt <- df2 %>% filter(pred == max(pred))
  df3_opt <- df3 %>% filter(pred == max(pred))
  df4_opt <- df4 %>% filter(pred == max(pred))
  df5_opt <- df5 %>% filter(pred == max(pred))
  df6_opt <- df6 %>% filter(pred == max(pred))
  df7_opt <- df7 %>% filter(pred == max(pred))
  df8_opt <- df8 %>% filter(pred == max(pred))
  df9_opt <- df9 %>% filter(pred == max(pred))
  df10_opt <- df10 %>% filter(pred == max(pred))

#Combine these values into a single dataframe
  pop_opts <- rbind(df1_opt, df2_opt, df3_opt, df4_opt, df5_opt, df6_opt, df7_opt,
                    df8_opt, df9_opt, df10_opt)

#Round the values to three decimal places and rename/reorder the columns 
  pop_opts$temp <- round(pop_opts$temp, digits = 3)
  pop_opts$pred <- round(pop_opts$pred, digits = 3)
  pop_opts <- pop_opts[, c(2,1,3)]
  colnames(pop_opts)[1] <- "Population"
  colnames(pop_opts)[2] <- "Optimum 
Temperature(C)"
  colnames(pop_opts)[3] <- "Questing 
Proportion"

#Create a table with matching colours to those in the plot 
  table.p <- ggtexttable(pop_opts, rows = NULL, theme = ttheme(colnames.style = colnames_style(size = 17), padding = unit(c(3,3), "mm"), tbody.style = tbody_style(size = 27)))
  table.p <- table_cell_bg(table.p, row = 2, column = 1:3, fill = "cornflowerblue", alpha = 0.7)
  table.p <- table_cell_bg(table.p, row = 3, column = 1:3, fill = "grey", alpha = 0.7)
  table.p <- table_cell_bg(table.p, row = 4, column = 1:3, fill = "gold3", alpha = 0.7)
  table.p <- table_cell_bg(table.p, row = 5, column = 1:3, fill = "mediumblue", alpha = 0.7)
  table.p <- table_cell_bg(table.p, row = 6, column = 1:3, fill = "purple", alpha = 0.7)
  table.p <- table_cell_bg(table.p, row = 7, column = 1:3, fill = "lawngreen", alpha = 0.7)
  table.p <- table_cell_bg(table.p, row = 8, column = 1:3, fill = "maroon4", alpha = 0.7)
  table.p <- table_cell_bg(table.p, row = 9, column = 1:3, fill = "orange", alpha = 0.7)
  table.p <- table_cell_bg(table.p, row = 10, column = 1:3, fill = "red3", alpha = 0.7)
  table.p <- table_cell_bg(table.p, row = 11, column = 1:3, fill = "yellow1", alpha = 0.7)
  table.p <- table_cell_font(table.p, row = 2:11, column = 1:3, size = 20)
  table.p <- thead_add_border(table.p,
                              from.row = 1,
                              to.row = tab_nrow(table.p),
                              from.column = 1,
                              to.column = tab_ncol(table.p),
                              linetype = 1,
                              linewidth = 1, linecolor = "black")
#Check the table 
  table.p

#Combine the two plots, using the table as a key for use in the report 
  ggarrange(pop_plot, table.p, widths = c(2,1), heights = c(1, 20))
  
  