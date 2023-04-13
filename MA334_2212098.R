# import libraries
library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)

# set the working directory
setwd("S:/Statistics for Data Science/MA334_Data analysis and Statistics with R/Assignment/MA334_2212098/Biodiversity_Aanalysis")

# read the data file
data_v3 <- read.csv("proportional_species_richness_V3.csv")

# get the count of missing values in each column
(sapply(data_v3, function(x) sum(is.na(x))))

# convert categorical variable to factor
data_v3$dominantLandClass <- as.factor(data_v3$dominantLandClass)

data_v3$period <- as.factor(data_v3$period)

# select assigned taxonomic groups
my_data <- data_v3 %>%
  select(c(2,4,5,7,9,10,12:17))

# summarise and print first few records of the data
head(my_data)
str(my_data)

my_data %>%
  select(c(1:7)) %>%
  summary()

# compute the mean ecological status for the 7 chosen species
mean_species_7 <- rowMeans(my_data[1:7], na.rm = TRUE)
my_data$eco_status_7 <- mean_species_7

# check the correlation of all 7 species with the computed mean eco score
cor(my_data[1:7], my_data$eco_status_7)

# compute the mean of only the priority species (high corr with eco score)
mean_priority_species <- rowMeans(my_data[c(1,3:6)], na.rm = T)

# compute the mean of remaining species (low corr with eco score)
mean_rem_species <- rowMeans(my_data[c(2,7)], na.rm = T)

# check the correlation for priority and remaining species 
# with the total eco score
cor(mean_priority_species, mean_species_7)
cor(mean_rem_species, mean_species_7)

my_data$eco_status_priority <- mean_priority_species
my_data$eco_status_rem <- mean_rem_species

priority_species_eco_change <- my_data %>%
  group_by(Easting, Northing, period) %>%
  summarise(species_mean = mean(eco_status_priority), .groups = 'drop') %>%
  pivot_wider(names_from = period, values_from = species_mean, values_fill = 0) %>%
  mutate(difference_eco = Y00 - Y70) %>%
  arrange(difference_eco)

# low_eco_status_priority_species <- priority_species_eco_change %>%
#   filter(difference_eco <= 0) %>%
#   count()
# 
# high_eco_status_priority_species <- priority_species_eco_change %>%
#   filter(difference_eco > 0) %>%
#   count()
# 
# 
# rem_species_eco_change <- my_data %>%
#   group_by(Easting, Northing, period) %>%
#   summarise(species_mean = mean(eco_status_rem), .groups = 'drop') %>%
#   pivot_wider(names_from = period, values_from = species_mean, values_fill = 0) %>%
#   mutate(difference_eco = Y00 - Y70) %>%
#   arrange(difference_eco)

# low_eco_status_rem_species <- rem_species_eco_change %>%
#   filter(difference_eco <= 0) %>%
#   count()
# 
# high_eco_status_rem_species <- rem_species_eco_change %>%
#   filter(difference_eco > 0) %>%
#   count()

# all_species_eco_change <- my_data %>%
#   group_by(Easting, Northing, period) %>%
#   summarise(species_mean = mean(eco_status_7), .groups = 'drop') %>%
#   pivot_wider(names_from = period, values_from = species_mean, values_fill = 0) %>%
#   mutate(difference_eco = Y00 - Y70) %>%
#   arrange(difference_eco)

# low_eco_status_all_species <- all_species_eco_change %>%
#   filter(difference_eco <= 0) %>%
#   count()
# 
# high_eco_status_all_species <- all_species_eco_change %>%
#   filter(difference_eco > 0) %>%
#   count()


t.test(priority_species_eco_change$difference_eco,
       alternative = "greater",
       mu = 0,
       conf.level = 0.95)

eco_change_zone_7 <- my_data %>%
  group_by(dominantLandClass, period) %>%
  summarise(species_mean = mean(eco_status_priority), .groups = 'drop') %>%
  pivot_wider(names_from = period, values_from = species_mean, values_fill = 0) %>%
  mutate(diff = Y00 - Y70) %>%
  arrange(diff) %>% print(n=45)

t.test(eco_change_zone_7$diff,
       alternative = "greater",
       mu = 0,
       conf.level = 0.99)


# Linear regression

# function to plot the linear regression results for each period
regression_plot <- function(lm_model) {
  ggplot(lm_model$model,
         aes_string(x = names(lm_model$model)[2], y = names(lm_model$model)[1])) +
    geom_point() +
    geom_smooth(method = lm, col = "#159895", se = F) +
    labs(title = paste(
      "Period = ", lm_model$terms[[2]][[2]], ",",
      "Intercept = ", round(lm_model$coefficients[[1]], 3), ",",
      "Slope = ", round(lm_model$coefficients[[2]], 3), ",",
      "p-value = ", round(summary(lm_model)$coefficients[2,4], 3), ",",
      "Adjusted R-squared = ", round(summary(lm_model_00)$adj.r.squared, 4)
    ),
    x = "Ecological score of 7 species group",
    y = "Ecological score of 11 species group"
    ) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black")
          )
}


# function to plot the residuals of the regression model
residuals_plot <- function(lm_model) {
  ggplot(lm_model, aes(sample = lm_model$residuals)) +
    stat_qq() +
    stat_qq_line(col = "red", lwd = 1) +
    labs(x = "Normal distribution",
         y = "Residual's distribution",
         title = paste("Distribution of residuals of lm() for ", lm_model$terms[[2]][[2]])
         ) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black")
    )
}

# filter for Y00 period
period_00 <- my_data %>%
  filter(period=="Y00")

# plot the correlation between eco_status_7 and ecologicalStatus
period_00 %>%
  ggplot(aes(x=eco_status_7, y=ecologicalStatus)) +
  geom_point() +
  geom_smooth(method = lm, col = "red", se = F) +
  theme_bw()

# fit the linear regression model
lm_model_00 <- lm(period_00$ecologicalStatus~period_00$eco_status_7)

# plot the lm() results
regression_plot(lm_model = lm_model_00)

# qq-plot for the residuals of the lm()
residuals_plot(lm_model = lm_model_00)


# filter for Y70 period
period_70 <- my_data %>%
  filter(period=="Y70")

# plot the correlation between eco_status_7 and ecologicalStatus
period_70 %>%
  ggplot(aes(x=eco_status_7, y=ecologicalStatus)) +
  geom_point() +
  geom_smooth(method = lm, col = "red", se = F) +
  theme_bw()

# fit the linear regression model
lm_model_70 <- lm(period_70$ecologicalStatus~period_70$eco_status_7)

# plot the lm() results
regression_plot(lm_model = lm_model_70)

# qq-plot for the residuals of the lm()
residuals_plot(lm_model = lm_model_70)

# multiple linear regression


