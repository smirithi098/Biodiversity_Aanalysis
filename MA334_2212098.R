# import libraries
library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
library(AICcmodavg)
library(gridExtra)
library(reshape2)

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

# Data Exploration

# Exploration - 1: Correlation between the species

# compute correlation between all the 7 species
cor_species <- melt(cor(my_data[1:7], my_data[1:7]))

# plot the heatmap
cor_species %>%
  ggplot(aes(x=Var1, y=Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient(high = "#C85C8E", low = "#FFBABA") +
  labs(title = "Correlation between 7 species" , x = "Species", y ="Species") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

# Exploration - 2: Mean distribution of original 11 species and
# allocated 7 species across both periods

# plot for the mean of 7 species
plot_7 <- my_data %>%
  ggplot(aes(x=eco_status_7)) +
  geom_histogram(bins = 45,
                 aes(y=after_stat(density)),
                 colour = "black", fill = "lightgrey") +
  geom_vline(aes(xintercept=mean(eco_status_7)), 
             linetype = "dashed", size = 0.6) +
  geom_density(lwd=0.8, fill = "#FF6666", alpha = 0.18) +
  facet_wrap(vars(period), scales = "free", ncol = 2) +
  labs(title = "Mean distribution of 7 species for both periods",
       x = "Mean of 7 species") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank()
  )

# plot for the mean of 11 species
plot_11 <- my_data %>%
  ggplot(aes(x=ecologicalStatus)) +
  geom_histogram(bins = 45,
                 aes(y=after_stat(density)),
                 colour = "black", fill = "lightgrey") +
  geom_vline(aes(xintercept=mean(ecologicalStatus)), 
             linetype = "dashed", size = 0.6) +
  geom_density(lwd=0.8, fill = "#FF6666", alpha = 0.18) +
  facet_wrap(vars(period), scales = "free", ncol = 2) +
  labs(title = "Mean distribution of 11 species for both periods",
       x = "Mean of 11 species") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank()
  )

# arrange both plots to compare the distributions
grid.arrange(plot_7, plot_11, nrow = 2)


# Hypothesis testing

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

# Hypothesis test - 1

priority_species_eco_change <- my_data %>%
  group_by(Easting, Northing, period) %>%
  summarise(species_mean = mean(eco_status_priority), .groups = 'drop') %>%
  pivot_wider(names_from = period, values_from = species_mean, values_fill = 0) %>%
  mutate(difference_eco = Y00 - Y70) %>%
  arrange(difference_eco)


t.test(priority_species_eco_change$difference_eco,
       alternative = "greater",
       mu = 0,
       conf.level = 0.95)

my_data %>%
  filter(str_detect(dominantLandClass, 'e')) %>%
  ggplot(aes(x=Easting, y=eco_status_7, colour = period)) +
  geom_point(show.legend = F) +
  geom_smooth(method = lm, se=F, col = "darkgrey") +
  facet_wrap(vars(dominantLandClass), ncol = 3)

# Hypothesis test - 2

england_species_7 <- my_data %>%
  filter(str_detect(dominantLandClass, 'e')) %>%
  group_by(dominantLandClass, period) %>%
  summarise(species_mean=mean(eco_status_7), .groups = 'drop') %>%
  pivot_wider(names_from = period, values_from = species_mean, values_fill = 0) %>%
  mutate(diff=Y00-Y70) %>% print(n=21)

england_species_11 <- my_data %>%
  filter(str_detect(dominantLandClass, 'e')) %>%
  group_by(dominantLandClass, period) %>%
  summarise(species_mean=mean(ecologicalStatus), .groups = 'drop') %>%
  pivot_wider(names_from = period, values_from = species_mean, values_fill = 0) %>%
  mutate(diff=Y00-Y70) %>% print(n=21)

t.test(england_species_7$diff,
       england_species_11$diff,
       alternative = "greater",
       mu = 0,
       conf.level = 0.95)

# Linear regression

# function to plot the linear regression results for each period
regression_plot <- function(data, lm_model) {
  ggplot(data,
         aes(x = predicted, y = observed)) +
    geom_point() +
    geom_smooth(method = lm, col = "#159895", se = F) +
    labs(title = paste(
      "Period = ", lm_model$terms[[2]][[2]], ",",
      "Intercept = ", round(lm_model$coefficients[[1]], 3), ",",
      "Slope = ", round(lm_model$coefficients[[2]], 3), ",",
      "p-value = ", round(summary(lm_model)$coefficients[2,4], 3), ",",
      "Adjusted R-squared = ", round(summary(lm_model_00)$adj.r.squared, 4)
    ),
    x = "Observed Values (Y)",
    y = "Predicted Values (Yhat)"
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
lm_model_00_df <- data.frame(predicted=fitted(lm_model_00), 
                             observed=period_00$ecologicalStatus)

# plot the lm() results
regression_plot(data = lm_model_00_df, lm_model = lm_model_00)

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
lm_model_70_df <- data.frame(predicted=fitted(lm_model_70), 
                             observed=period_70$ecologicalStatus)

# plot the lm() results
regression_plot(data = lm_model_70_df, lm_model = lm_model_70)

# qq-plot for the residuals of the lm()
residuals_plot(lm_model = lm_model_70)

# multiple linear regression

# select the remaining 4 species and compute the mean ecological score
species_4 <- data_v3 %>% select(c(3, 6, 8, 11))

my_data <- my_data %>%
  mutate(eco_status_4 = rowMeans(species_4[1:4], na.rm = T))

# function to plot the correlation for the prediction on test data by all mlr models
mlr_model_plot <- function(test_pred, model) {
  subtitle <- paste(str_remove_all(as.character(model$terms[[3]][2]), '()'),
                    "+", 
                    str_remove_all(as.character(model$terms[[3]][3]), '()'))
  ggplot(test_pred, aes(x=predicted, y=observed)) +
    geom_point() +
    geom_smooth(method = lm, se=F, col = "#EB455F") +
    labs(title = paste("Correlation = ", 
                       round(cor(test_pred$predicted, test_pred$observed), 4)),
         subtitle = paste("Species = ", subtitle),
         x = "Predicted test data (Yhat)",
         y = "Observed test data (Y)") +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black")
    )
}

# split the data into training and test data
train_set <- sample(1:nrow(my_data), 0.8*nrow(my_data))
training_data <- my_data[train_set, ]
testing_data <- my_data[-train_set, ]

# Model - 1

# run the multiple linear regression model with all 7 species
mlr_model_all <- lm(eco_status_4~.,
                data = training_data[c("Bees", "Bryophytes", "Butterflies", 
                                       "Hoverflies", "Ladybirds", "Macromoths",
                                       "Vascular_plants", "eco_status_4")],
                na.action=na.omit,
                y=TRUE)

summary(mlr_model_all)

# check the correlation between the predicted and observed values of the response variable
cor(mlr_model_all$fitted.values, mlr_model_all$y)

# run model prediction for test data
mlr_model_pred_all <- predict(mlr_model_all, testing_data)

# create a data frame containing the test prediction and original variables
pred_all <- data.frame(predicted=mlr_model_pred_all, observed=testing_data$eco_status_4)

# Model - 2

# run the mlr model with only species with p-value < 2e-16
mlr_model_4 <- lm(eco_status_4~.,
                  data = training_data[c("Bryophytes", "Hoverflies", "Ladybirds",
                                         "Vascular_plants", "eco_status_4")],
                  na.action=na.omit,
                  y=TRUE)

summary(mlr_model_4)

# check the correlation between the predicted and observed values of the response variable
cor(mlr_model_4$fitted.values, mlr_model_4$y)

# run model prediction for test data
mlr_model_pred_4 <- predict(mlr_model_4, testing_data)

# create a data frame containing the test prediction and original variables
pred_4 <- data.frame(predicted=mlr_model_pred_4, observed=testing_data$eco_status_4)

# Model - 3

# run the multiple linear regression model with remaining 3 species
mlr_model_3 <- lm(eco_status_4~.,
                    data = training_data[c("Bees", "Butterflies", "Macromoths",
                                           "eco_status_4")],
                    na.action=na.omit,
                    y=TRUE)

summary(mlr_model_3)

# check the correlation between the predicted and observed values of the response variable
cor(mlr_model_3$fitted.values, mlr_model_3$y)

# run model prediction for test data
mlr_model_pred_3 <- predict(mlr_model_3, testing_data)

# create a data frame containing the test prediction and original variables
pred_3 <- data.frame(predicted=mlr_model_pred_3, observed=testing_data$eco_status_4)

# Model - 4

# run the multiple linear regression model with 5 species with optimal p-values
mlr_model_5 <- lm(eco_status_4~.,
                    data = training_data[c("Butterflies", "Hoverflies", 
                                           "Ladybirds", "Macromoths", "Vascular_plants",
                                           "eco_status_4")],
                    na.action=na.omit,
                    y=TRUE)

summary(mlr_model_5)

# check the correlation between the predicted and observed values of the response variable
cor(mlr_model_5$fitted.values, mlr_model_5$y)

# run model prediction for test data
mlr_model_pred_5 <- predict(mlr_model_5, testing_data)

# create a data frame containing the test prediction and original variables
pred_5 <- data.frame(predicted=mlr_model_pred_5, observed=testing_data$eco_status_4)

# Model - 5

# run the multiple linear regression model with 2 species with high correlation
mlr_model_2 <- lm(eco_status_4~.,
                    data = training_data[c("Hoverflies", "Ladybirds", "eco_status_4")],
                    na.action=na.omit,
                    y=TRUE)

summary(mlr_model_2)

# check the correlation between the predicted and observed values of the response variable
cor(mlr_model_2$fitted.values, mlr_model_2$y)

# run model prediction for test data
mlr_model_pred_2 <- predict(mlr_model_2, testing_data)

# create a data frame containing the test prediction and original variables
pred_2 <- data.frame(predicted=mlr_model_pred_2, observed=testing_data$eco_status_4)

# plot the scatter plot for the test prediction of all mlr models
plot_7 <- mlr_model_plot(pred_all, mlr_model_all)
plot_5 <- mlr_model_plot(pred_5, mlr_model_5)
plot_4 <- mlr_model_plot(pred_4, mlr_model_4)
plot_3 <- mlr_model_plot(pred_3, mlr_model_3)
plot_2 <- mlr_model_plot(pred_2, mlr_model_2)

# arrange all the corr plots in a grid
grid.arrange(plot_7, plot_5, plot_4, plot_3, plot_2, ncol = 2)


# compute the AIC for all the mlr models and select the best fit one

models <- list(mlr_model_all, mlr_model_5, mlr_model_4, mlr_model_3, mlr_model_2)
model_names <- c('7 species', '5 species', '4 species', '3 species', '2 species')

aictab(cand.set = models, modnames = model_names) 
# the mlr model with all 7 species as predictors has the lowest AIC value 
# thereby the best fit model to predict BD4

# Open analysis

my_data %>%
  filter(period=="Y70") %>%
  group_by(dominantLandClass) %>%
  ggplot(aes(x=eco_status_7)) +
  geom_histogram(aes(y=after_stat(density)), colour = "black", fill = "lightgrey") +
  geom_density(lwd=0.8, fill = "#FF6666", alpha = 0.18) +
  facet_wrap(vars(dominantLandClass), scales = "free") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank()
  )
