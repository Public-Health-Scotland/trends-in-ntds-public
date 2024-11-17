###### TRENDS IN NEURAL TUBE DEFECTS (NTD) R CODE SCRIPT

#This script conducts the analyses looking at the time trends in NTDs

# ---- Section 1: Housekeeping ----

# clear the global environment
rm(list = ls())

# load required packages
require(tidyverse)
require(ggplot2)
require(lubridate)
require(dplyr)
require(broom)
require(VGAM) # for VGLM model
require(rms) # for restricted cubic splines 'rcs'
require(glmmTMB) # for Poisson COM model
require(segmented)

# Read in prepared data require for descriptive analyses 
df <- read.csv("/conf/Congenital-Descriptive-Epi/Trends in NTDs/Data/analysis_ready.csv")

# ---- Section 2: Prepare data to look at trends ----

# Calculate total births, and total number of outcomes for each year

#Calculate total number of NTDs and create a year variable coded from 0-21 
#for use with some models
df_sum <- df %>%
  group_by(year) %>%
  summarise(
    births_total = sum(count_all),
    n_ntd = sum(NTD_all),
    n_anen = sum(anen),
    n_sb = sum(SB),
    n_ntgng = sum(NTD_no_gen)) %>%
  ungroup() %>%
  mutate(year2=year-2000)

# ---- Section 3: Time trends in all NTDs ----

# First use Poisson model for association of total birth prevalence with year 
# this model assumes a linear trend over time
fit_p_l<-vglm(n_ntd ~ year, poissonff(bred=TRUE), offset=log(births_total), 
              data=df_sum)
summary(fit_p_l)
# Produce AIC to compare Poisson models
AIC(fit_p_l) 
# check dispersion parameter for this model
dp_p_l <-  deviance(fit_p_l)/df.residual(fit_p_l, type = "vlm")
dp_p_l

## Then use alternative Poisson model for association of total birth prevalence with year 
# this model assumes a non-linear trend over time using restricted cubic splines
# define knots for this model based on trends observed on plot
knots <- c(2000, 2007, 2014, 2021)
fit_p_nl<-vglm(n_ntd ~ rcs(year,knots), poissonff(bred=TRUE), offset=log(births_total), 
               data=df_sum)
summary(fit_p_nl)
# Produce AIC to compare Poisson models
AIC(fit_p_nl) 
# check dispersion parameter for this model
dp_p_nl <-  deviance(fit_p_nl)/df.residual(fit_p_nl, type = "vlm")
dp_p_nl

## Choosing the first model (fit_p_l) based on AIC criteria and after checking dispersion parameter

# calculate prevalence rate ratio (PRR) by exponentiating the d$time_period coefficient estimate
# and then create a dataframe to read exponentiated values from model
# Extract coefficients
coefficients <- coef(fit_p_l)
# Extract standard errors
std_errors <- sqrt(diag(vcov(fit_p_l)))
# Calculate z-values
z_values <- coefficients / std_errors
# Calculate p-values
p_values <- 2 * (1 - pnorm(abs(z_values)))
# Create a tidy data frame
tidy_fit <- data.frame(
  term = names(coefficients),
  estimate = coefficients,
  std.error = std_errors,
  z.value = z_values,
  p.value = p_values
)
#exponentiate estimates and calculate lower and upper confidence intervals
tidy_fit_all <- tidy_fit %>%
  mutate(PRR = round(exp(estimate), 3),
         LCI = round(exp(estimate - ((qnorm(0.975))*std.error)),3),
         UCI = round(exp(estimate + ((qnorm(0.975))*std.error)),3)) %>%
  subset(term=="year") %>%
  mutate(outcome = "All NTD") %>%
  dplyr::select(term, outcome, PRR, LCI, UCI, p.value) 

## For selected Poisson model (fit_p_l in this case)
# Predict values
predicted_values <- predict(fit_p_l, newdata = df_sum, type = "link", se.fit = TRUE)

# Extract predicted values and standard errors
predicted_means <- predicted_values$fit
standard_errors <- predicted_values$se.fit

# Calculate confidence intervals
z_value <- qnorm(0.975)  # 1.96 for 95% confidence interval
lower_bound <- exp(predicted_means - z_value * standard_errors)
upper_bound <- exp(predicted_means + z_value * standard_errors)

df_sum_exp <- df_sum %>%
mutate(predicted_mean1 = (exp(predicted_means) / births_total) * 10000,
       predicted_lci1 = (lower_bound / births_total) * 10000,
       predicted_uci1 = (upper_bound / births_total) * 10000,
       prevs1 = (n_ntd / births_total) * 10000)

# ---- Section 4: Time trends in anencephaly ----

# First use Poisson model for association of total birth prevalence and year which assumes a linear trend over time
fit_p_l<-vglm(n_anen ~ year, poissonff(bred=TRUE), offset=log(births_total), 
              data=df_sum)
summary(fit_p_l)
# check AIC
AIC(fit_p_l) 
# check dispersion parameter
dp_p_l <-  deviance(fit_p_l)/df.residual(fit_p_l, type = "vlm")
dp_p_l

# Then use Poisson model for association of total birth prevalence and year which assumes a non-linear trend over time
# set knots based on plot
knots <- c(2000, 2008, 2021)
fit_p_nl<-vglm(n_anen ~ rcs(year,knots), poissonff(bred=TRUE), offset=log(births_total), 
               data=df_sum)
summary(fit_p_nl)
# check AIC
AIC(fit_p_nl)
# check dispersion parameter
dp_p_nl <-  deviance(fit_p_nl)/df.residual(fit_p_nl, type = "vlm")
dp_p_nl

# Both models are overdispersed so we will use COM-Poisson

# COM-Poisson model for association of total birth prevalence with year assuming a linear trend over time
fit_pc_l <- glmmTMB(n_anen ~ year2 + offset(log(births_total)), family=compois, data=df_sum)
summary(fit_pc_l)
# check AIC
AIC(fit_pc_l) 

# COM-Poisson model for association of total birth prevalence with year assuming a non-linear trend over time
# set knots based on plot
knots <- c(0, 13, 21)
fit_pc_nl <- glmmTMB(n_anen ~ rcs(year2,knots) + offset(log(births_total)), family=compois, data=df_sum)
summary(fit_pc_nl)
# check AIC
AIC(fit_pc_nl) 

# Piecewise regression models, capturing trend in anencephaly within each time period

#### # Period 1: 2000-2012 ####
d1 <- df_sum %>%
  subset(year2<13)
  
fit_pc_l_p1 <- glmmTMB(n_anen ~ year2 + offset(log(births_total)), family=compois, data=d1)
summary(fit_pc_l_p1)

coefficients <- as.numeric(fixef(fit_pc_l_p1)$cond)
std_errors <- sqrt(diag(vcov(fit_pc_l_p1)$cond))
z_values <- coefficients / std_errors
p_values <- 2 * (1 - pnorm(abs(z_values)))
exp_coefficients <- exp(coefficients)
lower_ci <- exp(coefficients - 1.96 * std_errors)
upper_ci <- exp(coefficients + 1.96 * std_errors)

# Create a tidy data frame
tidy_fit <- data.frame(
  term = names(fixef(fit_pc_l_p1)$cond),
  estimate = coefficients,
  PRR = exp_coefficients,
  std.error = std_errors,
  lower_ci = lower_ci,
  upper_ci = upper_ci,
  z.value = z_values,
  p.value = p_values
)

#exponentiate estimates and calculate lower and upper confidence intervals
tidy_fit_anen <- tidy_fit %>%
  mutate(PRR = round(exp(estimate), 3),
         LCI = round(exp(estimate - ((qnorm(0.975))*std.error)),3),
         UCI = round(exp(estimate + ((qnorm(0.975))*std.error)),3)) %>%
  subset(term=="year2") %>%
  mutate(outcome = "Anencephaly") %>%
  mutate(term = if_else(term == "year2", "year 2000-12", "NA")) %>%
  dplyr::select(term, outcome, PRR, LCI, UCI, p.value) 

#### # Period 2: 2013-2021 ####
d2 <- df_sum%>%
  subset(year2>=13) %>%
  mutate(row = row_number())

fit_pc_l_p2 <- glmmTMB(n_anen ~ row + offset(log(births_total)), family=compois, data=d2)
summary(fit_pc_l_p2)

coefficients <- as.numeric(fixef(fit_pc_l_p2)$cond)
std_errors <- sqrt(diag(vcov(fit_pc_l_p2)$cond))
z_values <- coefficients / std_errors
p_values <- 2 * (1 - pnorm(abs(z_values)))
exp_coefficients <- exp(coefficients)
lower_ci <- exp(coefficients - 1.96 * std_errors)
upper_ci <- exp(coefficients + 1.96 * std_errors)

# Create a tidy data frame
tidy_fit <- data.frame(
  term = names(fixef(fit_pc_l_p2)$cond),
  estimate = coefficients,
  PRR = exp_coefficients,
  std.error = std_errors,
  lower_ci = lower_ci,
  upper_ci = upper_ci,
  z.value = z_values,
  p.value = p_values
)
#exponentiate estimates and calculate lower and upper confidence intervals
tidy_fit_anen_2 <- tidy_fit %>%
  mutate(PRR = round(exp(estimate), 3),
         LCI = round(exp(estimate - ((qnorm(0.975))*std.error)),3),
         UCI = round(exp(estimate + ((qnorm(0.975))*std.error)),3)) %>%
  subset(term=="row") %>%
  mutate(outcome = "Anencephaly") %>%
  mutate(term = if_else(term == "row", "year 2013-21", "NA")) %>%
  dplyr::select(term, outcome, PRR, LCI, UCI, p.value) 

### CALCULATE PREDICTED VALUES FOR GRAPHING FROM MODELS (Using splines)
predicted_values <- predict(fit_pc_nl, newdata = df_sum_exp, type = "link", se.fit = TRUE)

# Extract predicted values and standard errors
predicted_means <- predicted_values$fit
standard_errors <- predicted_values$se.fit

# Calculate confidence intervals
z_value <- qnorm(0.975)  # 1.96 for 95% confidence interval
lower_bound <- exp(predicted_means - z_value * standard_errors)
upper_bound <- exp(predicted_means + z_value * standard_errors)

df_sum_exp <- df_sum_exp %>%
  mutate(predicted_mean2 = (exp(predicted_means) / births_total) * 10000,
         predicted_lci2 = (lower_bound / births_total) * 10000,
         predicted_uci2 = (upper_bound / births_total) * 10000,
         prevs2 = (n_anen / births_total) * 10000)

# ---- Section 5: Time trends in spina bifida ----

# First use model which assumes a linear time trend
fit_p_l <-vglm(n_sb ~ year, poissonff(bred=TRUE), offset=log(births_total), 
              data=df_sum)
summary(fit_p_l)
# check AIC
AIC(fit_p_l) 
# check dispersion parameter
dp_p_l <-  deviance(fit_p_l)/df.residual(fit_p_l, type = "vlm")
dp_p_l

# Then use Poisson model which assumes non-linear trend over time with selected knots
knots <- c(2000, 2008, 2021)
fit_p_nl<-vglm(n_sb ~ rcs(year,knots), poissonff(bred=TRUE), offset=log(births_total), 
               data=df_sum)
summary(fit_p_nl)
# check AIC
AIC(fit_p_nl) 

# check dispersion parameter
dp_p_nl <-  deviance(fit_p_nl)/df.residual(fit_p_nl, type = "vlm")
dp_p_nl

# As these models are overdispersed, we will use a COM-Poisson model

# COM-Poisson model for year assuming a linear time trend
fit_pc_l <- glmmTMB(n_sb ~ year2 + offset(log(births_total)), family=compois, data=df_sum)
summary(fit_pc_l)
# check AIC
AIC(fit_pc_l) 

# COM-Poisson model for year assuming a non-linear trend over time with selected knots
knots <- c(0, 8, 21)
fit_pc_nl <- glmmTMB(n_sb ~ rcs(year2, knots) + offset(log(births_total)), family=compois, data=df_sum)
summary(fit_pc_nl)
# check AIC
AIC(fit_pc_nl)

# select the COM-Poisson model assuming a linear trend over time based on AIC value

# Create a dataframe to read exponentiated values from model
# Extract coefficients (conditional model) as numeric
coefficients <- as.numeric(fixef(fit_pc_l)$cond)

# Extract standard errors (conditional model) as numeric
std_errors <- sqrt(diag(vcov(fit_pc_l)$cond))

# Calculate z-values
z_values <- coefficients / std_errors

# Calculate p-values
p_values <- 2 * (1 - pnorm(abs(z_values)))

# Exponentiate coefficients to get PRRs (prevalence rate ratios)
exp_coefficients <- exp(coefficients)

# Calculate 95% confidence intervals
lower_ci <- exp(coefficients - 1.96 * std_errors)
upper_ci <- exp(coefficients + 1.96 * std_errors)

# Create a tidy data frame
tidy_fit <- data.frame(
  term = names(fixef(fit_pc_l)$cond),
  estimate = coefficients,
  PRR = exp_coefficients,
  std.error = std_errors,
  lower_ci = lower_ci,
  upper_ci = upper_ci,
  z.value = z_values,
  p.value = p_values
)

tidy_fit_sb <- tidy_fit %>%
  mutate(PRR = round(exp(estimate), 3),
         LCI = round(exp(estimate - ((qnorm(0.975))*std.error)),3),
         UCI = round(exp(estimate + ((qnorm(0.975))*std.error)),3)) %>%
  subset(term=="year2") %>%
  mutate(outcome = "spina bifida") %>%
  dplyr::select(term, outcome, PRR, LCI, UCI, p.value) 


### CALCULATE PREDICTED VALUES FOR GRAPHING FROM MODELS
predicted_values <- predict(fit_pc_l, newdata = df_sum_exp, type = "link", se.fit = TRUE)

# Extract predicted values and standard errors
predicted_means <- predicted_values$fit
standard_errors <- predicted_values$se.fit

# Calculate confidence intervals
z_value <- qnorm(0.975)  # 1.96 for 95% confidence interval
lower_bound <- exp(predicted_means - z_value * standard_errors)
upper_bound <- exp(predicted_means + z_value * standard_errors)

df_sum_exp <- df_sum_exp %>%
  mutate(predicted_mean3 = (exp(predicted_means) / births_total) * 10000,
         predicted_lci3 = (lower_bound / births_total) * 10000,
         predicted_uci3 = (upper_bound / births_total) * 10000,
         prevs3 = (n_sb / births_total) * 10000)

# ---- Section 6: Time trends in all NTDs without a genetic condition ----

# First use Poisson model for association of total birth prevalence and year 
# this model assumes a linear trend over time
fit_p_l<-vglm(n_ntgng ~ year, poissonff(bred=TRUE), offset=log(births_total), 
              data=df_sum)
summary(fit_p_l)
# Produce AIC to compare Poisson models
AIC(fit_p_l) 
# check dispersion parameter for this model
dp_p_l <-  deviance(fit_p_l)/df.residual(fit_p_l, type = "vlm")
dp_p_l

## Then use alternative Poisson model for association of total birth prevalence and year 
# this model assumes a non-linear trend over time using restricted cubic splines
# define knots for this model based on trends observed on plot
knots <- c(2000, 2007, 2014, 2021)
fit_p_nl<-vglm(n_ntgng ~ rcs(year, knots), poissonff(bred=TRUE), offset=log(births_total), 
               data=df_sum)
summary(fit_p_nl)
# Produce AIC to compare Poisson models
AIC(fit_p_nl) 
# check dispersion parameter for this model
dp_p_nl <-  deviance(fit_p_nl)/df.residual(fit_p_nl, type = "vlm")
dp_p_nl

## Choosing the first model (fit_p_l) based on AIC criteria and after checking dispersion parameter

# calculate prevalence rate ratio (PRR) by exponentiating the d$time_period coefficient estimate
# and then create a dataframe to read exponentiated values from model
# Extract coefficients
coefficients <- coef(fit_p_l)
# Extract standard errors
std_errors <- sqrt(diag(vcov(fit_p_l)))
# Calculate z-values
z_values <- coefficients / std_errors
# Calculate p-values
p_values <- 2 * (1 - pnorm(abs(z_values)))
# Create a tidy data frame
tidy_fit <- data.frame(
  term = names(coefficients),
  estimate = coefficients,
  std.error = std_errors,
  z.value = z_values,
  p.value = p_values
)
#exponentiate estimates and calculate lower and upper confidence intervals
tidy_fit_nongenetic <- tidy_fit %>%
  mutate(PRR = round(exp(estimate), 3),
         LCI = round(exp(estimate - ((qnorm(0.975))*std.error)),3),
         UCI = round(exp(estimate + ((qnorm(0.975))*std.error)),3)) %>%
  subset(term=="year") %>%
  mutate(outcome = "All NTD with GC") %>%
  dplyr::select(term, outcome, PRR, LCI, UCI, p.value) 

## For selected Poisson model (fit_p_l in this case)
# Predict values
predicted_values <- predict(fit_p_l, newdata = df_sum, type = "link", se.fit = TRUE)

# Extract predicted values and standard errors
predicted_means <- predicted_values$fit
standard_errors <- predicted_values$se.fit

# Calculate confidence intervals
z_value <- qnorm(0.975)  # 1.96 for 95% confidence interval
lower_bound <- exp(predicted_means - z_value * standard_errors)
upper_bound <- exp(predicted_means + z_value * standard_errors)

df_sum_exp <- df_sum_exp %>%
  mutate(predicted_mean4 = (exp(predicted_means) / births_total) * 10000,
         predicted_lci4 = (lower_bound / births_total) * 10000,
         predicted_uci4 = (upper_bound / births_total) * 10000,
         prevs4 = (n_ntgng / births_total) * 10000)

# ---- Section 7: Graph the results ----

#Prepare dataset for graphing 
df_sum_graph <- df_sum_exp %>%
  dplyr::select(year, 
         predicted_mean1, predicted_lci1, predicted_uci1, prevs1,
         predicted_mean2, predicted_lci2, predicted_uci2, prevs2,
         predicted_mean3, predicted_lci3, predicted_uci3, prevs3,
         predicted_mean4, predicted_lci4, predicted_uci4, prevs4) %>%
  pivot_longer(
    cols=c(-year),
    names_to = c(".value", "ntd_group"),
    names_pattern="(.*)(\\d+)") 

#### # Graph the observed versus the modelled estimates of total birth prevalence of NTDs, anencephaly and spina bifida ####
df_sum_graph1 <- df_sum_graph %>%
  subset(ntd_group<4) 

facet_labels <- c("1"="All NTDs", "2"="Anencephaly", "3"="Spina Bifida")

p <- ggplot(df_sum_graph1, aes(x = year)) +
  geom_point(aes(y = prevs, shape="Observed total birth prevalence"), color = "black", size = 2) +
  geom_line(aes(y = predicted_mean, color = "Modelled prevalence estimates")) +
  geom_ribbon(aes(ymin = predicted_lci, ymax = predicted_uci, fill="95% confidence interval"), alpha = 0.2) +
  facet_wrap(~ntd_group, labeller=as_labeller(facet_labels), nrow=3) +
  labs(x = "Year",
       y = "Total birth prevalence per 10,000 total births") + 
  scale_y_continuous(limits = c(0, 14), breaks = seq(0, max(14), by = 1), expand = c(0, 0)) +  
  scale_shape_manual(name=NULL, values=c("Observed total birth prevalence"=1), breaks=c("Observed total birth prevalence")) +
  scale_color_manual(name=NULL, values=c("Modelled prevalence estimates"="blue"), breaks=c("Modelled prevalence estimates")) +
  scale_fill_manual(name=NULL, values=c("95% confidence interval"="blue"), breaks=c("95% confidence interval")) +
  guides(
    shape=guide_legend(order=1),
    color=guide_legend(order=2),
    fill=guide_legend(order=3)) +
  theme_minimal() +
  theme(
    legend.spacing.y = unit(0.05, "cm"),
    panel.border = element_rect(color="black", fill=NA, linewidth = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(), 
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.line.x = element_line(color = "grey"),  
    axis.line.y = element_line(color = "grey"))

ggsave("/conf/Congenital-Descriptive-Epi/Trends in NTDs/Data/Figure1.tiff", plot=p, width=7.5, height=9, units="in", dpi=300)

#### # Graph the observed versus the modelled estimates of total birth prevalence of NTDs without an associated genetic condition ####

df_sum_graph2 <- df_sum_graph %>%
  subset(ntd_group==4) 

p2 <- ggplot(df_sum_graph2, aes(x = year)) +
  geom_point(aes(y = prevs, shape="Observed total birth prevalence"), color = "black", size = 2) +
  geom_line(aes(y = predicted_mean, color = "Modelled prevalence estimates")) +
  geom_ribbon(aes(ymin = predicted_lci, ymax = predicted_uci, fill="95% confidence interval"), alpha = 0.2) +
  labs(x = "Year",
       y = "Total birth prevalence per 10,000 total births") + 
  scale_y_continuous(limits = c(0, 14), breaks = seq(0, max(14), by = 1), expand = c(0, 0)) +  
  scale_shape_manual(name=NULL, values=c("Observed total birth prevalence"=1), breaks=c("Observed total birth prevalence")) +
  scale_color_manual(name=NULL, values=c("Modelled prevalence estimates"="blue"), breaks=c("Modelled prevalence estimates")) +
  scale_fill_manual(name=NULL, values=c("95% confidence interval"="blue"), breaks=c("95% confidence interval")) +
  guides(
    shape=guide_legend(order=1),
    color=guide_legend(order=2),
    fill=guide_legend(order=3)) +
  theme_minimal() +
  theme(
    legend.spacing.y = unit(0.05, "cm"),
    panel.border = element_rect(color="black", fill=NA, linewidth = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(), 
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.line.x = element_line(color = "grey"),  
    axis.line.y = element_line(color = "grey"))

ggsave("/conf/Congenital-Descriptive-Epi/Trends in NTDs/Data/FigureS1.tiff", plot=p2, width=7.5, height=4, units="in", dpi=300)

# ---- Section 8: Format and save model results for final tables ----

#Append together model results for different outcomes
time_trend_results <- bind_rows(tidy_fit_all, tidy_fit_anen, tidy_fit_anen_2, tidy_fit_sb, tidy_fit_nongenetic)
write.csv(time_trend_results, "/conf/Congenital-Descriptive-Epi/Trends in NTDs/Data/time_trend_model_results.csv", row.names = TRUE)

# ---- End of file ----
