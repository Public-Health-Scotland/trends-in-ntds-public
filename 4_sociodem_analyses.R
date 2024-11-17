###### TRENDS IN NEURAL TUBE DEFECTS (NTD) R CODE SCRIPT

# This script conducts the analyses looking at the association between sociodemographic characteristics and NTDs

# ---- Section 1: Housekeeping ----

# clear the global environment
rm(list = ls())

# load required packages
require(tidyverse)
require(ggplot2)
require(lubridate)
require(broom)
require(VGAM) # for VGLM model
require(rms) # for restricted cubic splines 'rcs'
require(glmmTMB) # for Poisson COM model

# Read in prepared data require for descriptive analyses 
df <- read.csv("/conf/Congenital-Descriptive-Epi/Trends in NTDs/Data/analysis_ready.csv")

# ---- Section 2: Set up functions for formatting data from models ----

#Function model_tidy takes output from poisson models and creates
#a dataframe with the PRR, 95% CI and p-value

model_tidy <- function(modelx, outcomex) {
  
  # Create a dataframe to read exponentiated values from model
  coefficients <- coef(modelx)
  std_errors <- sqrt(diag(vcov(modelx)))
  z_values <- coefficients / std_errors
  p_values <- 2 * (1 - pnorm(abs(z_values)))
  
  # Create a tidy data frame
  tidy_fit <- data.frame(
    term = names(coefficients),
    estimate = coefficients,
    std.error = std_errors,
    z.value = z_values,
    p.value = p_values
  )
  # exponentiate estimates and calculate lower and upper confidence intervals
  dataframe <- tidy_fit %>%
    mutate(PRR = round(exp(estimate), 3),
           LCI = round(exp(estimate - ((qnorm(0.975))*std.error)),3),
           UCI = round(exp(estimate + ((qnorm(0.975))*std.error)),3),
           p.value = round(p.value,3)) %>%
    subset(term!="(Intercept)") %>%
    mutate(outcome = outcomex) %>%
    select(term, outcome, PRR, LCI, UCI, p.value) 
}

#Function model_tidy_com takes output from COM-poisson models and creates
#a dataframe with the PRR, 95% CI and p-value

model_tidy_com <- function(modelx, outcomex) {
  
  coefficients <- as.numeric(fixef(modelx)$cond)
  std_errors <- sqrt(diag(vcov(modelx)$cond))
  z_values <- coefficients / std_errors
  p_values <- 2 * (1 - pnorm(abs(z_values)))
  exp_coefficients <- exp(coefficients)
  lower_ci <- exp(coefficients - 1.96 * std_errors)
  upper_ci <- exp(coefficients + 1.96 * std_errors)
  
  # Create a tidy data frame
  tidy_fit <- data.frame(
    term = names(fixef(modelx)$cond),
    estimate = coefficients,
    PRR = exp_coefficients,
    std.error = std_errors,
    lower_ci = lower_ci,
    upper_ci = upper_ci,
    z.value = z_values,
    p.value = p_values
  )
  
  #exponentiate estimates and calculate lower and upper confidence intervals
  dataframe <- tidy_fit %>%
    mutate(PRR = round(exp(estimate), 3),
           LCI = round(exp(estimate - ((qnorm(0.975))*std.error)),3),
           UCI = round(exp(estimate + ((qnorm(0.975))*std.error)),3), 
           p.value = round(p.value,3)) %>%
    subset(term!="(Intercept)") %>%
    mutate(outcome = outcomex) %>%
    select(term, outcome, PRR, LCI, UCI, p.value)
}

# ---- Section 3: Prepare data frames for models ----

# Set up age group and specify 25-29 as baseline
df$maternal_age_group2 <- factor(df$maternal_age_group,
                                 levels = c(3, 1, 2, 4, 5, 6),
                                 labels = c("25-29", "<20", "20-24", "30-34", "35-39", "40+"))

# Set up SIMD and specify group 5 as baseline
df$simd2 <- factor(df$simd,
                   levels = c(5, 1, 2, 3, 4),
                   labels = c("5", "1", "2", "3", "4"))

# Prepare data for looking at association between maternal age and outcomes
df_age <- df %>%
  group_by(maternal_age_group2) %>%
  summarise(
    births_total = sum(count_all),
    n_ntd = sum(NTD_all),
    n_anen = sum(anen),
    n_sb = sum(SB),
    n_ntdng = sum(NTD_no_gen)) %>%
  subset(!is.na(maternal_age_group2)) %>%
  ungroup() 
sum(df_age$births_total)
sum(df_age$n_ntd)

df_age_year <- df %>%
  group_by(maternal_age_group2, year) %>%
  summarise(
    births_total = sum(count_all),
    n_ntd = sum(NTD_all),
    n_anen = sum(anen),
    n_sb = sum(SB),
    n_ntdng = sum(NTD_no_gen)) %>%
  subset(!is.na(maternal_age_group2)) %>%
  ungroup() %>%
  mutate(year2=year-2000)
sum(df_age_year$births_total)
sum(df_age_year$n_ntd)

# Prepare data for looking at association between SIMD and outcomes
df_simd <- df %>%
  group_by(simd2) %>%
  summarise(
    births_total = sum(count_all),
    n_ntd = sum(NTD_all),
    n_anen = sum(anen),
    n_sb = sum(SB),
    n_ntdng = sum(NTD_no_gen)) %>%
  subset(!is.na(simd2)) %>%
  ungroup() 
sum(df_simd$births_total)
sum(df_simd$n_ntd)

df_simd_year <- df %>%
  group_by(simd2, year) %>%
  summarise(
    births_total = sum(count_all),
    n_ntd = sum(NTD_all),
    n_anen = sum(anen),
    n_sb = sum(SB),
    n_ntdng = sum(NTD_no_gen)) %>%
  subset(!is.na(simd2)) %>%
  ungroup() %>%
  mutate(year2=year-2000)
sum(df_simd_year$births_total)
sum(df_simd_year$n_ntd)

# Prepare data for looking at association between maternal age, SIMD and year and outcomes
df_age_simd_year <- df %>%
  group_by(maternal_age_group2, simd2, year) %>%
  summarise(
    births_total = sum(count_all),
    n_ntd = sum(NTD_all),
    n_anen = sum(anen),
    n_sb = sum(SB),
    n_ntdng = sum(NTD_no_gen)) %>%
  subset(!is.na(simd2)) %>%
  subset(!is.na(maternal_age_group2)) %>%
  ungroup() %>%
  mutate(year2=year-2000)
sum(df_age_simd_year$births_total)
sum(df_age_simd_year$n_ntd)

# ---- Section 4: Models for associations with all NTDs ----

#### Maternal age 

# Poisson model for association of maternal age with total birth prevalence
fit_ntd_age<-vglm(n_ntd ~ maternal_age_group2, poissonff(bred=TRUE), offset=log(births_total), 
              data=df_age)
summary(fit_ntd_age)

# create tidy output for merging to create table
tidy_fit_ntd_age <- model_tidy(fit_ntd_age, "All NTD")

# Poisson model for association of maternal age with total birth prevalence adjusted for year
fit_ntd_age_year<-vglm(n_ntd ~ maternal_age_group2 + year, poissonff(bred=TRUE), offset=log(births_total), 
              data=df_age_year)
summary(fit_ntd_age_year)

# dispersion parameter
dp_ntd_age_year <-  deviance(fit_ntd_age_year)/df.residual(fit_ntd_age_year, type = "vlm")
dp_ntd_age_year

# create tidy output for merging to create table
tidy_fit_ntd_age_year <- model_tidy(fit_ntd_age_year, "All NTD")
tidy_fit_ntd_age_year <- tidy_fit_ntd_age_year %>%
  rename(PRR2=PRR,
         LCI2=LCI,
         UCI2=UCI,
         p.value2=p.value)

#### SIMD 

# Poisson model for association between SIMD and total birth prevalence
fit_ntd_simd<-vglm(n_ntd ~ simd2, poissonff(bred=TRUE), offset=log(births_total), 
                   data=df_simd)
summary(fit_ntd_simd)

# create tidy output for merging to create table
tidy_fit_ntd_simd <- model_tidy(fit_ntd_simd, "All NTD")

# Poisson model for association of SIMD and year with total birth prevalence
fit_ntd_simd_year<-vglm(n_ntd ~ simd2 + year, poissonff(bred=TRUE), offset=log(births_total), 
                        data=df_simd_year)
summary(fit_ntd_simd_year)

# dispersion parameter
dp_ntd_simd_year <-  deviance(fit_ntd_simd_year)/df.residual(fit_ntd_simd_year, type = "vlm")
dp_ntd_simd_year

# create tidy output for merging to create table
tidy_fit_ntd_simd_year <- model_tidy(fit_ntd_simd_year, "All NTD")
tidy_fit_ntd_simd_year <- tidy_fit_ntd_simd_year %>%
  rename(PRR3=PRR,
         LCI3=LCI,
         UCI3=UCI,
         p.value3=p.value)

#### Age, SIMD and Year 

# Poisson model for association of total birth prevalence with simd, maternal age, and year
fit_ntd_age_simd_year<-vglm(n_ntd ~ maternal_age_group2 + simd2 + year, poissonff(bred=TRUE), offset=log(births_total), 
              data=df_age_simd_year)
summary(fit_ntd_age_simd_year)

# dispersion parameter
dp_p_l <-  deviance(fit_ntd_age_simd_year)/df.residual(fit_ntd_age_simd_year, type = "vlm")
dp_p_l

# create tidy output for merging to create table
tidy_fit_ntd_age_simd_year <- model_tidy(fit_ntd_age_simd_year, "All NTD")
tidy_fit_ntd_age_simd_year <- tidy_fit_ntd_age_simd_year %>%
  rename(PRR4=PRR,
         LCI4=LCI,
         UCI4=UCI,
         p.value4=p.value)

# create Table 3 by joining together final models results for all NTDs
table3 <- rbind(tidy_fit_ntd_age, tidy_fit_ntd_simd)
table3 <- table3 %>%
  full_join(tidy_fit_ntd_age_year, by = c("term", "outcome")) %>%
  full_join(tidy_fit_ntd_simd_year, by = c("term", "outcome")) %>%
  full_join(tidy_fit_ntd_age_simd_year, by = c("term", "outcome"))

# ---- Section 5: Models for associations with anencephaly ----

#### Maternal age 

# Poisson model for association of maternal age with anencephaly
fit_anen_age<-vglm(n_anen ~ maternal_age_group2, poissonff(bred=TRUE), offset=log(births_total), 
              data=df_age)
summary(fit_anen_age)

tidy_fit_anen_age <- model_tidy(fit_anen_age, "Anencephaly")

# COM-Poisson model for association of maternal age and year with anencephaly 
# Set knots for model assuming a non-linear trend over time
knots <- c(0, 13, 21)
# COM-Poisson model for association of total birth prevalence with maternal age adjusted for year
fit_anen_age_year <- glmmTMB(n_anen ~ maternal_age_group2 + rcs(year2, knots) + offset(log(births_total)), family=compois, data=df_age_year)
summary(fit_anen_age_year)

# create tidy output for merging to create table
tidy_fit_anen_age_year <- model_tidy_com(fit_anen_age_year, "Anencephaly")
tidy_fit_anen_age_year <- tidy_fit_anen_age_year %>%
  rename(PRR2=PRR,
         LCI2=LCI,
         UCI2=UCI,
         p.value2=p.value)

#### SIMD 

# Poisson regression model for association of SIMD and anencephaly
fit_anen_simd<-vglm(n_anen ~ simd2, poissonff(bred=TRUE), offset=log(births_total), 
                    data=df_simd)
summary(fit_anen_simd)

# create tidy output for merging to create table
tidy_fit_anen_simd <- model_tidy(fit_ntd_simd, "Anencephaly")

# COM-Poisson model for association between SIMD and year with anencephaly 
# This model assumes a non-linear trend over time based on knots 
# Set knots for model 
knots <- c(0, 13, 21)
fit_anen_simd_year <- glmmTMB(n_anen ~ simd2 + rcs(year2,knots) + offset(log(births_total)), family=compois, data=df_simd_year)
summary(fit_anen_simd_year)

# create tidy output for merging to create table
tidy_fit_anen_simd_year <- model_tidy_com(fit_anen_simd_year, "Anencephaly")
tidy_fit_anen_simd_year <- tidy_fit_anen_simd_year %>%
  rename(PRR3=PRR,
         LCI3=LCI,
         UCI3=UCI,
         p.value3=p.value)

#### Age, SIMD and Year 

# COM-Poisson COM model for association of total anencephaly birth prevalence with year, simd, and maternal age
# Set knots for model assuming non-linear trend over time
knots <- c(0, 13, 21)
fit_anen_age_simd_year <- glmmTMB(n_anen ~ maternal_age_group2 + simd2 + rcs(year2,knots) + offset(log(births_total)), family=compois, 
                     data=df_age_simd_year)
summary(fit_anen_age_simd_year)

# create tidy output for merging to create table
tidy_fit_anen_age_simd_year <- model_tidy_com(fit_anen_age_simd_year, "Anencephaly")
tidy_fit_anen_age_simd_year <- tidy_fit_anen_age_simd_year %>%
  rename(PRR4=PRR,
         LCI4=LCI,
         UCI4=UCI,
         p.value4=p.value)

# create Supplementary Table 3 by joining together final models results for anencephaly
# drop terms for year are these are not interpretable 
tableS3 <- rbind(tidy_fit_anen_age, tidy_fit_anen_simd)
tableS3 <- tableS3 %>%
  full_join(tidy_fit_anen_age_year, by = c("term", "outcome")) %>%
  full_join(tidy_fit_anen_simd_year, by = c("term", "outcome")) %>%
  full_join(tidy_fit_anen_age_simd_year, by = c("term", "outcome")) %>%
  filter(!grepl("knots", term, ignore.case=TRUE))

# ---- Section 6: Models for associations with spina bifida ----

#### Maternal age 

# Poisson model for association of maternal age and spina bifida
fit_sb_age<-vglm(n_sb ~ maternal_age_group2, poissonff(bred=TRUE), offset=log(births_total), 
              data=df_age)
summary(fit_sb_age)

# create tidy output for merging to create table
tidy_fit_sb_age <- model_tidy(fit_sb_age, "Spina Bifida")

# COM-Poisson for association of total birth prevalence with maternal age and year assuming a linear trend over time
fit_sb_age_year <- glmmTMB(n_sb ~ maternal_age_group2 + year2 + offset(log(births_total)), family=compois, data=df_age_year)
summary(fit_sb_age_year)

# create tidy output for merging to create table
tidy_fit_sb_age_year <- model_tidy_com(fit_sb_age_year, "Anencephaly")
tidy_fit_sb_age_year <- tidy_fit_sb_age_year %>%
  rename(PRR2=PRR,
         LCI2=LCI,
         UCI2=UCI,
         p.value2=p.value)

#### SIMD 

# Poisson model for association of SIMD with spina bifida
fit_sb_simd<-vglm(n_sb ~ simd2, poissonff(bred=TRUE), offset=log(births_total), 
                  data=df_simd)
summary(fit_sb_simd)

# create tidy output for merging to create table
tidy_fit_sb_simd <- model_tidy(fit_sb_simd, "Spina Bifida")

# COM-Poisson model for association between SIMD and year with spina bifida
fit_sb_simd_year <- glmmTMB(n_sb ~ simd2 + year2 + offset(log(births_total)), family=compois, data=df_simd_year)
summary(fit_sb_simd_year)

# create tidy output for merging to create table
tidy_fit_sb_simd_year <- model_tidy_com(fit_sb_simd_year, "Spina Bifida")
tidy_fit_sb_simd_year <- tidy_fit_sb_simd_year %>%
  rename(PRR3=PRR,
         LCI3=LCI,
         UCI3=UCI,
         p.value3=p.value)

#### Age, SIMD and Year 

# COM-Poisson model for association of total birth prevalence with year assuming a linear time trend
fit_sb_age_simd_year <- glmmTMB(n_sb ~ maternal_age_group2 + simd2 + year2 + offset(log(births_total)), family=compois, 
                    data=df_age_simd_year)
summary(fit_sb_age_simd_year)

# create tidy output for merging to create table
tidy_fit_sb_age_simd_year <- model_tidy_com(fit_sb_age_simd_year, "Spina Bifida")
tidy_fit_sb_age_simd_year <- tidy_fit_sb_age_simd_year %>%
  rename(PRR3=PRR,
         LCI3=LCI,
         UCI3=UCI,
         p.value3=p.value)

# create Supplementary Table 4 by joining together final models results for spina bifida
tableS4 <- rbind(tidy_fit_sb_age, tidy_fit_sb_simd)
tableS4 <- tableS4 %>%
  full_join(tidy_fit_sb_age_year, by = c("term", "outcome")) %>%
  full_join(tidy_fit_sb_simd_year, by = c("term", "outcome")) %>%
  full_join(tidy_fit_sb_age_simd_year, by = c("term", "outcome"))

# ---- Section 7: Models for associations with all NTDS without an associated genetic condition ----

#### Maternal age 

# Poisson model for association of maternal age with NTD without an associated genetic condition
fit_ntdng_age<-vglm(n_ntdng ~ maternal_age_group2, poissonff(bred=TRUE), offset=log(births_total), 
              data=df_age)
summary(fit_ntdng_age)

# create tidy output for merging to create table
tidy_fit_ntdng_age <- model_tidy(fit_ntdng_age, "NTD without GC")

# Poisson model for association of maternal age and year with NTD without an associated genetic condition
fit_ntdng_age_year<-vglm(n_ntdng ~ maternal_age_group2 + year, poissonff(bred=TRUE), offset=log(births_total), 
              data=df_age_year)
summary(fit_ntdng_age_year)

# dispersion parameter
dp_ntdng_age_year <-  deviance(fit_ntdng_age_year)/df.residual(fit_ntdng_age_year, type = "vlm")
dp_ntdng_age_year

# create tidy output for merging to create table
tidy_fit_ntdng_age_year <- model_tidy(fit_ntdng_age_year, "NTD without GC")
tidy_fit_ntdng_age_year <- tidy_fit_ntdng_age_year %>%
  rename(PRR2=PRR,
         LCI2=LCI,
         UCI2=UCI,
         p.value2=p.value)

#### SIMD 

# Poisson model for association between SIMD with NTDs without an associated genetic condition
fit_ntdng_simd<-vglm(n_ntdng ~ simd2, poissonff(bred=TRUE), offset=log(births_total), 
                     data=df_simd)
summary(fit_ntdng_simd)

# create tidy output for merging to create table
tidy_fit_ntdng_simd <- model_tidy(fit_ntdng_simd, "NTD without GC")

# Poisson model for association of SIMD and year with NTDs without an associated genetic condition
fit_ntdng_simd_year<-vglm(n_ntdng ~ simd2 + year, poissonff(bred=TRUE), offset=log(births_total), 
                          data=df_simd_year)
summary(fit_ntdng_simd_year)

# dispersion parameter
dp_ntdng_simd_year <-  deviance(fit_ntdng_simd_year)/df.residual(fit_ntdng_simd_year, type = "vlm")
dp_ntdng_simd_year

# create tidy output for merging to create table
tidy_fit_ntdng_simd_year <- model_tidy(fit_ntdng_simd_year, "NTD without GC")
tidy_fit_ntdng_simd_year <- tidy_fit_ntdng_simd_year %>%
  rename(PRR3=PRR,
         LCI3=LCI,
         UCI3=UCI,
         p.value3=p.value)

#### Age, SIMD and Year 

# Poisson model for association of simd, maternal age, and year with NTD within an associated genetic condition
fit_ntdng_age_simd_year<-vglm(n_ntdng ~ maternal_age_group2 + simd2 + year, poissonff(bred=TRUE), offset=log(births_total), 
              data=df_age_simd_year)
summary(fit_ntdng_age_simd_year)

# dispersion parameter
dp_ntdng_age_simd_year <-  deviance(fit_ntdng_age_simd_year)/df.residual(fit_ntdng_age_simd_year, type = "vlm")
dp_ntdng_age_simd_year

# create tidy output for merging to create table
tidy_fit_ntdng_age_simd_year <- model_tidy(fit_ntdng_age_simd_year, "NTD without GC")
tidy_fit_ntdng_age_simd_year <- tidy_fit_ntdng_age_simd_year %>%
  rename(PRR4=PRR,
         LCI4=LCI,
         UCI4=UCI,
         p.value4=p.value)

# create Supplementary Table 5 by joining together final models results for NTD within an associated genetic condition
tableS5 <- rbind(tidy_fit_ntdng_age, tidy_fit_ntdng_simd)
tableS5 <- tableS5 %>%
  full_join(tidy_fit_ntdng_age_year, by = c("term", "outcome")) %>%
  full_join(tidy_fit_ntdng_simd_year, by = c("term", "outcome")) %>%
  full_join(tidy_fit_ntdng_age_simd_year, by = c("term", "outcome"))

# ---- Section 8: Save table outputs ----

write.csv(table3, "/conf/Congenital-Descriptive-Epi/Trends in NTDs/Data/table3.csv", row.names = TRUE)
write.csv(tableS3, "/conf/Congenital-Descriptive-Epi/Trends in NTDs/Data/tableS3.csv", row.names = TRUE)
write.csv(tableS4, "/conf/Congenital-Descriptive-Epi/Trends in NTDs/Data/tableS4.csv", row.names = TRUE)
write.csv(tableS5, "/conf/Congenital-Descriptive-Epi/Trends in NTDs/Data/tableS5.csv", row.names = TRUE)

# ---- End of file ----
