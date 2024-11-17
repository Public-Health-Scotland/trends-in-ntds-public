###### TRENDS IN NEURAL TUBE DEFECTS (NTD) R CODE SCRIPT

# This script prepares the descriptive tables

# ---- Section 1: Housekeeping ----

# clear the global environment
rm(list = ls())

# load required packages
require(tidyverse)
require(ggplot2)
require(lubridate)
require(broom)

# LOADING PACKAGES

# Read in prepared data require for descriptive analyses 
df <- read.csv("/conf/Congenital-Descriptive-Epi/Trends in NTDs/Data/analysis_ready.csv")
df2 <- read.csv("/conf/Congenital-Descriptive-Epi/Trends in NTDs/Data/sub_nogenetic_analysis_ready.csv")

# ---- Section 2: Total and subtype descriptive analysis ----

# ##### Table 1; Number of babies with, and birth prevalence of, NTDs by NTD subtype and pregnancy outcome ####

#Calculate total number of NTDs
df_sum <- df %>%
  summarise(
    count_all = sum(count_all),
    count_lb = sum(count_lb),
    NTD_all = sum(NTD_all),
    anen = sum(anen),
    ence = sum(ence),
    SB = sum(SB),
    NTD_all_lb = sum(NTD_all_lb),
    anen_lb = sum(anen_lb),
    ence_lb = sum(ence_lb),
    SB_lb = sum(SB_lb),
    NTD_all_sb = sum(NTD_all_sb),
    anen_sb = sum(anen_sb),
    ence_sb = sum(ence_sb),
    SB_sb = sum(SB_sb),
    NTD_all_top = sum(NTD_all_top),
    anen_top = sum(anen_top),
    ence_top = sum(ence_top),
    SB_top = sum(SB_top))

#Make births data long
df_birth <- df_sum %>%
  select(count_all, count_lb) %>%
  pivot_longer(cols = starts_with("count_"),
               names_to = "group",
               values_to = "births") %>%
  mutate(row = row_number())

#All NTDs
df_ntd <- df_sum %>%
  select(NTD_all, NTD_all_lb, NTD_all_sb, NTD_all_top) %>%
  pivot_longer(cols = starts_with("NTD_all"),
                names_to = "group",
                values_to = "count") %>%
  mutate(row = row_number())  

df_ntd$perc <- NA
df_ntd$perc[2:4] <- round((df_ntd$count[2:4]/df_ntd$count[1]*100))

#Anencephaly
df_anen <- df_sum %>%
  select(anen, anen_lb, anen_sb, anen_top) %>%
  pivot_longer(cols = starts_with("anen"),
               names_to = "group",
               values_to = "count") %>%
  mutate(row = row_number())  

df_anen$perc <- NA
df_anen$perc[2:4] <- round((df_anen$count[2:4]/df_anen$count[1]*100))

#Encephalocele
df_ence <- df_sum %>%
  select(ence, ence_lb, ence_sb, ence_top) %>%
  pivot_longer(cols = starts_with("ence"),
               names_to = "group",
               values_to = "count") %>%
  mutate(row = row_number())  

df_ence$perc <- NA
df_ence$perc[2:4] <- round((df_ence$count[2:4]/df_ence$count[1]*100))

#Spina bifida
df_sb <- df_sum %>%
  select(SB, SB_lb, SB_sb, SB_top) %>%
  pivot_longer(cols = starts_with("SB"),
               names_to = "group",
               values_to = "count") %>%
  mutate(row = row_number())  

df_sb$perc <- NA
df_sb$perc[2:4] <- round((df_sb$count[2:4]/df_sb$count[1]*100))

#Append together overall NTD information and NTD subtypes
table1 <- bind_rows(df_ntd, df_anen, df_ence, df_sb)

#Merge on information on total number of births/live births
table1 <- left_join(table1, df_birth, by=c("row"))

#Add on total and live birth prevalence and 95% CIs
table1 <- table1 %>%
  mutate(prev=round((count/births)*10000, 1)) %>%
  mutate(CI=map2(count, births, ~if(!is.na(.x) & !is.na(.y)) {
    binom.test(.x, .y)$conf.int*10000
  } else {
    c(NA_real_, NA_real_)
  }),
lower_CI = round(map_dbl(CI, 1), 1),
upper_CI = round(map_dbl(CI, 2), 1)) %>%
  mutate(number = if_else(is.na(perc), as.character(count), paste0(count, " (", perc, "%)")),
         prev_CI= if_else(is.na(prev), as.character("-"), paste0(prev, " (", lower_CI, ", ", upper_CI, ")"))) %>%
  select(group.x, number, births, prev_CI) 

# #### Table S1; Number of babies with, and birth prevalence of, NTDs without an associated genetic condition by NTD subtype and pregnancy outcome ####

#Calculate total number of NTDs
df_sum <- df2 %>%
  summarise(
    NTD_all = sum(NTD_all),
    anen = sum(anen),
    ence = sum(ence),
    SB = sum(SB),
    NTD_all_lb = sum(NTD_all_lb),
    anen_lb = sum(anen_lb),
    ence_lb = sum(ence_lb),
    SB_lb = sum(SB_lb),
    NTD_all_sb = sum(NTD_all_sb),
    anen_sb = sum(anen_sb),
    ence_sb = sum(ence_sb),
    SB_sb = sum(SB_sb),
    NTD_all_top = sum(NTD_all_top),
    anen_top = sum(anen_top),
    ence_top = sum(ence_top),
    SB_top = sum(SB_top))

#All NTDs without an associated congenital condition
df_ntd <- df_sum %>%
  select(NTD_all, NTD_all_lb, NTD_all_sb, NTD_all_top) %>%
  pivot_longer(cols = starts_with("NTD_all"),
               names_to = "group",
               values_to = "count") %>%
  mutate(row = row_number())  

df_ntd$perc <- NA
df_ntd$perc[2:4] <- round((df_ntd$count[2:4]/df_ntd$count[1]*100))

#Anencephaly without an associated congenital condition
df_anen <- df_sum %>%
  select(anen, anen_lb, anen_sb, anen_top) %>%
  pivot_longer(cols = starts_with("anen"),
               names_to = "group",
               values_to = "count") %>%
  mutate(row = row_number())  

df_anen$perc <- NA
df_anen$perc[2:4] <- round((df_anen$count[2:4]/df_anen$count[1]*100))

#Encephalocele without an associated congenital condition
df_ence <- df_sum %>%
  select(ence, ence_lb, ence_sb, ence_top) %>%
  pivot_longer(cols = starts_with("ence"),
               names_to = "group",
               values_to = "count") %>%
  mutate(row = row_number())  

df_ence$perc <- NA
df_ence$perc[2:4] <- round((df_ence$count[2:4]/df_ence$count[1]*100))

#Spina bifida without an associated congenital condition
df_sb <- df_sum %>%
  select(SB, SB_lb, SB_sb, SB_top) %>%
  pivot_longer(cols = starts_with("SB"),
               names_to = "group",
               values_to = "count") %>%
  mutate(row = row_number())  

df_sb$perc <- NA
df_sb$perc[2:4] <- round((df_sb$count[2:4]/df_sb$count[1]*100))

#Append together overall NTD information and NTD subtypes
tableS1 <- bind_rows(df_ntd, df_anen, df_ence, df_sb)

#Merge on information on total number of births/live births
tableS1 <- left_join(tableS1, df_birth, by=c("row"))

#Add on total and live birth prevalence and 95% CIs
tableS1 <- tableS1 %>%
  mutate(prev=round((count/births)*10000, 1)) %>%
  mutate(CI=map2(count, births, ~if(!is.na(.x) & !is.na(.y)) {
    binom.test(.x, .y)$conf.int*10000
  } else {
    c(NA_real_, NA_real_)
  }),
  lower_CI = round(map_dbl(CI, 1), 1),
  upper_CI = round(map_dbl(CI, 2), 1)) %>%
  mutate(number = if_else(is.na(perc), as.character(count), paste0(count, " (", perc, "%)")),
         prev_CI= if_else(is.na(prev), as.character("-"), paste0(prev, " (", lower_CI, ", ", upper_CI, ")"))) %>%
  select(group.x, number, births, prev_CI) 

# ---- Section 3: Sociodemographic descriptive analysis ----  

# ##### Table 2; Number of babies with, and birth prevalence of, NTDs by sociodemographic characteristics ####

# Maternal age
# Calculate total number of births and NTDs
df_age <- df %>%
  group_by(maternal_age_group) %>%
  summarise(
    births = sum(count_all),
    count1 = sum(NTD_all),
    count2 = sum(anen),
    count3 = sum(SB),
    count4 = sum(NTD_no_gen)) %>%
  mutate(births_perc=round(births/sum(births)*100),
         perc1=round(count1/sum(count1)*100),
         perc2=round(count2/sum(count2)*100),
         perc3=round(count3/sum(count3)*100),
         perc4=round(count4/sum(count4)*100)) %>%
  pivot_longer(
    cols = -c(births, births_perc, maternal_age_group),
    names_to = c(".value", "ntd_group"),
    names_pattern="(.*)(\\d+)",
  ) %>%
  mutate(prev=round((count/births)*10000, 1)) %>%
  mutate(CI=map2(count, births, ~if(!is.na(.x) & !is.na(.y)) {
    binom.test(.x, .y)$conf.int*10000
  } else {
    c(NA_real_, NA_real_)
  }),
  lower_CI = round(map_dbl(CI, 1), 1),
  upper_CI = round(map_dbl(CI, 2), 1)) %>%
  mutate(number_perc_ntd = paste0(count, " (", perc, "%)"),
         number_perc_birth = paste0(births, " (", births_perc, "%)"), 
         prev_CI= if_else(maternal_age_group=="Unknown", as.character("-"), paste0(prev, " (", lower_CI, ", ", upper_CI, ")")))

# SIMD
# Calculate total number of births and NTDs
df_simd <- df %>%
  group_by(simd) %>%
  summarise(
    births = sum(count_all),
    count1 = sum(NTD_all),
    count2 = sum(anen),
    count3 = sum(SB),
    count4 = sum(NTD_no_gen)) %>%
  mutate(births_perc=round(births/sum(births)*100),
         perc1=round(count1/sum(count1)*100),
         perc2=round(count2/sum(count2)*100),
         perc3=round(count3/sum(count3)*100),
         perc4=round(count4/sum(count4)*100)) %>%
  pivot_longer(
    cols = -c(births, births_perc, simd),
    names_to = c(".value", "ntd_group"),
    names_pattern="(.*)(\\d+)",
  ) %>%
  mutate(prev=round((count/births)*10000, 1)) %>%
  mutate(CI=map2(count, births, ~if(!is.na(.x) & !is.na(.y)) {
    binom.test(.x, .y)$conf.int*10000
  } else {
    c(NA_real_, NA_real_)
  }),
  lower_CI = round(map_dbl(CI, 1), 1),
  upper_CI = round(map_dbl(CI, 2), 1)) %>%
  mutate(number_perc_ntd = paste0(count, " (", perc, "%)"),
         number_perc_birth = paste0(births, " (", births_perc, "%)"), 
         prev_CI= if_else(simd=="Unknown", as.character("-"), paste0(prev, " (", lower_CI, ", ", upper_CI, ")")))

# Infant sex
# Calculate total number of live births and NTDs live births
df_sex <- df %>%
  group_by(sex) %>%
  summarise(
    births = sum(count_lb),
    count1 = sum(NTD_all_lb),
    count2 = sum(anen_lb),
    count3 = sum(SB_lb),
    count4 = sum(NTD_no_gen_lb)) %>%
  mutate(births_perc=round(births/sum(births)*100),
         perc1=round(count1/sum(count1)*100),
         perc2=round(count2/sum(count2)*100),
         perc3=round(count3/sum(count3)*100),
         perc4=round(count4/sum(count4)*100)) %>%
  pivot_longer(
    cols = -c(births, births_perc, sex),
    names_to = c(".value", "ntd_group"),
    names_pattern="(.*)(\\d+)",
  ) %>%
  mutate(prev=round((count/births)*10000, 1)) %>%
  mutate(CI=map2(count, births, ~if(!is.na(.x) & !is.na(.y)) {
    binom.test(.x, .y)$conf.int*10000
  } else {
    c(NA_real_, NA_real_)
  }),
  lower_CI = round(map_dbl(CI, 1), 1),
  upper_CI = round(map_dbl(CI, 2), 1)) %>%
  mutate(number_perc_ntd = paste0(count, " (", perc, "%)"),
         number_perc_birth = paste0(births, " (", births_perc, "%)"), 
         prev_CI= if_else(sex=="Unknown", as.character("-"), paste0(prev, " (", lower_CI, ", ", upper_CI, ")")))

#Extract data for table 2 

#maternal age
table2_mat_age <- df_age %>%
  subset(ntd_group==1) %>%
  select(maternal_age_group, number_perc_ntd, number_perc_birth, prev_CI) %>%
  rename(var=maternal_age_group)
new_row <- data.frame(
  var="Maternal age (years)",
  number_perc_ntd = " ",
  number_perc_birth = " ",
  prev_CI=" ")
table2_mat_age <- bind_rows(new_row, table2_mat_age)

#simd
table2_simd <- df_simd %>%
  subset(ntd_group==1) %>%
  select(simd, number_perc_ntd, number_perc_birth, prev_CI) %>%
  rename(var=simd)
new_row <- data.frame(
  var="Maternal deprivation status",
  number_perc_ntd = " ",
  number_perc_birth = " ",
  prev_CI=" ")
table2_simd <- bind_rows(new_row, table2_simd)

#infant sex  
table2_sex <- df_sex %>%
  subset(ntd_group==1) %>%
  select(sex, number_perc_ntd, number_perc_birth, prev_CI) %>%
  rename(var=sex)
new_row <- data.frame(
  var="Infant sex",
  number_perc_ntd = " ",
  number_perc_birth = " ",
  prev_CI=" ")
table2_sex <- bind_rows(new_row, table2_sex)

#Bring together different variables to make final Table 2
table2 <- bind_rows(table2_mat_age, table2_simd, table2_sex)

# ##### Table S2; Number of babies with, and birth prevalence of, NTDs by subtype and sociodemographic characteristics ####

#Extract data for table S2 based on dataframes that have already been prepared

#maternal age
tableS2_mat_age <- df_age %>%
  subset(ntd_group>1) %>%
  select(ntd_group, maternal_age_group, births, number_perc_ntd, prev_CI) %>%
  rename(var=maternal_age_group) %>%
  pivot_wider(
    names_from=ntd_group,
    values_from= c(number_perc_ntd, prev_CI),
    id_cols=c(var, births)
  ) %>%
  select(var, births, number_perc_ntd_2, prev_CI_2, number_perc_ntd_3, prev_CI_3, number_perc_ntd_4, prev_CI_4) %>%
  rename(anencephaly_n=number_perc_ntd_2,
         anencephaly_prev=prev_CI_2,
         spina_bifida_n=number_perc_ntd_3,
         spina_bifida_prev=prev_CI_3,
         NTD_nogen_n=number_perc_ntd_4,
         NTD_nogen_prev=prev_CI_4,
         )
new_row <- data.frame(
  var="Maternal age (years)")
tableS2_mat_age <- bind_rows(new_row, tableS2_mat_age)

#SIMD
tableS2_simd <- df_simd %>%
  subset(ntd_group>1) %>%
  select(ntd_group, simd, births, number_perc_ntd, prev_CI) %>%
  rename(var=simd) %>%
  pivot_wider(
    names_from=ntd_group,
    values_from= c(number_perc_ntd, prev_CI),
    id_cols=c(var, births)
  ) %>%
  select(var, births, number_perc_ntd_2, prev_CI_2, number_perc_ntd_3, prev_CI_3, number_perc_ntd_4, prev_CI_4) %>%
  rename(anencephaly_n=number_perc_ntd_2,
         anencephaly_prev=prev_CI_2,
         spina_bifida_n=number_perc_ntd_3,
         spina_bifida_prev=prev_CI_3,
         NTD_nogen_n=number_perc_ntd_4,
         NTD_nogen_prev=prev_CI_4,
  )
new_row <- data.frame(
  var="Maternal deprivation status")
tableS2_simd <- bind_rows(new_row, tableS2_simd)

#Infant sex
tableS2_sex <- df_sex %>%
  subset(ntd_group>1) %>%
  select(ntd_group, sex, births, number_perc_ntd, prev_CI) %>%
  rename(var=sex) %>%
  pivot_wider(
    names_from=ntd_group,
    values_from= c(number_perc_ntd, prev_CI),
    id_cols=c(var, births)
  ) %>%
  select(var, births, number_perc_ntd_2, prev_CI_2, number_perc_ntd_3, prev_CI_3, number_perc_ntd_4, prev_CI_4) %>%
  rename(anencephaly_n=number_perc_ntd_2,
         anencephaly_prev=prev_CI_2,
         spina_bifida_n=number_perc_ntd_3,
         spina_bifida_prev=prev_CI_3,
         NTD_nogen_n=number_perc_ntd_4,
         NTD_nogen_prev=prev_CI_4,
  )
new_row <- data.frame(
  var="Sex")
tableS2_sex <- bind_rows(new_row, tableS2_sex)

#Bring together different variables to make final Table S2
tableS2 <- bind_rows(tableS2_mat_age, tableS2_simd, tableS2_sex)

# ---- Section 4: Save tables ----

write.csv(table1, "/conf/Congenital-Descriptive-Epi/Trends in NTDs/Data/table1.csv", row.names = TRUE)
write.csv(tableS1, "/conf/Congenital-Descriptive-Epi/Trends in NTDs/Data/tableS1.csv", row.names = TRUE)
write.csv(table2, "/conf/Congenital-Descriptive-Epi/Trends in NTDs/Data/table2.csv", row.names = TRUE)
write.csv(tableS2, "/conf/Congenital-Descriptive-Epi/Trends in NTDs/Data/tableS2.csv", row.names = TRUE)

# ---- End of file ---
