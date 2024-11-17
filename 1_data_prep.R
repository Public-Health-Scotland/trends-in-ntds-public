###### TRENDS IN NEURAL TUBE DEFECTS (NTD) R CODE SCRIPT

# This script prepares the datasets for analysis

# ---- Section 1: Housekeeping ----

# load required packages
require(tidyverse)
require(ggplot2)
require(lubridate)
require(broom)

# Read in the National Records of Scotland (NRS) births and Scottish Linked Congenital Condition (SLICCD) datasets
NRSbirths <- readRDS(paste0("/conf/Congenital-Descriptive-Epi/Trends in NTDs/Data/nrs_births_data.rds"))
SLICCD <- readRDS(paste0("/conf/Congenital-Descriptive-Epi/Trends in NTDs/Data/sliccd_all_cc.rds"))

# ---- Section 2: Prepare denominator data (number of all births + live births) ----

#recategorise NA to be unknown so consistent across variables
NRSbirths$baby_sex <- ifelse(is.na(NRSbirths$baby_sex), "Unknown", NRSbirths$baby_sex) 
NRSbirths$maternal_age_group <- ifelse(is.na(NRSbirths$maternal_age_group), "Unknown", NRSbirths$maternal_age_group) 
NRSbirths$maternal_SIMD_end_preg[NRSbirths$maternal_SIMD_end_preg == "1 (most deprived)"] <- 1
NRSbirths$maternal_SIMD_end_preg[NRSbirths$maternal_SIMD_end_preg == "5 (least deprived)"] <- 5

#make the dataset long format
NRSbirth_long <- NRSbirths %>%
  pivot_longer(
    cols = starts_with("20"),
    names_to = "year",
    values_to = "count") %>%
  replace_na(list(count=0)) %>%
  mutate(year=as.numeric(year)) %>%
  rename(sex=baby_sex,
         simd=maternal_SIMD_end_preg)
sum(NRSbirth_long$count)

#create a dataset for live births
NRSbirth_live <- NRSbirth_long %>%
  filter(birth_type=="B") %>%
  rename(count_lb=count) %>%
  select(-birth_type)
sum(NRSbirth_live$count_lb)

#create a dataset for all births
NRSbirth_all <- NRSbirth_long %>%
  group_by(maternal_age_group, simd, sex, year) %>%
  summarise(count_all = sum(count)) 
sum(NRSbirth_all$count_all)

#merge together denominator datasets including count for all births and count for live births
Denom <- full_join(NRSbirth_all, NRSbirth_live, by = c("year", "sex", "simd", "maternal_age_group"))
sum(Denom$count_all)
sum(Denom$count_lb, na.rm=TRUE)

# ---- Section 3: Prepare numerator data (number of NTDs) ----

#Prepare the dataset by selecting relevant variables, and setting up our exposure
#variables (i.e., year, maternal age, maternal SIMD and baby sex) so they are formatted
#in a consistent way with the denominator data
All_NTD <- SLICCD %>%
filter(ALL_1_1_NEURAL_TUBE_DEFECTS == 1) %>%
select(ALL_1_1_NEURAL_TUBE_DEFECTS, 
       ALL_1_1_1_ANENCEPHALUS, 
       ALL_1_1_2_ENCEPHALOCELE, 
       ALL_1_1_3_SPINA_BIFIDA, 
       ALL_13_GENETIC_CONDITIONS, 
       time_period, 
       pregnancy_end_type, 
       mother_age,
       simd,
       sex) %>%
mutate(maternal_age_group=cut(
   mother_age,
   breaks = c(0, 20, 25, 30, 35, 40, Inf),
   labels = c("<20", "20-24", "25-29", "30-34", "35-39", "40+"),
   right = FALSE)) %>%
mutate(
  sex= case_when(sex==1 ~ "M",
                 sex==2 ~ "F",
                 sex==3 ~ "Unknown"),
  pregnancy_end_type = case_when(pregnancy_end_type=="Livebirth" ~ "lb",
                                 pregnancy_end_type=="Late Fetal Death" ~ "sb", 
                                 pregnancy_end_type=="Stillbirth" ~ "sb",
                                 pregnancy_end_type=="TOPFA" ~ "top")) %>%
rename(year=time_period)

All_NTD$simd <- ifelse(is.na(All_NTD$simd), "Unknown", All_NTD$simd) 
All_NTD$maternal_age_group <- ifelse(is.na(All_NTD$maternal_age_group), "Unknown", All_NTD$maternal_age_group) 

#Create aggregate count data which can be merged to denominator data

#Calculate total number of NTDs
All_NTD_num <- All_NTD %>%
  group_by(year, maternal_age_group, simd, sex) %>%
  summarise(
    NTD_all = sum(ALL_1_1_NEURAL_TUBE_DEFECTS, na.rm=TRUE),
    anen = sum(ALL_1_1_1_ANENCEPHALUS, na.rm=TRUE),
    ence = sum(ALL_1_1_2_ENCEPHALOCELE, na.rm=TRUE),
    SB = sum(ALL_1_1_3_SPINA_BIFIDA, na.rm=TRUE),
    NTD_no_gen = sum(ALL_13_GENETIC_CONDITIONS==0, na.rm=TRUE))

#Calculate NTD totals by pregnancy outcome 
All_NTD_num_pregend <- All_NTD %>%
  group_by(year, pregnancy_end_type, maternal_age_group, simd, sex) %>%
  summarise(
    NTD_all = sum(ALL_1_1_NEURAL_TUBE_DEFECTS, na.rm=TRUE),
    anen = sum(ALL_1_1_1_ANENCEPHALUS, na.rm=TRUE),
    ence = sum(ALL_1_1_2_ENCEPHALOCELE, na.rm=TRUE),
    SB = sum(ALL_1_1_3_SPINA_BIFIDA, na.rm=TRUE),
    NTD_no_gen = sum(ALL_13_GENETIC_CONDITIONS==0, na.rm=TRUE))

All_NTD_num_pregend_wide <- All_NTD_num_pregend %>%
  pivot_wider(
    names_from=pregnancy_end_type,
    values_from= c(NTD_all, anen, ence, SB, NTD_no_gen),
    id_cols=c(maternal_age_group, simd, sex, year)
  )

#Merge together numerator datasets (i.e. all NTDs, NTDs by subgroup and NTDs by pregnancy outcome)
Num <- full_join(All_NTD_num, All_NTD_num_pregend_wide, by = c("year", "sex", "simd", "maternal_age_group"))

# ---- Section 4: Join numerator and denominator data and save file ----

final_data <- full_join(Denom, Num, by = c("year", "sex", "simd", "maternal_age_group"))
final_data <- final_data %>%
  mutate(across(where(is.numeric), ~replace_na(., 0)))

#Check numbers are as we would expect
sum(final_data$count_all)
sum(final_data$count_lb)
sum(final_data$NTD_all)
sum(final_data$anen)
sum(final_data$ence)
sum(final_data$SB)
sum(final_data$NTD_no_gen)

####Save dataset for further analysis
write.csv(final_data, "/conf/Congenital-Descriptive-Epi/Trends in NTDs/Data/analysis_ready.csv", row.names = TRUE)

# ---- Section 5: Prepare analysis file for descriptive analysis of NTDs without an associated genetic condition ----

All_NTD <- SLICCD %>%
  filter(ALL_1_1_NEURAL_TUBE_DEFECTS == 1 & ALL_13_GENETIC_CONDITIONS==0) %>%
  select(ALL_1_1_NEURAL_TUBE_DEFECTS, 
         ALL_1_1_1_ANENCEPHALUS, 
         ALL_1_1_2_ENCEPHALOCELE, 
         ALL_1_1_3_SPINA_BIFIDA, 
         time_period, 
         pregnancy_end_type, 
         mother_age,
         simd,
         sex) %>%
  mutate(maternal_age_group=cut(
    mother_age,
    breaks = c(0, 20, 25, 30, 35, 40, Inf),
    labels = c("<20", "20-24", "25-29", "30-34", "35-39", "40+"),
    right = FALSE)) %>%
  mutate(
    sex= case_when(sex==1 ~ "M",
                   sex==2 ~ "F",
                   sex==3 ~ "Unknown"),
    pregnancy_end_type = case_when(pregnancy_end_type=="Livebirth" ~ "lb",
                                   pregnancy_end_type=="Late Fetal Death" ~ "sb", 
                                   pregnancy_end_type=="Stillbirth" ~ "sb",
                                   pregnancy_end_type=="TOPFA" ~ "top")) %>%
  rename(year=time_period)

All_NTD$simd <- ifelse(is.na(All_NTD$simd), "Unknown", All_NTD$simd) 
All_NTD$maternal_age_group <- ifelse(is.na(All_NTD$maternal_age_group), "Unknown", All_NTD$maternal_age_group)


#Create aggregate count data which can be merged to denominator data

#Calculate total number of NTDs
All_NTD_num <- All_NTD %>%
  group_by(year, maternal_age_group, simd, sex) %>%
  summarise(
    NTD_all = sum(ALL_1_1_NEURAL_TUBE_DEFECTS, na.rm=TRUE),
    anen = sum(ALL_1_1_1_ANENCEPHALUS, na.rm=TRUE),
    ence = sum(ALL_1_1_2_ENCEPHALOCELE, na.rm=TRUE),
    SB = sum(ALL_1_1_3_SPINA_BIFIDA, na.rm=TRUE))

#Calculate NTD totals by pregnancy outcome 
All_NTD_num_pregend <- All_NTD %>%
  group_by(year, pregnancy_end_type, maternal_age_group, simd, sex) %>%
  summarise(
    NTD_all = sum(ALL_1_1_NEURAL_TUBE_DEFECTS, na.rm=TRUE),
    anen = sum(ALL_1_1_1_ANENCEPHALUS, na.rm=TRUE),
    ence = sum(ALL_1_1_2_ENCEPHALOCELE, na.rm=TRUE),
    SB = sum(ALL_1_1_3_SPINA_BIFIDA, na.rm=TRUE))

All_NTD_num_pregend_wide <- All_NTD_num_pregend %>%
  pivot_wider(
    names_from=pregnancy_end_type,
    values_from= c(NTD_all, anen, ence, SB),
    id_cols=c(maternal_age_group, simd, sex, year)
  )

#Merge together numerator datasets (i.e. all NTDs, NTDs by subgroup and NTDs by pregnancy outcome)
Num <- full_join(All_NTD_num, All_NTD_num_pregend_wide, by = c("year", "sex", "simd", "maternal_age_group"))

final_data <- Num %>%
  mutate(across(where(is.numeric), ~replace_na(., 0)))

#Save file for descriptive analysis
write.csv(final_data, "/conf/Congenital-Descriptive-Epi/Trends in NTDs/Data/sub_nogenetic_analysis_ready.csv", row.names = TRUE)

# ---- End of file ----
