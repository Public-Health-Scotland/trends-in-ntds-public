###### TRENDS IN NEURAL TUBE DEFECTS (NTD) R CODE SCRIPT

# This script compares cases of NTDs in SLICCD and in CARDRISS

# ---- Section 1: Housekeeping ----

# clear the global environment
rm(list = ls())

# load required packages
require(tidyverse)
require(ggplot2)
require(janitor)

# Read in the SlICCD AND CARDRISS 2021 data
dfc <- readRDS(paste0("/conf/Congenital-Descriptive-Epi/Trends in NTDs/Data/Cohort 2 - CARDRISS/cardriss_extract.rds"))
dfs <- readRDS(paste0("/conf/Congenital-Descriptive-Epi/Trends in NTDs/Data/sliccd_all_cc.rds"))


# ---- Section 2: Comparison of all NTDS ----

# Cases classified as NTD on SLICCD
dfs_ntd <- dfs %>% filter(ALL_1_1_NEURAL_TUBE_DEFECTS == 1) %>% 
  filter(time_period == 2021) 
table(dfs_ntd$pregnancy_end_type)

# Cases classified as NTD on SLICCD but not on CARDRISS
dfc1 <- dfc %>% filter(ALL_1_1_NEURAL_TUBE_DEFECTS == 0) 
table(dfc1$outcome_of_pregnancy)

# Cases classified as NTD on CARDRISS
dfc2 <- dfc %>% filter(ALL_1_1_NEURAL_TUBE_DEFECTS == 1)
table(dfc2$outcome_of_pregnancy)

# On CARDRISS with an NTD but not classified as NTD on SLICCD
#SLICCD but not coded as NTD
dfc2 <- dfc %>% select(sliccd_id)
# number not on SLICCD at all
dfc2 <- dfc2 %>% filter(sliccd_id != is.na(sliccd_id))
# merge datasets to see how many on SLICCD but not as NTD
dfc2 <- dfc2 %>% mutate(sliccd_id2 = sliccd_id)
dfmerge <- left_join(dfc2, dfs, by = "sliccd_id")
dfmerge <- dfmerge %>% filter(ALL_1_1_NEURAL_TUBE_DEFECTS == 0)

# Input values for sensitivity and Positive Predictive Value (PPV) calculation from the above
TP <- 47  # True Positives
FP <- 1   # False Positives
FN <- 22  # False Negatives
TN <- 0   # True Negatives

# Sensitivity calculation
sensitivity <- TP / (TP + FN)
sensitivity_ci <- binom.test(TP, TP + FN)$conf.int

# PPV calculation
ppv <- TP / (TP + FP)
ppv_ci <- binom.test(TP, TP + FP)$conf.int

# Output results
cat("Sensitivity: ", sensitivity, "\n")
cat("95% CI for Sensitivity: [", sensitivity_ci[1], ", ", sensitivity_ci[2], "]\n")
cat("PPV: ", ppv, "\n")
cat("95% CI for PPV: [", ppv_ci[1], ", ", ppv_ci[2], "]\n")

# ---- Section 3: Comparison of anencephaly ----

# On CARDRISS with anencephaly
dfc2 <- dfc %>% filter(ALL_1_1_1_ANENCEPHALUS == 1)
table(dfc2$outcome_of_pregnancy)

# Cases classified as anencephaly on SLICCD
dfs1 <- dfs %>% filter(ALL_1_1_1_ANENCEPHALUS == 1) %>% 
  filter(time_period == 2021) 
table(dfs1$pregnancy_end_type)

# Input values for calculating sensitivity and PPV
TP <- 26  # True Positives
FP <- 0   # False Positives
FN <- 11  # False Negatives
TN <- 0   # True Negatives

# Sensitivity calculation
sensitivity <- TP / (TP + FN)
sensitivity_ci <- binom.test(TP, TP + FN)$conf.int

# PPV calculation
ppv <- TP / (TP + FP)
ppv_ci <- binom.test(TP, TP + FP)$conf.int

# Output results
cat("Sensitivity: ", sensitivity, "\n")
cat("95% CI for Sensitivity: [", sensitivity_ci[1], ", ", sensitivity_ci[2], "]\n")
cat("PPV: ", ppv, "\n")
cat("95% CI for PPV: [", ppv_ci[1], ", ", ppv_ci[2], "]\n")

# ---- Section 4: Comparison of spina bifida ----

# On CARDRISS with spina bifida
dfc2 <- dfc %>% filter(ALL_1_1_3_SPINA_BIFIDA == 1)
table(dfc2$outcome_of_pregnancy)

# Cases classified as spina bifida on SLICCD
dfs1 <- dfs %>% filter(ALL_1_1_3_SPINA_BIFIDA == 1) %>% 
  filter(time_period == 2021) 
table(dfs1$pregnancy_end_type)

# On CARDRISS with an NTD but not classified as NTD on SLICCD
# SLICCD but not coded as NTD
dfc2 <- dfc2 %>% select(sliccd_id)
# number not on SLICCD at all
dfc2 <- dfc2 %>% filter(sliccd_id != is.na(sliccd_id))
# merge datasets to see how many on SLICCD but not as NTD
dfc2 <- dfc2 %>% mutate(sliccd_id2 = sliccd_id)
dfmerge <- left_join(dfc2, dfs1, by = "sliccd_id")
dfmerge <- dfmerge %>% filter(ALL_1_1_3_SPINA_BIFIDA == 1)

# Input values based on the above
TP <- 20  # True Positives
FP <- 1   # False Positives
FN <- 7  # False Negatives
TN <- 0   # True Negatives

# Sensitivity calculation
sensitivity <- TP / (TP + FN)
sensitivity_ci <- binom.test(TP, TP + FN)$conf.int

# Positive predictive value (PPV) calculation
ppv <- TP / (TP + FP)
ppv_ci <- binom.test(TP, TP + FP)$conf.int

# Output results
cat("Sensitivity: ", sensitivity, "\n")
cat("95% CI for Sensitivity: [", sensitivity_ci[1], ", ", sensitivity_ci[2], "]\n")
cat("PPV: ", ppv, "\n")
cat("95% CI for PPV: [", ppv_ci[1], ", ", ppv_ci[2], "]\n")

# ---- Section 5: Comparison of NTDs with an associated genetic condition ----

# create dataframe for SLICCD and CARDRISS by NTD and genetic condition status
data <- dfc %>% 
  # selecting reg id, sliccd id and the relevant EUROCAT groups
  select(registration_id, 
         sliccd_id,
         # renamed to ensure we know it is the groups on CARDRISS
         cardriss_ntd = ALL_1_1_NEURAL_TUBE_DEFECTS,
         cardriss_genetic = ALL_13_GENETIC_CONDITIONS) %>% 
  # joining on the relevant sliccd records and their EUROCAT groups for NTD and genetic
  left_join(dfs %>% 
              select(sliccd_id, 
                     # renamed to ensure we know it is the groups on SLiCCD
                     sliccd_ntd = ALL_1_1_NEURAL_TUBE_DEFECTS,
                     sliccd_genetic = ALL_13_GENETIC_CONDITIONS), 
            by = "sliccd_id") %>% 
  # filling in sliccd groups with 0 in rows where there is no associated sliccd record
  mutate(sliccd_ntd = replace_na(sliccd_ntd, 0),
         sliccd_genetic = replace_na(sliccd_genetic, 0))

# number NTD without GC on both SLiCCD and CARDRISS
truepos <- data %>% 
    filter(cardriss_ntd == 1 & cardriss_genetic == 0 & sliccd_ntd == 1
           & sliccd_genetic == 0)

# total number NTD without GC on CARDRISS
totalcard <- data %>% 
  filter(cardriss_ntd == 1 & cardriss_genetic == 0)

# total number NTD without GC on SLiCCD
totalslic <- data %>% 
  filter(sliccd_ntd == 1
         & sliccd_genetic == 0)

# number with no NTD recorded on CARDRISS with NTD without GC on SLiCCD
falseposN1 <- data %>% 
  filter(cardriss_ntd == 0 & sliccd_ntd == 1
         & sliccd_genetic == 0)

# number with NTD with GC recorded on CARDRISS with NTD without GC on SLiCCD
falseposN2 <- data %>% 
  filter(cardriss_ntd == 1 & cardriss_genetic == 1 & sliccd_ntd == 1
         & sliccd_genetic == 0)

# number NTD without GC on CARDRISS with no NTD on SLiCCD
falsenegN1 <- data %>% 
  filter(cardriss_ntd == 1 & cardriss_genetic == 0 & sliccd_ntd == 0)

# number NTD without GC on CARDRISS with NTD with GC on SLiCCD
falsenegN2 <- data %>% 
  filter(cardriss_ntd == 1 & cardriss_genetic == 0 & sliccd_ntd == 1
         & sliccd_genetic == 1)

# Input values
TP <- 44  # True Positives
FP <- 2   # False Positives
FN <- 23  # False Negatives
TN <- 0   # True Negatives

# Sensitivity calculation
sensitivity <- TP / (TP + FN)
sensitivity_ci <- binom.test(TP, TP + FN)$conf.int

# PPV calculation
ppv <- TP / (TP + FP)
ppv_ci <- binom.test(TP, TP + FP)$conf.int

# Output results
cat("Sensitivity: ", sensitivity, "\n")
cat("95% CI for Sensitivity: [", sensitivity_ci[1], ", ", sensitivity_ci[2], "]\n")
cat("PPV: ", ppv, "\n")
cat("95% CI for PPV: [", ppv_ci[1], ", ", ppv_ci[2], "]\n")

# ---- End of file ----
