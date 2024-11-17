###### TRENDS IN NEURAL TUBE DEFECTS (NTD) R CODE SCRIPT

# ---- Section 1: Housekeeping ----

# clear the global environment
rm(list = ls())

# load required packages
require(tidyverse)
require(ggplot2)
require(lubridate)
require(pwr)

# ---- Section 2: Calculate cumulative number of births & NTDs for 10 years following implementation ----

power <- data.frame(year=1:10) %>%
  mutate(births=year*48306,
         n_ntd=round(births*(11.4/10000)),
         prev=(n_ntd/births),
         births_baseline=193224,
         n_ntd_baseline=276,
         prev_baseline=14.3/10000)

# ---- Section 3: Conduct power calculations ----

#Create function which does the power calculation based on providing data on the total number in the
#exposed group, the total number in the unexposed group, the proportion in the exposed group and 
#the proportion in the unexposed group
calculate_power <- function(exposed_total, exposed_proportion, unexposed_total, unexposed_proportion) {
  effect_size <- ES.h(exposed_proportion, unexposed_proportion)
  power_result <- pwr.2p2n.test(h=effect_size,
                              n1=exposed_total,
                              n2=unexposed_total,
                              sig.level = 0.05,
                              alternative="two.sided")  
                              
                              return(power_result$power)
}

#Apply the function to each row of the dataset created in section 2
power$power_est <- mapply(calculate_power,
                          power$births,
                          power$prev,
                          power$births_baseline,
                          power$prev_baseline)

#Round the estimate of power to two decimal places
power <- power %>%
  mutate(power_est=round(power_est, 2))

#Save the dataframe
write.csv(power, "/conf/Congenital-Descriptive-Epi/Trends in NTDs/Data/power.csv", row.names = TRUE)

# ---- End of file ----
