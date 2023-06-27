library(tidyverse)
library(lubridate)
library(janitor)



########################
#### IMPORT DATASET ####
########################


merge <- read.csv("chronic_acute_merged.csv") %>% 
  select(-c(invert, vert, X)) %>% 
  mutate(date = as.Date(date, format = "%Y-%m-%d"),
         year = year(date)) %>% 
  filter(year >= 2008)



######################
#### INVERTEBRATE ####
######################

#### SURVIVAL ##################################################################

## Invert survival percent
subset_survival <- merge %>% 
  filter(name == "invert_survival") %>% 
  mutate(value = as.numeric(value)) %>% 
  arrange(npdes)


## This creates a df that has only the npdes sites with both acute and chronic data
crossover <- subset_survival %>% 
  select(npdes, chronic_acute) %>% 
  distinct(npdes, chronic_acute) %>%
  mutate(list = 1) %>% 
  pivot_wider(names_from = chronic_acute, values_from = list) %>% 
  drop_na(acute, chronic) %>% 
  select(npdes)


## assign pass/fail based on criteria (80%+ for chronic, 90%+ for acute)
## count # of pass/fail per year per acute/chronic test (chronic should be 1, acute will be multiple)
i_survival_pf <- subset_survival %>% 
  semi_join(crossover, by = "npdes") %>% 
  select(npdes, year, chronic_acute, value) %>% 
  mutate(pass_fail = ifelse(chronic_acute == "acute", 
                            ifelse(value >= 90, "pass", "fail"),
                            ifelse(value >= 80, "pass", "fail"))) %>%
  group_by(npdes, year, chronic_acute, pass_fail) %>% 
  summarise(count = n()) %>% 
  ungroup() %>% 
  drop_na() %>% 
  unite(test_pass, chronic_acute, pass_fail, sep = "_") %>% 
  pivot_wider(names_from = test_pass, values_from = count) %>% 
  select(npdes, year, acute_pass, acute_fail, chronic_pass, chronic_fail) %>% 
  mutate(acute_fail = ifelse(!is.na(acute_fail), acute_fail, 
                             ifelse(!is.na(acute_pass), 0, NA)),
         chronic_fail = ifelse(!is.na(chronic_fail), chronic_fail, 
                             ifelse(!is.na(chronic_pass), 0, NA)),
         chronic_pass = ifelse(!is.na(chronic_pass), chronic_pass, 
                               ifelse(!is.na(chronic_fail), 0, NA))) %>% 
  mutate(survival_discrepancy = ifelse(
    (is.na(acute_pass) & is.na(acute_fail)) | (is.na(chronic_pass) & is.na(chronic_fail)),
    NA, 
    ifelse((acute_pass > 0) & (chronic_fail >= 1),
           "Yes",
           ifelse((chronic_pass > 0) & (acute_fail >= 1),
                  "Yes",
                  "No"))))


#### LC50 ######################################################################

## Keep only the LC50 data for inverts
## Assign pass/fail based on LC50 endpoint results
subset_lc50 <- merge %>% 
  filter(name == "i_lc50") %>% 
  mutate(value = ifelse(value == ">100", "pass", 
                        ifelse(value == "100", "pass", "fail")))

## Creates a list of sites that have both acute and chronic results at some point in the 10 years sampling period
filtered_sites_lc50 <- subset_lc50 %>% 
  group_by(npdes, chronic_acute) %>% 
  summarize(n = n()) %>% 
  arrange(npdes, chronic_acute) %>% 
  pivot_wider(names_from = chronic_acute, values_from = n) %>% 
  filter(!is.na(chronic) & !is.na(acute)) %>% 
  distinct(npdes)

## Filters the df to only the sites with both chronic and acute data
i_lc50_pf <- subset_lc50 %>% 
  semi_join(filtered_sites_lc50, by = "npdes") %>%
  group_by(npdes, year, chronic_acute, value) %>%
  summarise(count = n()) %>% 
  ungroup() %>% 
  drop_na() %>% 
  unite(test_pass, chronic_acute, value, sep = "_") %>% 
  pivot_wider(names_from = test_pass, values_from = count) %>% 
  select(npdes, year, acute_pass, acute_fail, chronic_pass, chronic_fail) %>% 
  mutate(acute_fail = ifelse(!is.na(acute_fail), acute_fail, 
                             ifelse(!is.na(acute_pass), 0, NA)),
         chronic_fail = ifelse(!is.na(chronic_fail), chronic_fail, 
                               ifelse(!is.na(chronic_pass), 0, NA)),
         chronic_pass = ifelse(!is.na(chronic_pass), chronic_pass, 
                               ifelse(!is.na(chronic_fail), 0, NA))) %>% 
  mutate(lc50_discrepancy = ifelse(
    (is.na(acute_pass) & is.na(acute_fail)) | (is.na(chronic_pass) & is.na(chronic_fail)),
    NA, 
    ifelse((acute_pass > 0) & (chronic_fail >= 1),
           "Yes",
           ifelse((chronic_pass > 0) & (acute_fail >= 1),
                  "Yes",
                  "No"))))


#### EC50 ######################################################################

## Keep only the EC50 data for inverts
subset_ec50 <- merge %>% 
  filter(name == "i_ec50") %>% 
  mutate(value = ifelse(value == ">100", "pass", 
                        ifelse(value == "100", "pass", "fail")))

## Filters the df to only the sites with both chronic and acute data
i_ec50_pf <- subset_ec50 %>% 
  group_by(npdes, year, chronic_acute, value) %>%
  summarise(count = n()) %>% 
  ungroup() %>% 
  drop_na() %>% 
  unite(test_pass, chronic_acute, value, sep = "_") %>% 
  pivot_wider(names_from = test_pass, values_from = count) %>% 
  rename(ec50_chronic_pass = chronic_pass,
         ec50_chronic_fail = chronic_fail)
  
  
#### COMBINE ECO50 & LC50 ######################################################

i_lc50_ec50_pf <- i_lc50_pf %>% 
  left_join(i_ec50_pf, by = c("npdes", "year")) %>% 
  mutate(ec50_chronic_fail = ifelse(!is.na(ec50_chronic_fail), ec50_chronic_fail,
                               ifelse(!is.na(ec50_chronic_pass), 0, NA)),
         ec50_chronic_pass = ifelse(!is.na(ec50_chronic_pass), ec50_chronic_pass,
                               ifelse(!is.na(ec50_chronic_fail), 0, NA))) %>% 
  mutate(ec50_discrepancy = ifelse(
    (is.na(acute_pass) & is.na(acute_fail)) | (is.na(ec50_chronic_pass) & is.na(ec50_chronic_fail)),
    NA,
    ifelse((acute_pass > 0) & (ec50_chronic_fail >= 1),
           "Yes",
           ifelse((ec50_chronic_pass > 0) & (acute_fail >= 1),
                  "Yes",
                  "No")))) %>% 
  rename(lc50_chronic_pass = chronic_pass,
         lc50_chronic_fail = chronic_fail)


## final tables
i_survival_pf
i_lc50_ec50_pf



######################
####     FISH     ####
######################

#### SURVIVAL ##################################################################

## Fish survival percent
subset_survival <- merge %>% 
  filter(name == "vert_survival") %>% 
  mutate(value = as.numeric(value)) %>% 
  arrange(npdes)


## This creates a df that has only the npdes sites with both acute and chronic data
crossover <- subset_survival %>% 
  select(npdes, chronic_acute) %>% 
  distinct(npdes, chronic_acute) %>%
  mutate(list = 1) %>% 
  pivot_wider(names_from = chronic_acute, values_from = list) %>% 
  drop_na(acute, chronic) %>% 
  select(npdes)


## assign pass/fail based on criteria (80%+ for chronic, 90%+ for acute)
## count # of pass/fail per year per acute/chronic test (chronic should be 1, acute will be multiple)
v_survival_pf <- subset_survival %>% 
  semi_join(crossover, by = "npdes") %>% 
  select(npdes, year, chronic_acute, value) %>% 
  mutate(pass_fail = ifelse(chronic_acute == "acute", 
                            ifelse(value >= 90, "pass", "fail"),
                            ifelse(value >= 80, "pass", "fail"))) %>%
  group_by(npdes, year, chronic_acute, pass_fail) %>% 
  summarise(count = n()) %>% 
  ungroup() %>% 
  drop_na() %>% 
  unite(test_pass, chronic_acute, pass_fail, sep = "_") %>% 
  pivot_wider(names_from = test_pass, values_from = count) %>% 
  select(npdes, year, acute_pass, acute_fail, chronic_pass, chronic_fail) %>% 
  mutate(acute_fail = ifelse(!is.na(acute_fail), acute_fail, 
                             ifelse(!is.na(acute_pass), 0, NA)),
         chronic_fail = ifelse(!is.na(chronic_fail), chronic_fail, 
                               ifelse(!is.na(chronic_pass), 0, NA)),
         chronic_pass = ifelse(!is.na(chronic_pass), chronic_pass, 
                               ifelse(!is.na(chronic_fail), 0, NA))) %>% 
  mutate(survival_discrepancy = ifelse(
    (is.na(acute_pass) & is.na(acute_fail)) | (is.na(chronic_pass) & is.na(chronic_fail)),
    NA, 
    ifelse((acute_pass > 0) & (chronic_fail >= 1),
           "Yes",
           ifelse((chronic_pass > 0) & (acute_fail >= 1),
                  "Yes",
                  "No"))))


#### LC50 ######################################################################

## Keep only the LC50 data for inverts
## Assign pass/fail based on LC50 endpoint results
subset_lc50 <- merge %>% 
  filter(name == "v_lc50") %>% 
  mutate(value = ifelse(value == ">100", "pass", 
                        ifelse(value == "100", "pass", "fail")))

## Creates a list of sites that have both acute and chronic results at some point in the 10 years sampling period
filtered_sites_lc50 <- subset_lc50 %>% 
  group_by(npdes, chronic_acute) %>% 
  summarize(n = n()) %>% 
  arrange(npdes, chronic_acute) %>% 
  pivot_wider(names_from = chronic_acute, values_from = n) %>% 
  filter(!is.na(chronic) & !is.na(acute)) %>% 
  distinct(npdes)

## Filters the df to only the sites with both chronic and acute data
v_lc50_pf <- subset_lc50 %>% 
  semi_join(filtered_sites_lc50, by = "npdes") %>%
  group_by(npdes, year, chronic_acute, value) %>%
  summarise(count = n()) %>% 
  ungroup() %>% 
  drop_na() %>% 
  unite(test_pass, chronic_acute, value, sep = "_") %>% 
  pivot_wider(names_from = test_pass, values_from = count) %>% 
  select(npdes, year, acute_pass, acute_fail, chronic_pass, chronic_fail) %>% 
  mutate(acute_fail = ifelse(!is.na(acute_fail), acute_fail, 
                             ifelse(!is.na(acute_pass), 0, NA)),
         chronic_fail = ifelse(!is.na(chronic_fail), chronic_fail, 
                               ifelse(!is.na(chronic_pass), 0, NA)),
         chronic_pass = ifelse(!is.na(chronic_pass), chronic_pass, 
                               ifelse(!is.na(chronic_fail), 0, NA))) %>% 
  mutate(lc50_discrepancy = ifelse(
    (is.na(acute_pass) & is.na(acute_fail)) | (is.na(chronic_pass) & is.na(chronic_fail)),
    NA, 
    ifelse((acute_pass > 0) & (chronic_fail >= 1),
           "Yes",
           ifelse((chronic_pass > 0) & (acute_fail >= 1),
                  "Yes",
                  "No"))))


#### EC50 ######################################################################

## Keep only the EC50 data for inverts
subset_ec50 <- merge %>% 
  filter(name == "v_ec50") %>% 
  mutate(value = ifelse(value == ">100", "pass", 
                        ifelse(value == "100", "pass", "fail")))

## Filters the df to only the sites with both chronic and acute data
v_ec50_pf <- subset_ec50 %>% 
  group_by(npdes, year, chronic_acute, value) %>%
  summarise(count = n()) %>% 
  ungroup() %>% 
  drop_na() %>% 
  unite(test_pass, chronic_acute, value, sep = "_") %>% 
  pivot_wider(names_from = test_pass, values_from = count) %>% 
  rename(ec50_chronic_pass = chronic_pass,
         ec50_chronic_fail = chronic_fail)


#### COMBINE ECO50 & LC50 ######################################################

v_lc50_ec50_pf <- v_lc50_pf %>% 
  left_join(v_ec50_pf, by = c("npdes", "year")) %>% 
  mutate(ec50_chronic_fail = ifelse(!is.na(ec50_chronic_fail), ec50_chronic_fail,
                                    ifelse(!is.na(ec50_chronic_pass), 0, NA)),
         ec50_chronic_pass = ifelse(!is.na(ec50_chronic_pass), ec50_chronic_pass,
                                    ifelse(!is.na(ec50_chronic_fail), 0, NA))) %>% 
  mutate(ec50_discrepancy = ifelse(
    (is.na(acute_pass) & is.na(acute_fail)) | (is.na(ec50_chronic_pass) & is.na(ec50_chronic_fail)),
    NA,
    ifelse((acute_pass > 0) & (ec50_chronic_fail >= 1),
           "Yes",
           ifelse((ec50_chronic_pass > 0) & (acute_fail >= 1),
                  "Yes",
                  "No")))) %>% 
  rename(lc50_chronic_pass = chronic_pass,
         lc50_chronic_fail = chronic_fail)


## final tables
v_survival_pf
v_lc50_ec50_pf




invert <- i_survival_pf %>% 
  left_join(i_lc50_ec50_pf, by = c("npdes", "year")) %>% 
  select(npdes, 
         year, 
         survival_discrepancy,
         lc50_discrepancy, 
         ec50_discrepancy) %>% 
  rename(i_survival_discrepancy = survival_discrepancy,
         i_lc50_discrepancy = lc50_discrepancy, 
         i_ec50_discrepancy = ec50_discrepancy) 

invert_disc <- invert %>%
  filter_at(vars(c(3:5)), all_vars(!is.na(.))) %>% 
  filter_at(vars(c(3:5)), any_vars(. == "Yes"))

invert_multiples <- invert_disc %>% 
  group_by(npdes) %>% 
  summarize(n_invert = n()) %>% 
  filter(n_invert > 1)


fish <- v_survival_pf %>% 
  left_join(v_lc50_ec50_pf, by = c("npdes", "year")) %>% 
  select(npdes, 
         year, 
         survival_discrepancy,
         lc50_discrepancy, 
         ec50_discrepancy) %>% 
  rename(v_survival_discrepancy = survival_discrepancy,
         v_lc50_discrepancy = lc50_discrepancy, 
         v_ec50_discrepancy = ec50_discrepancy) 

fish_disc <- fish %>%
  filter_at(vars(c(3:5)), all_vars(!is.na(.))) %>% 
  filter_at(vars(c(3:5)), any_vars(. == "Yes"))

fish_multiples <- fish_disc %>% 
  group_by(npdes) %>% 
  summarize(n_fish = n()) %>% 
  filter(n_fish > 1)

multiples <- invert_multiples %>% 
  left_join(fish_multiples, by = "npdes")




combine <- invert %>%
  left_join(fish, by = c("npdes", "year"))



#### save tables ####
# write.csv(i_survival_pf, "outputs/invert_survival.csv")
# write.csv(i_lc50_ec50_pf, "outputs/invert_lc50_ec50.csv")
# write.csv(v_survival_pf, "outputs/fish_survival.csv")
# write.csv(v_lc50_ec50_pf, "outputs/fish_lc50_ec50_pf.csv")
# write.csv(fish_disc, "outputs/fish_discrepancies.csv")
# write.csv(invert_disc, "outputs/invert_discrepancies.csv")
# write.csv(multiples, "outputs/site_multiple_discrepancies.csv")
















  
  
  
  
  
  