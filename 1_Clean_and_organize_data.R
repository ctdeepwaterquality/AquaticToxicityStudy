library(tidyverse)
library(lubridate)
library(psych)
library(pastecs)
library(janitor)

#rm(list = ls())


### Reading in datasets ###
## all chronic bio data
chronic <- read.csv("chronic_bio_2023update.csv", na.strings=c("","NA"))

## all summary river data
river <- read.csv("chronic_river_2023update.csv", na.strings=c("","NA")) %>% janitor::clean_names() %>% 
  mutate(date = as.Date(date, format = "%m/%d/%Y"))

## < or > 100 summary of chronic results
chronic_summary <- read.csv("chronic_summary_2023update.csv", na.strings=c("","NA")) %>% 
  janitor::clean_names() %>%
  mutate(date = as.Date(date, format = "%m/%d/%Y")) %>%
  mutate(across(everything(), gsub, pattern = "%", replacement = "")) %>% 
  mutate(across(c(6:10), gsub, pattern = " ", replacement = "")) 

## all acute bio data
acute <- read.csv("acute_bio_2023update.csv", na.strings=c("","NA"))




## only river names and npdes
river_name <- river %>% 
  select(facility, npdes) %>% 
  distinct()

river2 <- river %>% 
  mutate(invert_toxic = ifelse(is.na(invert_toxic), NA, paste("i", invert_toxic, sep = "_"))) %>% 
  mutate(vert_toxic = ifelse(is.na(vert_toxic), NA, paste("v", vert_toxic, sep = "_"))) %>% 
  unite(toxic, c(invert_toxic, vert_toxic)) %>% 
  filter(toxic != "NA_NA") %>% 
  mutate(across(toxic, gsub, pattern = "NA", replacement = "")) %>%
  select(date, npdes, toxic)

  # select(facility, npdes, date, invert_toxic, vert_toxic) %>% 
  # drop_na(invert_toxic, vert_toxic)
  

## lab control: test is invalid if lab control does not meet test acceptable criteria
invalid_control <- chronic_summary %>% 
  mutate(date = as.Date(date)) %>% 
  filter(grepl('Lab Control', facility)) %>% 
  filter(if_any(everything(), ~ grepl("<", .)))


## chronic results, effluent only
chronic2 <- chronic %>% 
  janitor::clean_names() %>% # cleans up the column header names
  mutate(date = as.Date(date, format = "%m/%d/%Y")) %>% 
  filter(!grepl('Lab Control | River | Creek | Brook | Dam', sites)) %>% 
  filter(grepl('WPCF', sites)) %>% 
  anti_join(invalid_control, by = c("npdes", "date")) %>% # remove tests from dataset that have invalid controls
  anti_join(river2, by = c("npdes", "date")) %>% # remove tests that had some toxicity in river shown in river tests
  # left_join(river_name, by = "npdes") %>% 
  rename(invert = invertebrate_spp, 
         vert = fish_spp, 
         invert_survival = invertebrate_survival, 
         vert_survival = fish_survival) %>% 
  select(-c(x3_brood_per_or_fecunditiy,
            young_per_female,
            mysid_growth_mg,
            fish_growth_mg)) %>% 
  filter(!is.na(as.numeric(invert_survival))) %>% # remove those with no survival info
  mutate(chronic_acute = "chronic",
         vert_survival = as.character(vert_survival)) %>% 
  pivot_longer(cols = invert_survival:vg_noec) %>% 
  mutate(across(value, gsub, pattern = "%", replacement = ""))



## acute dataset; convert survival to LC50 values
 acute2 <- acute %>% 
  janitor::clean_names() %>% 
  mutate(date = as.Date(testdate, format = "%m/%d/%Y")) %>% 
  filter(!if_any(note, ~ grepl("invalid", .,ignore.case = TRUE))) %>% 
  select(sites = stp,
         npdes,
         date,
         invert_survival = dsurv, 
         vert_survival = fsurv) %>%
   mutate(i_lc50 = ifelse(invert_survival < 50, "<100",
                         ifelse(invert_survival > 50, ">100", NA)),
         v_lc50 = ifelse(vert_survival < 50, "<100",
                         ifelse(vert_survival > 50, ">100", NA))) %>%
   mutate(invert_survival = as.character(invert_survival),
          vert_survival = as.character(vert_survival),
          chronic_acute = "acute") %>%
   pivot_longer(cols = invert_survival:v_lc50)


merge <- chronic2 %>% 
  bind_rows(acute2)

write.csv(merge, "chronic_acute_merged.csv")
