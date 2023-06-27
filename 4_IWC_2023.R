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
  filter(year >= 2008,
         chronic_acute == "chronic") %>% 
  filter(name %in% c("is_loec", ## filter for just NOEL and LOEL values
                     "is_noec",
                     "ir_loec",
                     "ir_noec",
                     "vs_loec",
                     "vs_noec",
                     "vg_loec",
                     "vg_noec")) %>% 
  select(npdes, year, name, value)

iwc <- read.csv("iwc_2023.csv") %>% 
  clean_names() %>% 
  select(npdes, iwc)


#######################
#### DATA ANALYSIS ####
#######################


merge_iwc <- merge %>% 
  left_join(iwc, by = "npdes") %>% 
  drop_na(value) %>% 
  filter(!value == ">100") %>% 
  mutate(value = ifelse(value == "<100", 99.9, value)) %>% 
  mutate(value = str_remove_all(value, ">|<"),
         value = str_remove_all(value, "/100|100/"),
         value = ifelse(grepl("6.25", value), 6.25, value),
         value = ifelse(grepl("12.5", value), 12.5, value),
         value = str_remove_all(value, "-|/"),
         value = as.numeric(value))

noec <- merge_iwc %>% 
  filter(grepl("noec", name)) %>% 
  mutate(impact = ifelse(value >= iwc, "not impacted", NA)) %>% 
  select(npdes, year, name, impact) %>% 
  drop_na(impact) %>% 
  pivot_wider(names_from = name, values_from = impact) %>% 
  mutate(across(everything(), ~na_if(.x, "NULL"))) %>% 
  unnest(cols = c(ir_noec, is_noec, vg_noec, vs_noec))



loec <- merge_iwc %>% 
  filter(grepl("loec", name)) %>% 
  mutate(impact = ifelse(value <= iwc, "impacted", NA)) %>% 
  select(npdes, year, name, impact) %>% 
  drop_na(impact) %>% 
  pivot_wider(names_from = name, values_from = impact) %>% 
  mutate(across(everything(), ~na_if(.x, "NULL"))) %>% 
  unnest(cols = c(vs_loec, vg_loec, ir_loec, is_loec))

noec_summary <- noec %>% 
  group_by(npdes) %>% 
  summarize(not_impacted = n())

loec_summary <- loec %>% 
  group_by(npdes) %>% 
  summarize(impacted = n())

summary <- noec_summary %>% 
  left_join(loec_summary, by = "npdes")

# write.csv(noec, "outputs/iwc_noec.csv")
# write.csv(loec, "outputs/iwc_loec.csv")
# write.csv(summary, "outputs/iwc_summary.csv")
