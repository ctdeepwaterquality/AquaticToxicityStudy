library(tidyverse)
library(lubridate)
library(janitor)

### Does out copper criteria protect against acute and chronic toxicity?

### This analysis looks at the copper data (dissolved and total) for the rivers in the chronic toxicity
### tests. Datasets used include chronic test endpoints and chronic chemistry results. We are looking at 
### Vertebrate acute, vertebrate chronic, invertebrate acute, and invertebrate chronic toxicity
### results and comparing the copper results from these tests to the WQ criteria. 

######################
#### READ IN DATA ####
######################

## All summary river data
river_sub <- read.csv("chronic_river_2023update.csv", na.strings=c("","NA")) %>%
  janitor::clean_names() %>%
  mutate(date = as.Date(date, format = "%m/%d/%Y")) 
  

## < or > 100 summary of all chronic results
chronic_summary <- read.csv("chronic_summary_2023update.csv", na.strings=c("","NA")) %>% 
  janitor::clean_names() %>%
  mutate(date = as.Date(date, format = "%m/%d/%Y")) %>%
  mutate(across(everything(), gsub, pattern = "%", replacement = "")) %>% 
  mutate(across(c(6:10), gsub, pattern = " ", replacement = "")) 

## Toxicity dataset
tox <- read.csv("chronic_acute_merged.csv")

## All chronic bio data
chronic <- read.csv("chronic_bio_2023update.csv", na.strings=c("","NA")) %>% 
  janitor::clean_names() %>%
  mutate(date = as.Date(date, format = "%m/%d/%Y"))

## Municipal chemistry
chem_mun <- read.csv("chronic_chem_wgt_2023update3.csv", na.strings=c("","NA")) %>% 
  janitor::clean_names() %>% 
  mutate(date = as.Date(date, format = "%m/%d/%Y"))




########################
#### TRANSFORM DATA ####
########################

## Saltwater sites
sw <- tox %>% 
  filter(chronic_acute == "chronic",
         invert == "Mb") %>% 
  distinct(npdes)

## Only sites/dates that have invalid controls
invalid_control_mun <- chronic_summary %>% 
  mutate(date = as.Date(date)) %>% 
  filter(grepl('Lab Control', facility)) %>% 
  filter(if_any(everything(), ~ grepl("<", .))) %>% 
  select(npdes, date)

## Select columns in rivers dataset
river_sub <- river_sub %>% 
  select(facility, date, npdes, i_r, v_g)


## Copper data
copper_mun <- chem_mun %>% 
  select(sites, ss, npdes, date, cu_diss, cu_tot) %>% 
  rename(cu_diss_ugl = cu_diss,
         cu_tot_ugl = cu_tot) %>% 
  slice(-1) %>% 
  anti_join(invalid_control_mun, by = c("npdes", "date")) %>% ## remove with invalid lab control
  anti_join(sw, by = "npdes") %>%  ## remove saltwater sites
  filter(year(date) %in% 2010:2020,
         !if_all(c(cu_diss_ugl, cu_tot_ugl), is.na)) %>%
  arrange(npdes, date) %>% 
  mutate(year = year(date),
         month = month(date))
  # select(npdes, date, month, year, cu_diss_ugl, cu_tot_ugl)
  
  # summarize_copper <- copper_mun %>%
  #   group_by(npdes, year, month) %>%
  #   summarise(n = n())
  



## Chronic data, subsetted for just rivers
chronic_rivers <- chronic %>%
  left_join(river_sub, by = c("npdes", "date")) %>% 
  filter(grepl('River|Creek|Brook|Dam|Harbor', sites),
         year(date) %in% 2010:2020,
         invertebrate_spp == "Cd") %>% 
  mutate(across(everything(), gsub, pattern = "%", replacement = ""),
         across(c(4:10), gsub, pattern = " ", replacement = ""),
         across(everything(), gsub, pattern = "NT", replacement = NA, fixed = TRUE),
         year = year(date),
         month = month(date)) %>% 
  select(sites,
         npdes,
         date,
         month,
         year,
         invert_survival = invertebrate_survival,
         invert_brood = x3_brood_per_or_fecunditiy, 
         invert_young = young_per_female,
         i_r_noec = i_r,
         fish_survival,
         fish_growth_mg,
         v_g_noec = v_g) %>% 
  mutate(site_specific = NA,
    a = ifelse(as.numeric(invert_survival) >= 90, "not toxic", "toxic"),
    b = ifelse(as.numeric(invert_brood) >= 6, "not toxic", "toxic"),
    c = ifelse(as.numeric(invert_young) >= 15, "not toxic", "toxic"),
    d = ifelse(i_r_noec == "<100", "toxic", 
                            ifelse(is.na(i_r_noec), NA, "not toxic")),
    x = ifelse(as.numeric(fish_survival) >= 90, "not toxic", "toxic"),
    y = ifelse(as.numeric(fish_survival) >= 0.23, "not toxic", "toxic"),
    z = ifelse(v_g_noec == "<100", "toxic", 
                         ifelse(is.na(i_r_noec), NA, "not toxic"))) %>% 
  mutate(date = as.Date(date)) %>% 
  left_join(copper_mun, by = c("npdes", "month", "year"))




box_df <- chronic_rivers %>% 
  mutate(invert_acute = a,
         invert_chronic = ifelse(((b == "not toxic" | is.na(b)) & (c == "not toxic" | is.na(c)) & (d == "not toxic" | is.na(d))),
                                 "not toxic", "toxic"),
         fish_acute = x,
         fish_chronic = ifelse(((y == "not toxic" | is.na(y)) & (z == "not toxic" | is.na(z))),
                               "not toxic", "toxic")) %>% 
  select(c(npdes,
           ss,
           invert_acute,
           invert_chronic,
           fish_acute,
           fish_chronic,
           cu_diss_ugl, 
           cu_tot_ugl)) %>% 
  mutate(cu_diss_ugl = as.numeric(cu_diss_ugl),
         cu_tot_ugl = as.numeric(cu_tot_ugl))

# write.csv(box_df, "outputs/copper_boxplot_dataset.csv")
# write.csv(chronic_rivers, "outputs/chronic_river_copper_data.csv")



###################
#### BOXPOLOTS ####
###################

names(box_df)

##### DISSOLVED COPPER #######

cat.labs <- c("Invertebrate Acute",
              "Invertebrate Chronic",
              "Vertebrate Acute",
              "Vertebrate Chronic")

names(cat.labs) <- c("invert_acute",
                     "invert_chronic",
                     "fish_acute",
                     "fish_chronic")


### Dissolved copper

box_df %>% 
  pivot_longer(cols = c(3:6)) %>% 
  drop_na(value) %>% 
  ggplot(aes(x = cu_diss_ugl, y = value)) +
  geom_boxplot() +
  geom_vline(xintercept = 14.3, color = "red") + #acute
  geom_vline(xintercept = 4.8, color = "orange") + #chronic
  facet_wrap(~name, 
             labeller = labeller(name = cat.labs),
             nrow=4)+
  labs(x = "Dissolved Copper (ug/L)", 
       y = "Toxicity",
       title = "Dissolved Copper Toxicity") + 
  theme_light()

# ggsave("outputs/diss_copper.png", width = 8.5, height = 11, units = "in")

### Total copper

box_df %>% 
  pivot_longer(cols = c(3:6)) %>% 
  drop_na(value) %>% 
  ggplot(aes(x = cu_tot_ugl, y = value)) +
  geom_boxplot() +
  geom_vline(xintercept = 14.3, color = "red") + #acute
  geom_vline(xintercept = 4.8, color = "orange") + #chronic
  facet_wrap(~name, 
             labeller = labeller(name = cat.labs),
             nrow=4)+
  labs(x = "Total Copper (ug/L)", 
       y = "Toxicity",
       title = "Total Copper Toxicity") + 
  theme_light()

# ggsave("outputs/tot_copper.png", width = 8.5, height = 11, units = "in")



# Site specific rivers/criteria
### Dissolved copper

box_df %>% 
  filter(ss == "ss") %>% 
  pivot_longer(cols = c(3:6)) %>% 
  drop_na(value) %>% 
  ggplot(aes(x = cu_diss_ugl, y = value)) +
  geom_boxplot() +
  geom_vline(xintercept = 25.7, color = "red") + #acute
  geom_vline(xintercept = 18.1, color = "orange") + #chronic
  facet_wrap(~name, 
             labeller = labeller(name = cat.labs),
             nrow=4)+
  labs(x = "Dissolved Copper (ug/L)", 
       y = "Toxicity",
       title = "Dissolved Copper Toxicity") + 
  theme_light()

# ggsave("outputs/ss_diss_copper.png", width = 8.5, height = 11, units = "in")

### Total copper

box_df %>% 
  filter(ss == "ss") %>% 
  pivot_longer(cols = c(3:6)) %>% 
  drop_na(value) %>% 
  ggplot(aes(x = cu_tot_ugl, y = value)) +
  geom_boxplot() +
  geom_vline(xintercept = 25.7, color = "red") + #acute
  geom_vline(xintercept = 18.1, color = "orange") + #chronic
  facet_wrap(~name, 
             labeller = labeller(name = cat.labs),
             nrow=4)+
  labs(x = "Total Copper (ug/L)", 
       y = "Toxicity",
       title = "Total Copper Toxicity") + 
  theme_light()

# ggsave("outputs/ss_tot_copper.png", width = 8.5, height = 11, units = "in")
