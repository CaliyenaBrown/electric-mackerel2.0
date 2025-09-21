
library(dplyr)
library(survival)
library(survminer)

####LOAD DATA---------------------------------
#load cleaned detections, tag metadata, and station metadata
det <- readRDS("NW_det_filt.rds")
tagmeta <- readRDS("tagmeta_clean.rds")
NW_meta <- readRDS("NW_meta_clean.rds")


####SURVIVAL ANALYSIS---------------------------------
#survival analysis of tagged fish using detection data

#study period
study_start <- as.Date("2024-07-24")
study_end <- as.Date("2024-11-11")

#extract tagged and detected fish ids
tagged_oids <- unique(tagmeta$oid)
detected_oids <- unique(det$oid)


####UNDETECTED FISH (day 0 mortalities)-----------------
#fish that were tagged but never detected are assumed mortalities on day 0
undetected_oids <- setdiff(tagged_oids, detected_oids)

undetected_fish_df <- data.frame(
  oid = undetected_oids,
  date_of_death = rep(study_start, length(undetected_oids)),
  fate = 0,     #0=dead
  event = 1     #1=mortality event
) %>%
  mutate(int = as.numeric(difftime(date_of_death, study_start, units = "days")))


####STATIONARY MORTALITY (known individual)-----------------
#example case: one fish (oid A69-1303-8980) identified as stationary for extended period
stationary_oid <- "A69-1303-8980"

#extract detections for this fish
stationary_fish <- det %>%
  filter(oid == stationary_oid) %>%
  arrange(dt) %>%
  mutate(date = as.Date(dt))

#count detections per day
daily_counts <- stationary_fish %>%
  count(date)

#identify consecutive days of detection
daily_counts <- daily_counts %>%
  arrange(date) %>%
  mutate(
    diff = as.integer(date - lag(date, default = first(date))),
    group = cumsum(coalesce(diff != 1, TRUE))
  )

#find first stationary run ≥ 5 days
stationary_run <- daily_counts %>%
  group_by(group) %>%
  filter(n() >= 5) %>%
  summarise(mortality_start = max(date), .groups = "drop") %>%
  slice(1)

#backdate mortality to be conservative (10 days earlier)
stationary_run$mortality_start <- stationary_run$mortality_start - 10

#create mortality dataframe
stationary_fish_df <- data.frame(
  oid = stationary_oid,
  date_of_death = stationary_run$mortality_start,
  fate = 0,
  event = 1
) %>%
  mutate(int = as.numeric(difftime(date_of_death, study_start, units = "days")))


####SURVIVORS (censored)-----------------
#fish that survived until study_end (no mortality observed)
survivor_oids <- setdiff(tagged_oids, c(undetected_oids, stationary_oid))

survivor_df <- det %>%
  filter(oid %in% survivor_oids) %>%
  group_by(oid) %>%
  summarise(last_dt = max(as.Date(dt)), .groups = "drop") %>%
  mutate(
    fate = 1,   #1=alive
    event = 0,  #0=censored
    int = as.numeric(difftime(study_end, study_start, units = "days"))
  )


####COMBINE ALL GROUPS-----------------
#merge undetected fish, stationary mortality, and survivors into one dataset
final_df <- bind_rows(
  undetected_fish_df %>% select(oid, fate, event, int),
  stationary_fish_df %>% select(oid, fate, event, int),
  survivor_df %>% select(oid, fate, event, int)
)

#confirm mortality count
print(table(final_df$event)) #expect: 0=survivors, 1=mortalities


####SURVIVAL ANALYSIS PLOT-----------------
#kaplan-meier survival estimate
ms_field <- survfit(Surv(int, event) ~ 1, data = final_df)

#plot survival curve
p_field <- ggsurvplot(
  ms_field,
  censor = TRUE,
  palette = "darkblue"
) +
  theme_survminer(font.x = c(16), font.y = c(16), font.legend = c(14)) +
  labs(x = "time (days)", y = "survival probability")

print(p_field)



####SURVIVAL CODE FOR LAB MACKEREL-----------------
#NEED TO LOAD MACKEREL LAB FILE AND RUN SCRIPT BEFORE NEXT STEPS

ms_lab <- survfit(Surv(int, m) ~ 1, data = d)

#Plot survival curve for all mackerel
p_lab <- ggsurvplot(ms_lab, 
                    censor = T, 
                    palette = "darkgreen") +
  theme_survminer(font.x = c(16), font.y = c(16), font.legend = c(14)) +
  labs(x = "Time (Days)", y = "Survival Probability")
p_lab

#combine the two survival plots
#use survival objects created for field and lab mackerel

p_combined <- ggsurvplot_combine(
  list(ms_field, ms_lab),
  data = list(last_detection, d),
  legend.title = "",
  legend.labs = c("Field", "Lab"),
  censor = T,
  conf.int = F,
  palette = c("darkblue", "forestgreen")
) +
  theme_survminer(font.x = c(16), font.y = c(16), font.legend = c(16)) +
  labs(x = "Time (Days)", y = "Survival Probability") 
p_combined

ggsave("survival_analysis_overall2.png", last_plot(), width = 8, height = 6, 
       units = "in", dpi = 300)

##LOG-RANK TEST-----------------
#prepare field data (event already correctly defined where 0 = dead, 1 = alive/censored)
field_data <- final_df %>%
  dplyr::select(int, event) %>%
  dplyr::mutate(group = "Field")

#prepare lab data (convert m to event so that 0 = dead, 1 = alive/censored)
lab_data <- d %>%
  select(int, m) %>%
  mutate(
    event = ifelse(m == 2, 1, 0),  # 2 means death (event=1), 1 means survived (event=0)
    group = "Lab"
  ) %>%
  select(int, event, group)

##CHECK-----------------
table(field_data$event)  #should show few deaths (1s)
table(lab_data$event)    #should show more deaths (1s)

#combine the datasets
combined_data <- bind_rows(field_data, lab_data)

#log-rank test
logrank_test <- survdiff(Surv(int, event) ~ group, data = combined_data)
logrank_test



####CHECKING SURVIVAL USING MATCHED DETECTIONS ON RECEIVERS OUTSIDE OF THE NWA---------------

#load detections from OTN
dets_out <- read_csv("~/Documents/MSc/Ongoing Papers/Mackerel /v2lnasti_matched_detections_2024.csv") %>%
  dplyr::select(tagname, station, datecollected, longitude, latitude) %>%
  rename(oid = tagname, date = datecollected, lon = longitude, lat = latitude)

#merge with clean tag metadata
dets_out <- merge(dets_out, tagmeta, by="oid")

#prepare NWA station metadata
NW_meta <- NW_meta %>% 
  dplyr::select(meta_lon = long, meta_lat = lat) %>%
  distinct()

#remove any detections that exactly match the lat/lon of NW_meta stations
dets_out_cleaned <- dets_out %>%
  anti_join(NW_meta, by = c("lon" = "meta_lon", "lat" = "meta_lat"))

#get list of unique fish (oids) detected outside the study site
fish_detected_outside <- unique(dets_out_cleaned$oid)

#view or write to file if needed
View(fish_detected_outside)
length(fish_detected_outside)
#[33]

#list of mortality fish to check
mortality_ids <- c('A69-1303-8979', 'A69-1303-8980', 'A69-1303-8982', 'A69-1303-9010')

#check which, if any, were detected outside the study site
detected_outside_mortality <- mortality_ids[mortality_ids %in% fish_detected_outside]

#display results
if (length(detected_outside_mortality) == 0) {
  cat("✅ None of the 4 mortality fish were detected outside the study site.\n")
} else {
  cat("⚠️ The following mortality fish were detected outside the study site:\n")
  print(detected_outside_mortality)
}



####PLOT DETECTIONS OUTSIDE STUDY SITE-----------------
#base map: North Atlantic region
atlantic_base <- rnaturalearthhires::countries10 %>%
  filter(NAME_IT == "Canada" | grepl("America", NAME_IT)) %>%
  ggplot() +
  geom_sf() +
  coord_sf(xlim = c(-53.0, -71), ylim = c(53.0, 40)) +  # Adjust to match your study region
  theme_classic(base_size = 25)
atlantic_base

#plot detections outside the study site (from cleaned dets_out)
outside_det <- atlantic_base +
  geom_point(
    data = dets_out_cleaned,aes(x = lon, y = lat, color = oid),size = 2.5, alpha = 0.8
  ) +
  labs(x = "Longitude",y = "Latitude", color = "Fish ID"
  ) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 14)
  )

outside_det


####UPDATED SURVIVAL ANALYSIS INCORPORATING DETECTIONS OUTSIDE STUDY SITE-----------------

#define study start and end
study_start <- as.Date("2024-07-24")
study_end <- as.Date("2024-11-11")

#all tagged fish
tagged_oids <- unique(tagmeta$oid)
detected_oids <- unique(det$oid)

#fish never detected = Day 0 mortality
undetected_oids <- setdiff(tagged_oids, detected_oids)

undetected_fish_df <- data.frame(
  oid = undetected_oids,
  date_of_death = rep(study_start, length(undetected_oids)),
  fate = 0,     # dead
  event = 1     # death
) %>%
  mutate(int = as.numeric(difftime(date_of_death, study_start, units = "days")))

#stationary mortality (known individual)
stationary_oid <- "A69-1303-8980"

stationary_fish <- det %>%
  filter(oid == stationary_oid) %>%
  arrange(dt) %>%
  mutate(date = as.Date(dt))

#daily counts
daily_counts <- stationary_fish %>%
  count(date) %>%
  arrange(date) %>%
  mutate(
    diff = as.integer(date - lag(date, default = first(date))),
    group = cumsum(coalesce(diff != 1, TRUE))
  )

#identify first 5+ day stationary run
stationary_run <- daily_counts %>%
  group_by(group) %>%
  filter(n() >= 5) %>%
  summarise(mortality_start = max(date), .groups = "drop") %>%
  slice(1)

#set death date
stationary_fish_df <- data.frame(
  oid = stationary_oid,
  date_of_death = stationary_run$mortality_start,
  fate = 0,  # dead
  event = 1,
  int = as.numeric(difftime(stationary_run$mortality_start, study_start, units = "days"))
)

#survivors (censored)
survivor_oids <- setdiff(tagged_oids, c(undetected_oids, stationary_fish_df$oid))

survivor_df <- det %>%
  filter(oid %in% survivor_oids) %>%
  group_by(oid) %>%
  summarise(last_dt = max(as.Date(dt)), .groups = "drop") %>%
  mutate(
    fate = 1,   # survived
    event = 0,  # censored
    int = as.numeric(difftime(last_dt, study_start, units = "days"))
  )

#combine all into final dataframe
final_df <- bind_rows(
  undetected_fish_df %>% select(oid, fate, event, int),
  stationary_fish_df %>% select(oid, fate, event, int),
  survivor_df %>% select(oid, fate, event, int)
)

#confirm breakdown
print(table(final_df$event))  # event: 1 = deaths, 0 = censored

#kaplan–Meier survival curve (field)
ms_field <- survfit(Surv(int, event) ~ 1, data = final_df)

p_field <- ggsurvplot(
  ms_field,
  censor = TRUE,
  palette = "darkblue"
) +
  theme_survminer(font.x = c(16), font.y = c(16), font.legend = c(14)) +
  labs(x = "Time (Days)", y = "Survival Probability")

print(p_field)

#lab dataset survival (from 'd')
ms_lab <- survfit(Surv(int, m) ~ 1, data = d)

p_lab <- ggsurvplot(ms_lab, 
                    censor = TRUE, 
                    palette = "darkgreen") +
  theme_survminer(font.x = c(16), font.y = c(16), font.legend = c(14)) +
  labs(x = "Time (Days)", y = "Survival Probability")

print(p_lab)

#combine field vs lab
p_combined <- ggsurvplot_combine(
  list(ms_field, ms_lab),
  data = list(final_df, d),
  legend.title = "",
  legend.labs = c("Field", "Lab"),
  censor = TRUE,
  conf.int = FALSE,
  palette = c("darkblue", "forestgreen")
) +
  theme_survminer(font.x = c(16), font.y = c(16), font.legend = c(16)) +
  labs(x = "Time (Days)", y = "Survival Probability")

print(p_combined)

#log-rank test
field_data <- final_df %>%
  select(int, event) %>%
  mutate(group = "Field")

lab_data <- d %>%
  select(int, m) %>%
  mutate(
    event = case_when(
      m == 0 ~ 1,  # dead
      m == 2 ~ 1,  # alternate coding
      m == 1 ~ 0,  # alive
      TRUE ~ NA_integer_
    ),
    group = "Lab"
  ) %>%
  select(int, event, group)

combined_data <- bind_rows(field_data, lab_data)

logrank_test <- survdiff(Surv(int, event) ~ group, data = combined_data)
print(logrank_test)
