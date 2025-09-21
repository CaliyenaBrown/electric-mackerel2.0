
#load data 
NW_det <- readRDS("NWArm_Mackerel (2).rds") %>%
  dplyr::select(-speed, -llon, -llat, -lag_station, -lagdt, -dist, -time_diff) %>%
  dplyr::rename(station = r) %>%
  dplyr::inner_join(NW_meta2, by = "station") %>%
  mutate(oid_short = substr(oid, nchar(oid) - 3, nchar(oid))) #shorten fish ID


#Filtering out false detections 
#calculate the minimum time difference (min_lag) between consecutive detections for each fish at each station
NW_det_filt <- NW_det %>%
  arrange(oid, station_number, dt) %>%
  group_by(oid, station_number) %>%
  mutate(min_lag = as.numeric(difftime(dt, lag(dt), units = "mins"))) %>%
  ungroup() %>%
  group_by(oid, station_number, detection_hour = floor(as.numeric(dt) / (60 * 60))) %>%
  mutate(count_per_hour = n()) %>%
  ungroup() %>%
  mutate(passed_filter = ifelse(count_per_hour == 1 & (is.na(min_lag) | min_lag > 60), 0, 1)) %>%
  filter(passed_filter == 1) %>%
  dplyr::select(-count_per_hour, -detection_hour, -min_lag, -passed_filter)


#count number of filered detections 
total_detections <- nrow(NW_det)  # Total original detections
filtered_detections <- nrow(NW_det_filt)  # Detections after filtering

#calculate number and proportion of false detections
false_detections <- total_detections - filtered_detections
proportion_false <- false_detections / total_detections

cat("False detections filtered out:", false_detections, "\n")
cat("Proportion of false detections:", round(proportion_false * 100, 2), "%\n")
