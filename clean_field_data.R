

library(dplyr)
library(readr)
library(gsheet)

####LOAD DATA---------------------------------

#load NWA and Halifax harbour receiver metadata
NW_meta <-read_csv("~/Documents/MSc/Ongoing Papers/Mackerel /northwestarm_recmetadata.csv") %>% 
  slice (-17) %>% 
  mutate(station = if_else(station == "HFXAPPROACHES003", "HFX-Approaches003", station))

station_order <- c("OTN-NWA-001", "OTN-NWA-002", "OTN-NWA-003", "OTN-NWA-004", 
                   "OTN-NWA-005", "OTN-NWA-006", "OTN-NWA-007", "OTN-NWA-008", 
                   "HFX-Approaches001", "HFX-Approaches002", "HFX-Approaches003",
                   "OTN-NWA-009", "OTN-NWA-010", "NSCC_OTN_Georges", "OTN-NSCC-MacDonald",
                   "OTN-NSCC-MacKay", "OTN-NSCC-Eastern Passage")

NW_meta2 <- NW_meta %>% 
  mutate(station_number = match(station, station_order))
view(NW_meta2)

saveRDS(NW_meta2, file = "NW_meta_clean.rds")


#load detection data
NW_det <- readRDS("~/Documents/MSc/Ongoing Papers/Mackerel /NWArm_Mackerel (2).rds") %>%
  dplyr::select(-speed, -llon, -llat, -lag_station, -lagdt, -dist, -time_diff) %>%
  dplyr::rename(station = r) %>%
  dplyr::inner_join(NW_meta2, by = "station") %>%
  mutate(oid_short = substr(oid, nchar(oid) - 3, nchar(oid))) #shorten fish ID

#Filtering out false detections 
# Calculate the minimum time difference (min_lag) between consecutive detections for each fish at each station
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

# View the cleaned dataset
head(NW_det_filt)

saveRDS(NW_det_filt, file = "NW_det_filt.rds")



#load tag metadata
tagmeta<-gsheet2tbl('https://docs.google.com/spreadsheets/d/1l8XHcmFLQvQJExbapCyESCTvF-1U7OqrWFQvYsfQSB8/edit?gid=1359375544#gid=1359375544')

tagmeta_filtered <- tagmeta %>%
  filter(grepl("mackerel", COMMON_NAME_E, ignore.case = TRUE)) %>%  #keep rows with "mackerel"
  filter(RELEASE_LOCATION == "Northwest Arm, Halifax, NS") %>% 
  mutate(
    codespace = TAG_CODE_SPACE,   #rename 
    tagid = TAG_ID_CODE,          #rename 
    oid = paste(TAG_CODE_SPACE, TAG_ID_CODE, sep = "-"),  #create oid by combining two columns
    TL = `LENGTH (m)` * 100  #convert TL from m to cm
  ) %>%
  select(
    TL,
    ht = SURGERY_DURATION,
    Project,                   # Keep "Project" as is
    oid,                       #include the new concatenated column
    UTC_RELEASE_DATE_TIME,  #keep "UTC_RELEASE_DATE_TIME" as is
    RELEASE_LOCATION
  )

saveRDS(tagmeta_filtered, file = "tagmeta_clean.rds")
