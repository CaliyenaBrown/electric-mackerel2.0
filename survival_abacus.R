#load packages 
library(tidyverse)    #ggplot2, dplyr, readr, etc.
library(lubridate)  
library(gsheet)       
library(rnaturalearth)      
library(rnaturalearthhires) 


####LOAD DATA---------------------------------
det <- readRDS("NW_det_filt.rds")
tagmeta <- readRDS("tagmeta_clean.rds")
NW_meta <- readRDS("NW_meta_clean.rds")


#####ABACUS---------------------------------
#abacus plot of all fish over time coloured by station number
#get the last detection per fish
last_detections <- det %>%
  group_by(oid) %>%
  filter(dt == max(dt)) %>%
  ungroup() %>%
  mutate(
    day = as.Date(dt)
  )

#make sure station_number is a factor (for consistent color mapping)
det$station_number <- factor(det$station_number)
last_detections$station_number <- factor(last_detections$station_number)

#potential mortality events from stationary detections
#Get the true last detection per fish
last_dt_per_fish <- det %>%
  group_by(oid) %>%
  summarise(last_dt = max(dt), .groups = "drop")

#Group consecutive detections at the same receiver
det1 <- det %>%
  arrange(oid, dt) %>%
  group_by(oid) %>%
  mutate(
    same_as_last = station == lag(station, default = first(station)),
    group_id = cumsum(!same_as_last)
  ) %>%
  ungroup()

#Summarize sequences of same-station detections
stationary_events <- det1 %>%
  group_by(oid, group_id, station) %>%
  summarise(
    start_time = min(dt),
    end_time = max(dt),
    duration_hr = as.numeric(difftime(max(dt), min(dt), units = "hours")),
    n_detections = n(),
    .groups = "drop"
  ) %>%
  filter(n_detections >= 5 & duration_hr >= 24)

#join with the last detection per fish
stationary_events_final <- stationary_events %>%
  left_join(last_dt_per_fish, by = "oid") %>%
  filter(end_time == last_dt)

potential_mortality_ids <- unique(stationary_events_final$oid)

#highlight mortality events
#create a new column indicating mortality status in last detections
last_detections <- last_detections %>%
  mutate(mortality = ifelse(oid %in% potential_mortality_ids, "Mortality", "Survived"))

#ensure station number and mortality are factors
det$station_number <- factor(det$station_number)
last_detections$station_number <- factor(last_detections$station_number)
last_detections$mortality <- factor(last_detections$mortality, levels = c("Survived", "Mortality"))

####FIRST ABACUS PLOT-----------------
ggplot(det, aes(x = dt, y = oid_short)) +
  geom_point(aes(colour = station_number), size = 3.5) + # All detections
  geom_point(  # Final detections for mortalities (black fill, colored border)
    data = filter(last_detections, mortality == "Mortality"),
    aes(x = dt, y = oid_short, colour = station_number),
    shape = 21, fill = "black", size = 3.5, stroke = 1,
    show.legend = FALSE
  ) +
  geom_text(
    data = last_detections,
    aes(x = dt + lubridate::days(4), y = oid_short, label = station_number),
    size = 3.5,
    color = "black"
  ) +
  labs(
    x = "Date",
    y = "Fish ID",
    colour = "Receiver Station"
  ) +
  scale_colour_viridis_d(option = "plasma") +
  scale_x_datetime(date_labels = "%b %d", breaks = "5 days") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.key.size = unit(0.6, "cm"),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "right",
    plot.background = element_rect(fill = "white", color = "black")
  )

#save plot
ggsave("abacus_plot_mortality.png", last_plot(), width = 11, height = 8, 
       units = "in", dpi = 300)



####CHECK FOR DETECTIONS OUTSIDE NWA---------------------------------
#identify and filter detections occurring outside the NWA study array

#rec meta from otn
geoserver_receivers <-readr::read_csv('https://members.oceantrack.org/geoserver/otn/ows?service=WFS&version=1.0.0&request=GetFeature&typeName=otn:stations_receivers&outputFormat=csv', guess_max = 13579)

#keep only the most recent download per station and select useful columns
geoserver_receivers_filtered<-geoserver_receivers %>% 
  group_by(station_name) %>%
  slice_max(order_by = last_download, n = 1) %>%
  dplyr::select(lat=stn_lat, long=stn_long, station=station_name, project=collectioncode)

#map to check location of receivers
atlantic <- rnaturalearthhires::countries10 %>% 
  dplyr::filter(NAME_IT == "Canada" | grepl("America", NAME_IT)) %>%
  ggplot() +
  geom_point(data = geoserver_receivers_filtered, aes(x = long, y = lat), color = "red", size = 2) +
  geom_sf()+
  coord_sf(xlim=c(-53.0, -71),
           ylim=c(53.0, 40))+
  theme_classic(base_size=25)
atlantic


#load detections from OTN
dets_out <- read_csv("v2lnasti_matched_detections_2024.csv") %>% 
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

#filter to only include potential moralities 
dets_out_filt <- dets_out_cleaned %>%
  filter(oid %in% potential_mortality_ids)

#remove all detections from NWA receivers
dets_out_filt <- dets_out_filt %>%
  filter(!station %in% NW_meta$station)

#confirm that all detections are from outside the study site
unique(dets_out_filt$station)

#define the stations you want to remove
stations_to_remove <- c("APPROACHES03", "APPROACHES02", "APPROACHES01", 
                        "MACDONALD", "MACKAY", "HFX001")

#filter them out of the detections
dets_out_filt <- dets_out_filt %>%
  filter(!station %in% stations_to_remove)

#update potential mortalities to not include any oid found outside the study site 
#oids with detections outside the study site
oids_to_remove <- unique(dets_out_filt$oid)

#remove them from potential mortalities
potential_mortality_ids <- setdiff(potential_mortality_ids, oids_to_remove)

#how many days was this final fish detected as receiver 10?
target_oid <- "A69-1303-8980"

# Filter to detections for this fish at station 10
receiver10_days <- det %>%
  filter(oid == target_oid, station_number == 10) %>%
  mutate(date = as.Date(dt)) %>%  #extract date from datetime
  distinct(date) %>%              #keep unique dates
  nrow()                          #count unique dates

receiver10_days #24 days


####FINAL ABACUS PLOT-----------------
ggplot(det, aes(x = dt, y = oid_short)) +
  geom_point(aes(colour = station_number), size = 3.5) +
  geom_point(
    data = filter(last_detections, oid %in% potential_mortality_ids),
    aes(x = dt, y = oid_short, colour = station_number),
    shape = 21, fill = "black", size = 3.5, stroke = 1,
    show.legend = FALSE) +
  geom_text(data = last_detections,
    aes(x = dt + lubridate::days(4), y = oid_short, label = station_number),
    size = 3.5, color = "black") +
    labs(
          x = "Date",
          y = "Fish ID",
          colour = "Receiver Station") +
  scale_colour_viridis_d(option = "plasma") +
  scale_x_datetime(date_labels = "%b %d", breaks = "5 days") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.key.size = unit(0.6, "cm"),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "right",
    plot.background = element_rect(fill = "white", color = "black")
  )

ggsave("abacus_plot_mortality.png", last_plot(), width = 12, height = 8, 
       units = "in", dpi = 300)

