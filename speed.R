#load packages 
library(ggplot2)
library(tidyverse)
library(dplyr)
library(readr)
library(lubridate)
library(gsheet)

####LOAD DATA
#load receiver metadata 
NW_rec <-read_csv("~/Documents/MSc/Ongoing Papers/Mackerel /northwestarm_recmetadata.csv") %>% 
  filter(project == "Northwest Arm") %>% 
  mutate(station_number = row_number())

#load tag metadata
tagmeta <-gsheet2tbl('https://docs.google.com/spreadsheets/d/1l8XHcmFLQvQJExbapCyESCTvF-1U7OqrWFQvYsfQSB8/edit?gid=2134606216#gid=2134606216')

tagmeta_filt <- tagmeta %>%
  filter(grepl("mackerel", COMMON_NAME_E, ignore.case = TRUE)) %>%  
  mutate(
    codespace = TAG_CODE_SPACE,   
    tagid = TAG_ID_CODE,          
    oid = paste(TAG_CODE_SPACE, TAG_ID_CODE, sep = "-"),
    TL = as.numeric(`TOTAL_LENGTH (m)`) * 100  # Convert TL from meters to cm
  ) %>%
  rename(ht = SURGERY_DURATION) %>%
  dplyr::select(
    TL,  # Now in cm
    Project,                  
    oid,
    ht
  )

tagmeta_filt <- tagmeta_filtered %>%
  filter(grepl("ARM-2024", Project, ignore.case = TRUE))  # Filter for project NW-2024

#load detection data
NW_det <- readRDS("~/Documents/MSc/Ongoing Papers/Mackerel /NWArm_Mackerel (2).rds")
NW_det <- NW_det %>% select (-speed, -llon, -llat, -lag_station, -lagdt, -dist, -time_diff) %>%
  rename(station = r) %>% 
  inner_join(NW_rec, by = "station")

####FILTERING FALSE DETECTIONS
#calculate the minimum time difference (min_lag) between consecutive detections for each fish at each station
NW_det_filt <- NW_det %>%
  arrange(oid, station, dt) %>%  #ensure data is sorted by fish ID, station, and datetime
  group_by(oid, station) %>%
  mutate(min_lag = as.numeric(difftime(dt, lag(dt), units = "mins"))) %>%  #compute time difference in minutes
  ungroup() %>%
  
  #group by fish, station, and detection hour
  group_by(oid, station, detection_hour = floor(as.numeric(dt) / (60 * 60))) %>%
  mutate(count_per_hour = n()) %>%
  ungroup() %>%
  
  #filter false detections: keep detections if they occur more than once in an hour 
  #OR if there's a previous detection within 60 minutes
  mutate(passed_filter = ifelse(count_per_hour == 1 & (is.na(min_lag) | min_lag > 60), 0, 1)) %>%
  filter(passed_filter == 1) %>% #only keep the detections that passed the filter
  select(-count_per_hour, -detection_hour, -min_lag, -passed_filter) #remove temporary columns

#view the cleaned dataset
head(NW_det_filt)


####MERGE DETECTION DATA WITH TAG METADATA
NW_det_filt<-merge(NW_det_filt, tagmeta_filt, by="oid")

#data check 
str(NW_det_filt)

####INDIVIDUAL ABACUS PLOT FOR EACH FISH COLOURED BY SPEED

#calculating speed for each fish
armfish_speed <-NW_det_filt %>%
  arrange(dt) %>%                               #arrange by datetime
  group_by(oid) %>%                             #group by 'oid'
  mutate(
    llon = lag(longitude),                      #previous longitude
    llat = lag(latitude),                       #previous latitude
    lag_station = lag(station),                 #previous station
    lagdt = lag(dt),                            #previous timestamp
    time_diff = as.numeric(difftime(dt, lagdt, units = "mins")) #time difference in minutes
  ) %>%
  filter(!is.na(llon) & !is.na(lagdt)) %>%      #remove rows with missing lagged values
  rowwise() %>%                                 #perform operations row-wise
  mutate(
    dist = argosfilter::distance(llat, latitude, llon, longitude), #calculate distance row-wise
    speed = dist / (time_diff / 60)             #speed (distance divided by time in hours)
  ) %>%
  ungroup() %>%                                 #ungroup to return data frame
  filter(speed <= 30) 


library(purrr)

#convert dt to Date format
armfish_speed <- armfish_speed %>%
  mutate(dt = as.Date(dt))  #ensures dt is in date format

min_dt <- min(armfish_speed$dt, na.rm = TRUE)
max_dt <- max(armfish_speed$dt, na.rm = TRUE)


#split data by fish ID and create individual plots
mack_abacus <- armfish_speed %>%
  as_tibble() %>%
  mutate(
    oid = gsub("A69-9007-", "", oid),
    station_number = factor(station_number, levels = rev(sort(unique(station_number))))  # Standardized y-axis station order
  ) %>%
  split(.$oid) %>%
  map(function(df) {  
    df <- df %>% arrange(speed)  # Order by speed to layer properly
    
    ggplot(df, aes(x = dt, y = station_number, color = speed)) +
      geom_point(size = 4) +  
      scale_color_gradientn(
        colors = c("blue", "cyan", "green", "yellow", "orange", "red", "darkred"),
        values = scales::rescale(c(0, 5, 10, 15, 20, 25, 30)),  
        limits = c(0, 30),  
        name = "Speed (km/h)"
      ) +
      scale_x_date(limits = c(as.Date(min_dt), as.Date(max_dt)+ 5),
                   date_breaks = "2 weeks",
                   date_labels = "%b %d") +
      scale_y_discrete(drop = FALSE) +  # Ensures all stations are displayed
      labs(
        title = paste("Fish ID:", unique(df$oid)),
        x = "Date",
        y = "Receiver Station"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.key.size = unit(0.8, "cm"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        title = element_text(size = 16)
      )
  })

mack_abacus

# Open PDF device
pdf("arm_abacus_speed.pdf", width = 7, height = 5) 

# Loop through the plots and print them to the PDF
for (plot in mack_abacus) {
  print(plot)
}

# Close the PDF device
dev.off()