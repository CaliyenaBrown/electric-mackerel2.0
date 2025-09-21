
#load data 
NW_det <- readRDS("NWArm_Mackerel (2).rds")
NW_det_filt <- readRDS("NW_det.rds")

#count number of filered detections 
total_detections <- nrow(NW_det)  # Total original detections
filtered_detections <- nrow(NW_det_filt)  # Detections after filtering

#Calculate number and proportion of false detections
false_detections <- total_detections - filtered_detections
proportion_false <- false_detections / total_detections

cat("False detections filtered out:", false_detections, "\n")
cat("Proportion of false detections:", round(proportion_false * 100, 2), "%\n")