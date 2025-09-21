# electric-mackerel2.0
 Code space for Brown et al. 2025 "Surgical Tagging of Atlantic Mackerel (Scomber scombrus): Electroanesthesia and Survival in Captivity and the Field"

This repository contains the code and data used in the analysis and figures for the paper "Surgical Tagging of Atlantic Mackerel (Scomber scombrus): Electroanesthesia and Survival in Captivity and the Field" by Brown et al. (2025). 
The study investigates the use of electroanesthesia for surgical tagging of Atlantic mackerel and evaluates the survival rates in both captive and field conditions.

All code works on R version 4.4.2.

Download the zipfile of the repository. Some of the files need to be ran before others. 
If you are interested in the analysis for the lab portion of the study, start with the
"mackerel_lab" file. This includes the logistic regressions and survival analyses for the captive portion of the study.

"Wound" includes the code for the wound condition boxplot included in the supplementary materials.

If you are interested in the field portion of the study, "false_det_calc" includes the code 
used to identify and remove false detections from the acoustic telemetry data, as well as
the code used to calculate the proportion of false detections

"survival_abacus", includes the process of identifying mortality events based on stationary horizontal detections 
and lack of detections outside of the Northwest Arm study area. 
It also includes the code used to create the first and final draft of the abacus plot. 

To run "lab_field_surv_analysis", you will need to have run "mackerel_lab" first.
"lab_field_surv_analysis" includes the Kaplan-Meier survival analysis for both the lab portion of the study, 
the log rank test, and code for the combined survival curves of the lab and field portions of the study.






