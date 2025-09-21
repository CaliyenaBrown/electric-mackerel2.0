#load packages 
library(dplyr)
library(ggplot2)
library(survival)
library(survminer)
library(tidyr)
library(patchwork)
library(boot)


####LOAD DATA---------------------------------
lab <- readRDS("lab_clean.rds")
mort <- readRDS("Mort_date.rds")

#visualize plot results as boxplot of total length
ggplot(lab, aes(factor(m), tl)) + geom_boxplot()

#visualize plot results as boxplot of handling time
ggplot(lab, aes(factor(m), ht)) + geom_boxplot()


####LOGISTIC REGRESSION---------------------------------
#fit logistic regression model 
m1<-glm(m ~ tl + ht, data = lab, family = "binomial") 
summary(m1)

summary_stats <- summary(m1)
summary(lab)

#calculate the odds ratios for the coefficients by exponentiating the coefficient estimates 
#extract coefficient estimates from the summary object
coefficients <- coef(summary(m1))

#extract odds ratios
odds_ratios <- exp(coefficients[, "Estimate"])
print(odds_ratios)

#95% confidence intervals 
#extract coefficient estimates and their confidence intervals
coef <- coef(summary_stats)
conf_int <- confint(m1)

#print the 95% confidence intervals
print(conf_int)

#exponentiate the conf intervals because we exponentiated the odds ratios 
exp_conf_int <- exp(conf_int)
print(exp_conf_int)



#### PREDICTIONS ACROSS TL -------------------------------------
#make a separate copy of lab for TL analysis
lab_tl <- lab

#generate predictions based on model (m1) across TL range
prediction_data_tl <- tibble(tl = 27:39) %>%
  mutate(ht = mean(lab_tl$ht, na.rm = TRUE)) %>%  #use average handling time
  predict(m1, newdata = ., type = "response", se.fit = TRUE)

#extract fitted values and confidence intervals
prediction_data_tl <- tibble(fit = prediction_data_tl$fit, se = prediction_data_tl$se.fit) %>%
  mutate(min = fit - (1.96 * se), max = fit + (1.96 * se)) %>%
  bind_cols(tibble(tl = 27:39))  # Add TL values back

#summarize lab dataset to count fish at each TL
lab_summary_tl <- lab_tl %>%
  group_by(tl) %>%
  summarise(count = n(), .groups = "drop")

#merge summarized count data back for plotting
lab_tl <- left_join(lab_tl, lab_summary_tl, by = "tl")

#convert count into a factor for discrete color mapping (TL uses 1:6 levels)
lab_tl$count <- factor(lab_tl$count, levels = 1:6)


#### LOGISTIC REGRESSION PLOT (TL) -----------------------------
#define custom blue color palette for TL
custom_blues_tl <- c("1" = "navyblue", "2" = "royalblue4", "3" = "royalblue",
                     "4" = "steelblue2", "5" = "dodgerblue", "6" = "lightskyblue")

#define line thickness mapping for TL
size_values_tl <- c(2, 3.2, 4.3, 5.2, 6, 15)

#order `lab_tl` so that larger counts are plotted first, and smaller points appear on top
lab_tl <- lab_tl %>% arrange(desc(count))

#plot logistic regression prediction with confidence intervals for TL
r1 <- ggplot(prediction_data_tl, aes(x = tl, y = fit)) +
  geom_ribbon(aes(ymin = min, ymax = max), alpha = 0.3, fill = "lightgrey") +  # Confidence Interval
  geom_line(size = 0.8, colour = "black") +  # Logistic regression line
  geom_point(data = lab_tl, aes(x = tl, y = m, size = count, color = count)) +  # Raw data points
  scale_size_manual(values = size_values_tl, guide = guide_legend(title = "# of Fish")) +
  scale_color_manual(values = custom_blues_tl, guide = guide_legend(title = "# of Fish")) +
  labs(x = "Total Length (cm)", y = "Survival Probability") +
  scale_x_continuous(breaks = seq(27, 39, by = 1)) +
  scale_y_continuous(limits = c(min(prediction_data_tl$min), max(prediction_data_tl$max)), 
                     breaks = seq(0, 1, by = 0.5)) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 16, colour = "black", margin = margin(t = 10)),
        axis.title.y = element_text(size = 16, colour = "black"),
        axis.text = element_text(size = 12, colour = "black"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))
r1

# Save the plot as PNG
ggsave("total_length_logistic_regression.png", last_plot(), width = 8, height = 6, 
       units = "in", dpi = 300)


#### PROBABILITY THRESHOLD (0.75) ------------------------------
#Solve for TL where predicted probability = 0.75
# Extract coefficients
b0 <- coef(m1)[1]
b1 <- coef(m1)[2]  # TL
b2 <- coef(m1)[3]  # HT

# Define target probability
target_prob <- 0.75
logit_p <- log(target_prob / (1 - target_prob))

# Use mean HT
mean_ht <- mean(lab$ht, na.rm = TRUE)

# Solve for TL
target_tl <- (logit_p - b0 - b2 * mean_ht) / b1
print(target_tl)

#add horizontal and vertical lines to the plot
r1 <- r1 +
  geom_hline(yintercept = 0.75, linetype = "dashed", color = "red") +  # Horizontal line at 0.75
  geom_vline(xintercept = target_tl, linetype = "dashed", color = "red") +  # Vertical line at TL
  annotate("text", x = target_tl + 0.5, y = 0.77, 
           label = paste0("TL = ", round(target_tl, 1), " cm"), 
           color = "red", size = 4.5, hjust = 0)
r1



#### REGRESSION PLOT (HANDLING TIME) ----------------------------
#generate predictions based on model (m1) across HT range
prediction_data_ht <- tibble(ht = seq(min(lab$ht, na.rm = TRUE), max(lab$ht, na.rm = TRUE), length.out = 100)) %>%
  mutate(tl = mean(lab$tl, na.rm = TRUE)) %>%  # Use average total length
  predict(m1, newdata = ., type = "response", se.fit = TRUE)

#extract fitted values and confidence intervals
prediction_data_ht <- tibble(fit = prediction_data_ht$fit, se = prediction_data_ht$se.fit) %>%
  mutate(min = fit - (1.96 * se), max = fit + (1.96 * se)) %>%
  bind_cols(tibble(ht = seq(min(lab$ht, na.rm = TRUE), max(lab$ht, na.rm = TRUE), length.out = 100)))  # Add HT values back

#summarize lab dataset to count fish at each HT
lab_summary_ht <- lab %>%
  group_by(ht) %>%
  summarise(count = n(), .groups = "drop")

#merge summarized count data back for plotting
lab <- left_join(lab, lab_summary_ht, by = "ht")

#convert count into a factor for discrete color mapping
lab$count <- factor(lab$count, levels = 1:3)

#### LOGISTIC REGRESSION PLOT (HT) ------------------------------

#define color palette and point size mapping
custom_blues <- c("1" = "darkorange4", "2" = "darkorange", "3" = "gold1")
size_values <- c(2, 4, 6)

#order `lab` so that larger counts are plotted first, and smaller points appear on top
lab <- lab %>% arrange(desc(count))

#plot logistic regression prediction with confidence intervals
r2 <- ggplot(prediction_data_ht, aes(x = ht, y = fit)) +
  geom_ribbon(aes(ymin = min, ymax = max), alpha = 0.3, fill = "lightgrey") +  # Confidence Interval
  geom_line(size = 0.8, colour = "black") +  # Logistic regression line
  # Plot larger points first (from `lab` after sorting)
  geom_point(data = lab, aes(x = ht, y = m, size = count, color = count)) +  
  scale_size_manual(values = size_values, guide = guide_legend(title = "# of Fish")) +
  scale_color_manual(values = custom_blues, guide = guide_legend(title = "# of Fish")) +
  labs(x = "Handling Time (s)", y = "Survival Probability") +
  scale_y_continuous(limits = c(min(prediction_data_ht$min), max(prediction_data_ht$max)), 
                     breaks = seq(0, 1, by = 0.5)) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 16, colour = "black", margin = margin(t = 10)),
        axis.title.y = element_text(size = 16, colour = "black"),
        axis.text = element_text(size = 12, colour = "black"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))
r2



#### COMBINE TL & HT REGRESSION PLOTS --------------------------------------

#arrange plots side by side with equal width  
final_rplot <- r2 + r1 + plot_layout(ncol = 2, widths = c(1, 1)) +
  plot_annotation(tag_levels = 'a')
final_rplot

ggsave("log_regression_2panel.png", final_rplot, width = 12, height = 6, units = "in", dpi = 300)





####SURVIVAL ANALYSES---------------------------------

#HT SURVIVAL ANALYSIS
#bin handling time into two groups: ≤103 sec & >103 sec
d <- lab %>%
  mutate(int = case_when(is.na(int) ~ 40, TRUE ~ int)) %>%  # Assign 40 to survivors
  mutate(Handling_Time = case_when(ht > 103 ~ "> 103 sec", TRUE ~ "≤ 103 sec")) %>%  
  dplyr::select(ht, m, int, Handling_Time) %>%
  mutate(m = case_when(m == 0 ~ 2, TRUE ~ 1)) %>%  # Convert mortality to binary survival indicator
  left_join(Mort_date, by = c("int" = "Days"))

#convert mortality counts into binary survival data
d <- d %>%
  mutate(Control_Surv = case_when(Control == 0 ~ 2, TRUE ~ 1),
         Treatment_Surv = case_when(Treatment == 0 ~ 2, TRUE ~ 1))

#create survival objects for handling time groups, control, and treatment
ms_handling <- survfit(Surv(int, m) ~ Handling_Time, data = d)

#define custom colors for survival plot
custom_colours <- c("gold1", "darkorange")

#plot survival curves for handling time groups
p1 <- ggsurvplot(ms_handling, combine = TRUE, legend.title = "", palette = custom_colours, 
           legend.labs = c("Handling Time > 103 sec", "Handling Time ≤ 103 sec"), 
           censor = FALSE) +
  theme_survminer(font.x = c(16), font.y = c(16), font.legend = c(14)) +
  labs(x = "Time (Days)", y = "Survival Probability")
p1

ggsave("survival_analysis_handling_time.png", last_plot(), width = 8, height = 6, 
       units = "in", dpi = 300)

#Cox Proportional Hazard Model for HT Groups
mcph_ht <- coxph(Surv(int, m) ~ Handling_Time, data = d) 
summary(mcph_ht)


#TL SURVIVAL ANALYSIS
#bin total length into two groups: ≤ median & > median
tl_median <- median(lab$tl, na.rm = TRUE)  # Calculate median total length (= 32cm)

#count the number of fish above and below the median 
lab %>%
  mutate(TL_Group = case_when(tl > tl_median ~ "Above Median",
                              tl <= tl_median ~ "Below Median")) %>%
  count(TL_Group)
#Above = 20
#Below = 21

d <- lab %>%
  mutate(int = case_when(is.na(int) ~ 40, TRUE ~ int)) %>%  # Assign 40 to survivors
  mutate(TL_Group = case_when(tl > tl_median ~ paste("> ", tl_median, "cm", sep = ""),
                              TRUE ~ paste("≤ ", tl_median, "cm", sep = ""))) %>%
  dplyr::select(tl, m, int, TL_Group) %>%
  mutate(m = case_when(m == 0 ~ 2, TRUE ~ 1)) %>%  # Convert mortality to binary survival indicator
  left_join(Mort_date, by = c("int" = "Days"))

#convert mortality counts into binary survival data
d <- d %>%
  mutate(Control_Surv = case_when(Control == 0 ~ 2, TRUE ~ 1),
         Treatment_Surv = case_when(Treatment == 0 ~ 2, TRUE ~ 1))

#create survival objects for total length groups, control, and treatment
ms_tl <- survfit(Surv(int, m) ~ TL_Group, data = d)

#define custom colors for survival plot
custom_colours1 <- c("navyblue","lightskyblue")

#plot survival curves for total length groups
p2 <- ggsurvplot(ms_tl, combine = TRUE, legend.title = "", palette = custom_colours1, 
           legend.labs = c(paste("Total Length > ", tl_median, "cm", sep = ""),
                           paste("Total Length ≤ ", tl_median, "cm", sep = "")), 
           censor = FALSE) +
  theme_survminer(font.x = c(16), font.y = c(16), font.legend = c(14)) +
  labs(x = "Time (Days)", y = "Surival Probability")
p2

ggsave("survival_analysis_total_lenth.png", last_plot(), width = 8, height = 6, 
       units = "in", dpi = 300)

#Cox Proportional Hazard Model for Handling Time Groups
mcph_tl <- coxph(Surv(int, m) ~ TL_Group, data = d) 
summary(mcph_tl)


####COMBINE TL & HT SURVIVAL PLOTS --------------------------------------
#extract ggplot objects  
plot1 <- p1$plot +  
  theme(legend.position = "top",  
        plot.margin = margin(10, 10, 10, 10))  

plot2 <- p2$plot +  
  theme(legend.position = "top",  
        plot.margin = margin(10, 10, 10, 10))  

#Arrange plots side by side with equal width  
final_plot <- plot1 + plot2 + plot_layout(ncol = 2, widths = c(1, 1)) +
  plot_annotation(tag_levels = 'a')
final_plot

#Save the combined plot
ggsave("survival_analysis_2panel.png", final_plot, width = 12, height = 6, units = "in", dpi = 300)




####MEDIAN SIZE OF SURVIVORS --------------------------------------

#filter for survivors 
survivors <- lab %>% filter(m == 1)

#median TL of survivors
median_tl <- median(survivors$tl, na.rm = TRUE)
#range of all TL values
range_tl <- range(lab$tl, na.rm = TRUE)

#bootstrap function to get median (used here because of the small and skewed sample)
set.seed(123)  # for reproducibility
boot_median <- boot::boot(survivors$tl, function(data, i) median(data[i], na.rm = TRUE), R = 1000)

#get 95% CI from bootstrap
ci_median <- boot::boot.ci(boot_median, type = "perc")$percent[4:5]

cat("Median TL of survivors:", round(median_tl, 1), "cm\n")
cat("95% CI for median TL:", round(ci_median[1], 1), "–", round(ci_median[2], 1), "cm\n")
cat("Overall TL range in dataset:", range_tl[1], "–", range_tl[2], "cm\n")

#same thing for handling time 
#median HT of survivors
median_ht <- median(survivors$ht, na.rm = TRUE)
range_ht <- range(lab$ht, na.rm = TRUE)

set.seed(123)
boot_median_ht <- boot::boot(survivors$ht, function(data, i) median(data[i], na.rm = TRUE), R = 1000)

ci_median <- boot::boot.ci(boot_median_ht, type = "perc")$percent[4:5]

cat("Median HT of survivors:", round(median_ht, 1), "seconds\n")
cat("95% CI for median HT:", round(ci_median[1], 1), "–", round(ci_median[2], 1), "seconds\n")
cat("Overall HT range in dataset:", range_ht[1], "–", range_ht[2], "seconds\n")





####WOUND HEALING ANALYSES --------------------------------------
#ingest wound healing data

wound <- readRDS("wound.rds") %>% 
  select(2, 4) %>% 
  rename(id = Tag_num) %>% 
  mutate(id = as.character(id)) 

#join wound data into lab
lab <- lab %>% mutate(id = as.character(id)) %>% 
  left_join(wound, by = "id")

#plot
#set the monitoring period start date
monitoring_start <- as.Date("2023-11-11")

#filter only fish that died (m == 0), and calculate days until death
lab_dead <- lab %>%
  filter(m == 0) %>%
  mutate(
    redness = factor(Wound_condition,
                     levels = c(1, 2, 3),
                     labels = c("a", "b", "c")),
    days_until_death = as.numeric(as.Date(Date) - monitoring_start)
  ) %>%
  select(redness, tl, ht, days_until_death)

#pivot to long format for plotting
lab_long_dead <- lab_dead %>%
  pivot_longer(cols = c(tl, ht, days_until_death),
               names_to = "variable",
               values_to = "value") %>%
  mutate(variable = recode(variable,
                           tl = "Total Length (cm)",
                           ht = "Handling Time (sec)",
                           days_until_death = "Days Until Death"))

#boxplots of each variable by redness category (for dead fish)
ggplot(lab_long_dead, aes(x = redness, y = value)) +
  geom_boxplot(outlier.shape = NA, fill = "grey85", color = "black") +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7) +
  facet_wrap(~variable, scales = "free_y") +
  labs(x = "Wound Condition", y = NULL) +
  theme_classic() +
  theme(strip.text = element_text(size = 13),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))

ggsave("wound_condition_boxplots.png", last_plot(), width = 10, height = 6, 
       units = "in", dpi = 300)
