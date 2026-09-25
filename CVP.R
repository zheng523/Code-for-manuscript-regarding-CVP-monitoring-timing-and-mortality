# Load packages
library(DBI)
library(RPostgres)
library(tidyr)
library(dplyr)
library(ggplot2)
library(gtsummary)
library(rms)
library(missForest)
library(patchwork)
library(lubridate)
library(survey)

# Remove all objects except icustays
rm(list = ls()[ls() != "icustays"])

#### Retrieve data ####
# Establish a connection
con <- dbConnect(
  drv = RPostgres::Postgres(),  # specify the PostgreSQL driver
  host = "localhost",           # host address (localhost for local; remote needs IP/domain)
  port = 5432,                  # PostgreSQL default port 
  dbname = "mimiciv",          # database name (shown as mimiciv3 in the figure)
  user = "postgres",           # database user name
  password = "postgres"         # database password
)
# Verify whether the connection was successful
dbIsValid(con)
# Read the stay4 table
icustays <- dbReadTable(con, Id(schema = "public", table = "stay4"))
# Close the connection
dbDisconnect(con)
rm(con)

#### Variable processing ####
# Check whether survival days have negative values: none, because they were already removed during data extraction
range(icustays$survival_day)
# Remove records with negative survival days
stay <- icustays %>%
  filter(survival_day >= 0)
# Remove records with empty disease_group
stay <- stay %>%
  drop_na(disease_group) 
# Convert certain variables to factors; note that survival_status must not be converted to a factor. Set the reference group for some variables.
str(stay)
stay <- stay %>%
  mutate(across(c(thirty_day_death, one_year_death, hospital_expire_flag, gender, insurance, renal_disease, severe_liver_disease, metastatic_cancer, aids, vasoactive, mechanical_ventilation, heart_failure, sepsis, disease_group), factor)) %>% 
  mutate(
    insurance = relevel(insurance, ref = "Medicare"),
    disease_group = relevel(disease_group, ref = "sepsis")
  )
# Keep only HFsepsis patients and remove some variables
hfsepsis_raw <- stay %>%
  filter(disease_group == "HFsepsis") %>%
  select(-c(heart_failure, sepsis, sepsis_time, disease_group)) %>%
  select(-c(hadm_id, outtime, one_year_death)) 
# Remove duplicate ICU records for the same patient
hfsepsis_raw <- hfsepsis_raw %>%
  group_by(subject_id) %>%
  slice_min(intime, n = 1) %>%
  ungroup()
# Handle the issue that the insurance variable No charge is treated as 0
table(hfsepsis_raw$insurance)
hfsepsis_raw$insurance <- droplevels(hfsepsis_raw$insurance)
table(hfsepsis_raw$insurance)

# Convert the lvef variable to numeric
hfsepsis_raw <- hfsepsis_raw %>% 
  mutate(lvef = as.numeric(lvef)) 
range(hfsepsis_raw$lvef, na.rm=TRUE)
# Finally decided to convert records with lvef<5% and lvef>95% to missing values
hfsepsis_raw <- hfsepsis_raw %>%
  mutate(lvef = if_else(lvef < 5 | lvef > 95, NA, lvef))
# Adjust the order of variables
hfsepsis <- hfsepsis_raw %>%
  relocate(lvef, .after = bilirubin)
# Patient count statistics
nrow(hfsepsis)
summary(hfsepsis)
sum(!is.na(hfsepsis$first_cvp_interval))

# Define the names of all covariates
name_cov_all <- hfsepsis %>%
  select(age:mechanical_ventilation) %>%
  names()
# View the covariates
name_cov_all

#### Missing value processing ####
# Remove covariates with a missing proportion greater than 30%
# Calculate the missing rate
cov_missrate <-
  hfsepsis %>%
  select(age:mechanical_ventilation) %>%
  summarise(across(everything(), ~ mean(is.na(.)))) %>%
  pivot_longer(
    everything(),
    names_to = "variable",
    values_to = "missing_rate"
  )
# Calculate the number of missing cases
cov_missing_n <-
  hfsepsis %>%
  select(age:mechanical_ventilation) %>%
  summarise(across(everything(), ~ sum(is.na(.)))) %>%
  pivot_longer(
    everything(),
    names_to = "variable",
    values_to = "missing_n"
  )
# Merge the two results
cov_missrate <- cov_missing_n %>%
  left_join(cov_missrate, by = "variable") %>%
  arrange(desc(missing_rate))
View(cov_missrate)
# Export the missing rate table
write.csv(cov_missrate, file = "协变量缺失率_用首日值.csv")

# Names of variables with missing rate greater than 30%
name_cov_misshigh <- cov_missrate %>% 
  filter(missing_rate > 0.3) %>%
  pull(variable)
name_cov_misshigh # 4 variables in total have high missing rates
# Remove covariates with high missing rates
hfsepsis <- hfsepsis %>%
  select(-any_of(name_cov_misshigh))

# Check the current missing situation of covariates
cov_missrate_new <-
  hfsepsis %>% select(age:mechanical_ventilation) %>%
  summarise(across(everything(), ~mean(is.na(.)))) %>% 
  pivot_longer(
    everything(), 
    names_to = "variable", 
    values_to = "missing_rate" 
  ) %>%
  arrange(desc(missing_rate)) 
View(cov_missrate_new) 

# Set stay_id as the rownames for convenient merging later
hfsepsis <- as.data.frame(hfsepsis)
rownames(hfsepsis) <- hfsepsis$stay_id
# The dataset to be imputed, hfsepsis_cov, only includes covariates
hfsepsis_cov <- hfsepsis %>%
  select(age:mechanical_ventilation)
# Covariate names remaining after removing covariates with missing rate above 30%
name_cov <- names(hfsepsis_cov)

##### Random forest imputation #####
# Use random forest to impute covariates with missing proportion not exceeding 30%
library(missForest)
library(doParallel)
# Check the current CPU core count: 16 cores (actually 8 cores with 16 threads)
detectCores()
# Register a parallel cluster
cl <- makeCluster(12) # use 12 cores
registerDoParallel(cl)
# Start imputation
set.seed(123)
hfsepsis_cov_imputed <- missForest(
  xmis = hfsepsis_cov,
  maxiter = 15,
  ntree = 500,
  verbose = TRUE,
  parallelize = "forests"
)
# Close the parallel cluster
stopCluster(cl)
# Extract the imputed data
hfsepsis_cov_imputed <- hfsepsis_cov_imputed$ximp

# Restore the rownames back into stay_id
hfsepsis_cov_imputed$stay_id <- as.integer(rownames(hfsepsis_cov_imputed))
# Merge other variables by stay_id, then remove the stay_id variable after merging
hfsepsis_imputed <- hfsepsis_cov_imputed %>%
  full_join(
    (hfsepsis %>% 
       select(stay_id, los, survival_status, survival_day, thirty_day_death, hospital_expire_flag, first_cvp, first_cvp_interval )),
    by = "stay_id") %>%
  select(-stay_id)

#### Subgroup RCS analysis: the subgroup with CVP monitoring within 14 days ####
# Filter the subgroup hfsepsiscvp as those with CVP monitoring records and CVP monitoring within 14 days
hfsepsiscvp <- hfsepsis_imputed %>%
  filter(first_cvp_interval <= 336)

# Check the missing situation
library(visdat)
vis_dat(hfsepsis_imputed)
# vis_miss(hfsepsis_imputed)
vis_dat(hfsepsiscvp)

# Fit the RCS model
library(rms)
library(ggplot2)
# Define the reference value ranges for all variables in the dataset
dd <- datadist(hfsepsiscvp) 
# View the reference limits of each variable in different dimensions
dd$limits
options(datadist='dd') 
# Which model has the smallest AIC when the number of knots is 3:7
rcs_fit3 <- lrm(as.formula(
  paste("thirty_day_death ~ rcs(first_cvp_interval, 3) + ", paste(name_cov, collapse = " + "))),
  data = hfsepsiscvp)
rcs_fit4 <- lrm(as.formula(
  paste("thirty_day_death ~ rcs(first_cvp_interval, 4) + ", paste(name_cov, collapse = " + "))),
  data = hfsepsiscvp)
rcs_fit5 <- lrm(as.formula(
  paste("thirty_day_death ~ rcs(first_cvp_interval, 5) + ", paste(name_cov, collapse = " + "))),
  data = hfsepsiscvp)
rcs_fit6 <- lrm(as.formula(
  paste("thirty_day_death ~ rcs(first_cvp_interval, 6) + ", paste(name_cov, collapse = " + "))),
  data = hfsepsiscvp)
rcs_fit7 <- lrm(as.formula(
  paste("thirty_day_death ~ rcs(first_cvp_interval, 7) + ", paste(name_cov, collapse = " + "))),
  data = hfsepsiscvp)
AIC(rcs_fit3)
AIC(rcs_fit4)
AIC(rcs_fit5)
AIC(rcs_fit6)
AIC(rcs_fit7)

# Considering the AIC and model complexity under different numbers of knots, finally decided to choose 5 knots
rcs_fit <- lrm(as.formula(
  paste("thirty_day_death ~ rcs(first_cvp_interval, 5) + ", paste(name_cov, collapse = " + "))),
  data = hfsepsiscvp
)
# View the model
rcs_fit
# Non-linearity test
anova(rcs_fit)
# View the positions of the 5 knots of first_cvp_interval in the RCS model
rcs_fit$Design$parms$first_cvp_interval
# Use the model to generate predicted values and confidence intervals
OR <- rms::Predict(rcs_fit, 
              first_cvp_interval=seq(from=0, to=336, by=0.01),
              type = "predictions", 
              fun = exp, 
              ref.zero = TRUE) 
# first_cvp_interval at the lowest OR
lowest_x <- OR$first_cvp_interval[which.min(OR$yhat)]
lowest_y <- min(OR$yhat)
lowest_point <- data.frame(
  first_cvp_interval = lowest_x,
  OR = lowest_y
)
lowest_point

# Set the reference point to the nadir
dd$limits["Adjust to", "first_cvp_interval"] <- lowest_x
options(datadist='dd') 
# Refit the final model under the new reference point
rcs_fit <- lrm(as.formula(
  paste("thirty_day_death ~ rcs(first_cvp_interval, 5) + ", paste(name_cov, collapse = " + "))),
  data = hfsepsiscvp
)
anova(rcs_fit)
# Use the model to generate predicted values and confidence intervals
OR <- rms::Predict(rcs_fit, 
              first_cvp_interval=seq(from=0, to=336, by=0.01),
              type = "predictions", 
              fun = exp, 
              ref.zero = TRUE) 
# Redefine the nadir
lowest_x <- OR$first_cvp_interval[which.min(OR$yhat)]
lowest_y <- min(OR$yhat)
lowest_point <- data.frame(
  first_cvp_interval = lowest_x,
  OR = lowest_y
)
lowest_point

##### Find inflection points #####
# Step 1: generate a dense sequence + calculate the OR values
x_seq <- seq(0, 336, by = 0.01)  
pred_or <- Predict(rcs_fit, first_cvp_interval = x_seq, fun = exp, ref.zero = TRUE)
pred_df <- data.frame(
  first_cvp_interval = x_seq,
  OR = pred_or$yhat
)
# Step 2: use lead/lag to compute the first derivative
pred_df <- pred_df %>%
  mutate(
    # take the OR and x of the next row
    OR_next = lead(OR, 1),
    x_next = lead(first_cvp_interval, 1),
    # take the OR and x of the previous row
    OR_prev = lag(OR, 1),
    x_prev = lag(first_cvp_interval, 1),
    # first derivative (central difference method, using lead/lag instead of direct indexing)
    dOR_dx = case_when(
      # first row: forward difference (OR_next - OR)/(x_next - x)
      row_number() == 1 ~ (OR_next - OR) / (x_next - first_cvp_interval),
      # last row: backward difference (OR - OR_prev)/(x - x_prev)
      row_number() == n() ~ (OR - OR_prev) / (first_cvp_interval - x_prev),
      # middle rows: central difference (OR_next - OR_prev)/(x_next - x_prev)
      TRUE ~ (OR_next - OR_prev) / (x_next - x_prev)
    )
  ) %>%
  # remove the temporary columns (keep the data frame tidy)
  select(-OR_next, -x_next, -OR_prev, -x_prev)
# Step 3: use lead/lag to compute the second derivative
pred_df <- pred_df %>%
  mutate(
    # take the first derivative of the next and previous rows
    dOR_dx_next = lead(dOR_dx, 1),
    dOR_dx_prev = lag(dOR_dx, 1),
    # take the x of the next and previous rows
    x_next = lead(first_cvp_interval, 1),
    x_prev = lag(first_cvp_interval, 1),
    # second derivative
    d2OR_dx2 = case_when(
      row_number() == 1 ~ (dOR_dx_next - dOR_dx) / (x_next - first_cvp_interval),
      row_number() == n() ~ (dOR_dx - dOR_dx_prev) / (first_cvp_interval - x_prev),
      TRUE ~ (dOR_dx_next - dOR_dx_prev) / (x_next - x_prev)
    )
  ) %>%
  select(-dOR_dx_next, -dOR_dx_prev, -x_next, -x_prev)
# Step 4: detect inflection points
pred_df <- pred_df %>%
  mutate(
    d2_sign = sign(d2OR_dx2),
    sign_change = abs(d2_sign - lag(d2_sign, default = first(d2_sign)))
  )
# Filter the inflection points
inflection_points <- pred_df %>%
  filter(
    sign_change != 0, 
    !is.na(d2OR_dx2), 
    !is.infinite(d2OR_dx2),
    first_cvp_interval >= 0
  ) 
inflection_points


##### RCS plotting #####
library(ggplot2)
library(devEMF)

# General parameters (consistent with before, no modification needed)
bin_width <- 12  # bin width 6h
y_main_max <- 16 # maximum tick of OR on the main Y-axis
x_full_min <- 0
x_full_max <- 336

# First temporarily plot a histogram to extract the maximum density automatically computed by ggplot (for dual Y-axis scaling)
temp_hist <- ggplot(hfsepsiscvp, aes(x = first_cvp_interval)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = bin_width, boundary = 0)
max_density <- max(ggplot_build(temp_hist)$data[[1]]$density, na.rm = TRUE) # extract the internal maximum density
scale_factor <- y_main_max / max_density # calculate the scale factor (main Y-axis OR / maximum density)

# Figure 1 shows the full 0-14 day curve
rcs1 <- 
  ggplot() +
  geom_histogram(
    data = hfsepsiscvp, aes(x = first_cvp_interval, y = after_stat(density) * scale_factor),
    binwidth = bin_width, boundary = 0, # align the bins with the x-axis
    fill = "gray80", color = "gray60", linewidth = 0.2, alpha = 0.6, na.rm = TRUE) +
  geom_ribbon(data=OR, aes(x = first_cvp_interval, ymin = lower, ymax = upper), fill = '#B30000', alpha = 0.2)+
  geom_line(data = OR, aes(x = first_cvp_interval, y = yhat), linetype = 'solid', linewidth = 0.7, alpha = 1, color = '#B30000')+
  geom_hline(yintercept=1, linetype=2, color="grey") +
  scale_x_continuous(name = "Time to first CVP measurement (h)", breaks = seq(0, 336, by = 48)) +
  scale_y_continuous(name = "OR (95% CI)", breaks = seq(0, 16, by = 1),
                     sec.axis = sec_axis(
                       transform = ~ . / scale_factor, 
                       name = "Probability Density",
                       breaks = seq(0, round(max_density, 4) + 0.01, by = 0.01))) +
  geom_point(data = lowest_point, aes(x = first_cvp_interval, y = OR), color = "#0052CC", size = 1.5, shape = 16) +
  geom_text(aes(x = lowest_x, y = lowest_y-0.3, label = "Nadir"), hjust = 0.2, vjust = 1, size = 3, color = "#0052CC") +
  geom_point(data = lowest_point, aes(x = inflection_points[1,"first_cvp_interval"], y = inflection_points[1,"OR"]), color = "#6B00B3", size = 1.5, shape = 16) +
  geom_text(aes(x = inflection_points[1,"first_cvp_interval"]+7, y = inflection_points[1,"OR"], label = "Inflection point 1"), hjust = 0, vjust = 0, size = 3, color = "#6B00B3") +
  geom_point(data = lowest_point, aes(x = inflection_points[2,"first_cvp_interval"], y = inflection_points[2,"OR"]), color = "#6B00B3", size = 1.5, shape = 16) +
  geom_text(aes(x = inflection_points[2,"first_cvp_interval"], y = inflection_points[2,"OR"]-0.45, label = "Inflection point 2 (86.60 h)"), hjust = 0.1, vjust = 1, size = 3, color = "#6B00B3") +
  ggtext::geom_richtext(aes(x=48, y=15), label = "Overall <i>p</i> < 0.0001<br>Nonlinearity <i>p</i> < 0.0001", hjust=0, size=3.5, fill=NA, label.colour=NA) + 
  theme_bw() +
  theme(
    axis.line=element_line(),
    panel.grid=element_blank(),
    panel.border=element_blank()
  ) 

# Figure 2 shows the 0-12 hour curve
rcs2 <- ggplot() + 
  coord_cartesian(xlim = c(0, 12), ylim = c(0, 3)) +
  geom_ribbon(data=OR, aes(x = first_cvp_interval, ymin = lower, ymax = upper), fill = '#B30000', alpha = 0.2)+
  geom_line(data = OR, aes(x = first_cvp_interval, y = yhat), linetype = 'solid', linewidth = 0.7, alpha = 1, color = '#B30000')+
  geom_hline(yintercept=1, linetype=2, color="grey") +
  scale_x_continuous(name = "Time to first CVP measurement (h)", breaks = seq(0, 12, by = 2)) +
  scale_y_continuous(name = "OR (95% CI)", breaks = seq(0, 3, by = 0.5)) +
  geom_point(data = lowest_point, aes(x = first_cvp_interval, y = OR), color = "#0052CC", size = 1.5, shape = 16) +
  geom_text(aes(x = lowest_x, y = lowest_y-0.15, label = "Nadir (2.89 h)"), hjust = 0.5, vjust = 0.5, size = 3, color = "#0052CC") +
  geom_point(data = lowest_point, aes(x = inflection_points[1,"first_cvp_interval"], y = inflection_points[1,"OR"]), color = "#6B00B3", size = 1.5, shape = 16) +
  geom_text(aes(x = inflection_points[1,"first_cvp_interval"]+1.2, y = inflection_points[1,"OR"]-0.1, label = "Inflection point 1 (5.21 h)"), hjust = 0.3, vjust = 0.5, size = 3, color = "#6B00B3") +
  ggtext::geom_richtext(aes(x=1.6, y=2.54), label = "Overall <i>p</i> < 0.0001<br>Nonlinearity <i>p</i> < 0.0001", hjust=0, size=3.5, fill=NA, label.colour=NA) + 
  theme_bw() +
  theme(
    axis.line=element_line(),
    panel.grid=element_blank(),
    panel.border=element_blank()
  ) 

# Display side by side
library(patchwork)
rcs_plot <- rcs1 + rcs2 + 
  plot_annotation(tag_levels = 'A') +   # add labels
  theme(
    plot.tag = element_text(size = 12, face = 'bold', color = 'black'),   # 12pt bold black font
    plot.tag.position = c(0.01, 0.99)    # label position in the top-left corner
  )
ggsave(filename = "Figure 3 - RCS.emf", plot = rcs_plot, device = devEMF::emf, width = 4, height = 2, dpi = 300, scale = 2.25)
rcs_plot # manually export 800*400
ggsave(
  filename = "Figure 3 - RCS.eps",
  plot = rcs_plot,
  device = cairo_ps, 
  width = 8,
  height = 4,
  dpi = 300 
)

#### Grouping ####
# cvp_group1: divided into 2 groups: with or without CVP monitoring
hfsepsis_imputed <- hfsepsis_imputed %>%
  mutate(cvp_group1 = if_else(
    is.na(first_cvp_interval), 0, 1), 
  cvp_group1 = factor(
    cvp_group1, 
    levels = c(0, 1))
  )
table(hfsepsis_imputed$cvp_group1)
# cvp_group2: divided into 2 groups: early CVP monitoring, non-early CVP monitoring
hfsepsis_imputed <- hfsepsis_imputed %>%
  mutate(cvp_group2 = case_when(
    is.na(first_cvp_interval) ~ 0,
    first_cvp_interval<=6 ~ 1,
    first_cvp_interval>6 ~ 0
  ),
  cvp_group2 = factor(
    cvp_group2, 
    levels = c(0, 1))
  )
table(hfsepsis_imputed$cvp_group2)
# cvp_group4: divided into 4 groups: none, early, intermediate, late
hfsepsis_imputed <- hfsepsis_imputed %>%
  mutate(cvp_group4 = case_when(
    is.na(first_cvp_interval) ~ 0,
    first_cvp_interval<=6 ~ 1,
    first_cvp_interval>6 & first_cvp_interval<=86 ~ 2,
    first_cvp_interval>86 ~ 3
  ),
  cvp_group4 = factor(
    cvp_group4, 
    levels = c(0, 1, 2, 3))
  )
table(hfsepsis_imputed$cvp_group4)


#### Group 1 ####
library(gtsummary)
##### Crude model #####
crude_fit1 <- glm(
  thirty_day_death ~ cvp_group1,
  data = hfsepsis_imputed,
  family = binomial()
)
tbl_regression(
  crude_fit1, 
  exponentiate = TRUE,
  estimate_fun = function(x) style_number(x, digits = 3),  # OR and CI rounded to 3 decimal places
  pvalue_fun = function(x) style_pvalue(x, digits = 3)  # p value rounded to 3 decimal places
)

##### Multivariable logistic regression #####
glm_fit1 <- glm(
  formula = as.formula(paste("thirty_day_death ~ cvp_group1 +", paste(name_cov, collapse = " + "))),
  data = hfsepsis_imputed,
  family = binomial()
)
summary(glm_fit1)
tbl_regression(
  glm_fit1, 
  exponentiate = TRUE,
  estimate_fun = function(x) style_number(x, digits = 3),  # OR and CI rounded to 3 decimal places
  pvalue_fun = function(x) style_pvalue(x, digits = 3)  # p value rounded to 3 decimal places
)

##### IPTW #####
library(WeightIt)
library(cobalt)
library(survey)
library(pROC)
# Covariate balance check on the original data
bal.tab(as.formula(paste("cvp_group1 ~", paste(name_cov, collapse = " + "))),
        data = hfsepsis_imputed,
        estimand = "ATE",
        stats = c("m","v"),
        abs = TRUE,
        thresholds = c(m = 0.1, v = 2))
# Calculate the IPTW weights
iptw_obj1 <- weightit(
  formula = as.formula(paste("cvp_group1 ~", paste(name_cov, collapse = " + "))),
  data = hfsepsis_imputed,
  estimand = "ATE",
  method = "glm"
)
# Covariate balance check after weighting
bal.tab(iptw_obj1, 
        stats = c("m", "v"),
        abs = TRUE,
        thresholds = c(m = 0.1, v = 2))
# All variables satisfy balance
# Extract the weights and add them to the original data
hfsepsis_imputed$iptw_w1 <- weights(iptw_obj1)
# Extract the propensity scores and add them to the original data
hfsepsis_imputed$iptw_ps1 <- iptw_obj1$ps
# View the ranges of weights and PS
range(hfsepsis_imputed$iptw_w1)
range(hfsepsis_imputed$iptw_ps1)
# Define the weighted survey design (id=~1 when there is no clustering)
iptw_design1 <- svydesign(
  id = ~1,
  weights = ~iptw_w1,
  data = hfsepsis_imputed
)
# Weighted logistic regression
iptw_fit1 <- svyglm(
  formula = as.formula("thirty_day_death ~ cvp_group1"),
  design = iptw_design1,
  family = quasibinomial(link = "logit")
)
tbl_regression(
  iptw_fit1, 
  exponentiate = TRUE,
  estimate_fun = function(x) style_number(x, digits = 3), 
  pvalue_fun = function(x) style_pvalue(x, digits = 3)
)

#### Group 2 ####
library(gtsummary)
##### Crude model #####
crude_fit2 <- glm(
  thirty_day_death ~ cvp_group2,
  data = hfsepsis_imputed,
  family = binomial()
)
tbl_regression(
  crude_fit2, 
  exponentiate = TRUE,
  estimate_fun = function(x) style_number(x, digits = 3),  # OR and CI rounded to 3 decimal places
  pvalue_fun = function(x) style_pvalue(x, digits = 3)  # p value rounded to 3 decimal places
)

##### Multivariable logistic regression #####
glm_fit2 <- glm(
  formula = as.formula(paste("thirty_day_death ~ cvp_group2 +", paste(name_cov, collapse = " + "))),
  data = hfsepsis_imputed,
  family = binomial()
)
summary(glm_fit2)
tbl_regression(glm_fit2, exponentiate = TRUE)
tbl_regression(
  glm_fit2, 
  exponentiate = TRUE,
  estimate_fun = function(x) style_number(x, digits = 3),  # OR and CI rounded to 3 decimal places
  pvalue_fun = function(x) style_pvalue(x, digits = 3)  # p value rounded to 3 decimal places
)

##### LOS as outcome: multivariable lm regression #####
lm_fit2_los <- lm(
  formula = as.formula(paste("los ~ cvp_group2 +", paste(name_cov, collapse = " + "))),
  data = hfsepsis_imputed
)
summary(lm_fit2_los)
tbl_regression(
  lm_fit2_los, 
  estimate_fun = function(x) style_number(x, digits = 3),  # OR and CI rounded to 3 decimal places
  pvalue_fun = function(x) style_pvalue(x, digits = 3)  # p value rounded to 3 decimal places
)

##### In-hospital death as outcome: multivariable logistic regression #####
glm_fit2_ihm <- glm(
  formula = as.formula(paste("hospital_expire_flag ~ cvp_group2 +", paste(name_cov, collapse = " + "))),
  data = hfsepsis_imputed,
  family = binomial()
)
summary(glm_fit2_ihm)
tbl_regression(
  glm_fit2_ihm, 
  exponentiate = TRUE,
  estimate_fun = function(x) style_number(x, digits = 3),  # OR and CI rounded to 3 decimal places
  pvalue_fun = function(x) style_pvalue(x, digits = 3)  # p value rounded to 3 decimal places
)


##### IPTW #####
library(WeightIt)
library(cobalt)
library(survey)
library(pROC)
# Covariate balance check on the original data
bal.tab(as.formula(paste("cvp_group2 ~", paste(name_cov, collapse = " + "))),
        data = hfsepsis_imputed,
        estimand = "ATE",
        stats = c("m","v"),
        abs = TRUE,
        thresholds = c(m = 0.1, v = 2))
# Imbalanced
# Calculate the IPTW weights
iptw_obj2_glm <- weightit(
  formula = as.formula(paste("cvp_group2 ~", paste(name_cov, collapse = " + "))),
  data = hfsepsis_imputed,
  estimand = "ATE",
  method = "glm"
)
# Covariate balance check after weighting
bal.tab(iptw_obj2_glm, 
        stats = c("m", "v"),
        abs = TRUE,
        thresholds = c(m = 0.1, v = 2))
# Found that several variables still do not satisfy SMD<0.1, and spo2, lactate do not satisfy VR<2

# Use the gbm method to recalculate the weights (remember to set the random seed in advance to ensure reproducible results)
set.seed(123)
iptw_obj2 <- weightit(
  formula = as.formula(paste("cvp_group2 ~", paste(name_cov, collapse = " + "))),
  data = hfsepsis_imputed,
  estimand = "ATE",
  method = "gbm", 
  stop.method = "es.mean"
)
# View the weight distribution
plot(summary(iptw_obj2), xlim=c(0,10))
# Covariate balance check after weighting
bal.tab(iptw_obj2, 
        stats = c("m", "v"),
        abs = TRUE,
        thresholds = c(m = 0.1, v = 2))
# All balanced
love.plot(
  iptw_obj2, 
  stats = c("m", "v"), 
  thresholds = c(m = 0.1, v = 2),   # use 0.1 as the cutoff for SMD
  abs = TRUE,  # SMD shows the absolute value; the variance ratio divides the larger value by the smaller value
  var.order = "unadjusted",  # sort by the SMD of unadjusted variables from largest to smallest
  stars = "none"  # show the standardized mean difference (SMD) for both continuous and binary variables
)
# Extract the weights and add them to the original data
hfsepsis_imputed$iptw_w2 <- weights(iptw_obj2)
# Extract the propensity scores and add them to the original data
hfsepsis_imputed$iptw_ps2 <- iptw_obj2$ps
# View the ranges of weights and PS
range(hfsepsis_imputed$iptw_w2)
range(hfsepsis_imputed$iptw_ps2)
# Define the weighted survey design (id=~1 when there is no clustering)
iptw_design2 <- svydesign(
  id = ~1,
  weights = ~iptw_w2,
  data = hfsepsis_imputed
)
##### Weighted logistic regression #####
iptw_fit2 <- svyglm(
  formula = as.formula("thirty_day_death ~ cvp_group2"),
  design = iptw_design2,
  family = quasibinomial(link = "logit")
)
tbl_regression(
  iptw_fit2, exponentiate = TRUE,
  estimate_fun = function(x) style_number(x, digits = 3),
  pvalue_fun = function(x) style_pvalue(x, digits = 3) 
)
##### LOS as outcome: weighted lm regression #####
iptw_fit2_los <- svyglm(
  formula = as.formula("los ~ cvp_group2"), 
  design = iptw_design2,
  family = gaussian() 
)
summary(iptw_fit2_los)
tbl_regression(
  iptw_fit2_los, 
  estimate_fun = function(x) style_number(x, digits = 3), 
  pvalue_fun = function(x) style_pvalue(x, digits = 3) 
)
##### In-hospital death as outcome: weighted lm regression #####
iptw_fit2_ihm <- svyglm(
  formula = as.formula("hospital_expire_flag ~ cvp_group2"),
  design = iptw_design2,
  family = quasibinomial(link = "logit")
)
tbl_regression(
  iptw_fit2_ihm, 
  exponentiate = TRUE,
  estimate_fun = function(x) style_number(x, digits = 3), 
  pvalue_fun = function(x) style_pvalue(x, digits = 3) 
)

#### Output Table 1 ####
# Table 1 for the original data
library(tableone)
dput(names(hfsepsis_imputed))
tab1 <- CreateTableOne(
  strata = "cvp_group2",
  vars = c("thirty_day_death", "hospital_expire_flag", "los", name_cov), 
  data = hfsepsis_imputed)
tab1csv <- print(tab1, showAllLevels = TRUE, smd = TRUE, nonnormal = "los")
write.csv(tab1csv, "Table1-原始数据.csv")
# Weighted data
tab2 <- svyCreateTableOne(
  strata = "cvp_group2",
  vars = name_cov, 
  data = iptw_design2)
tab2csv <- print(tab2, showAllLevels = TRUE, smd = TRUE)
write.csv(tab2csv, "Table1-IPTW_GBM加权数据.csv")
svyvar(~spo2, subset(iptw_design2, cvp_group2 == 0))
sqrt(21.644)
# When finally creating Table 1, SMD and VR are not taken from the tableone package results, but from the earlier bal.tab() results

#### Forest plot ####
library(forestplot)
# Enter the data
fp1 <- read.table(
  header = TRUE,
  sep = "\t", 
  text = "
var	or	pvalue	or_mean	or_1	or_2
Model	OR (95% CI)	p-value			
Univariate-LR	0.540 (0.475 ~ 0.612) 	< 0.001	0.540	0.475	0.612
Multivariate-LR	0.576 (0.497 ~ 0.668)	< 0.001	0.576	0.497	0.668
IPTW	0.816 (0.688 ~ 0.967)	0.019	0.816	0.688	0.967
")
# Draw the forest plot
forestplot(
  labeltext=as.matrix(fp1[,1:3]),
  mean=fp1$or_mean,
  lower=fp1$or_1,
  upper=fp1$or_2,
  zero=1,  # position of the vertical line
  boxsize=0.2,  # box size
  graph.pos=3,  # the column in which the graph is positioned
  lwd.ci = 2,   # bolden the confidence interval lines and the box borders
  col = fpColors(line = "black", zero = "black"),   # set the colors
  ci.vertices = TRUE,  # show the vertices
  ci.vertices.height=0.1,  # vertex size
  xlim = c(0.4, 1),  # scale range
  xticks = c(0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1),  # position of the scale ticks
  txt_gp = fpTxtGp(
    ticks = gpar(cex = 0.8)  # size of the scale text
  )  
) # save as an EPS image at 800*250 resolution

#### Group 4 Landmark analysis ####
# Filter records with survival time > 86 hours
hfsepsis_imputed_2 <- hfsepsis_imputed %>%
  filter(survival_day > 3.59)
nrow(hfsepsis_imputed_2)

#### Group 4 logistic regression ####
library(gtsummary)
##### Crude model #####
crude_fit4 <- glm(
  thirty_day_death ~ cvp_group4,
  data = hfsepsis_imputed_2,
  family = binomial()
)
tbl_regression(
  crude_fit4, 
  exponentiate = TRUE,
  estimate_fun = function(x) style_number(x, digits = 3),  # OR and CI rounded to 3 decimal places
  pvalue_fun = function(x) style_pvalue(x, digits = 3)  # p value rounded to 3 decimal places
)

##### Multivariable logistic regression #####
glm_fit4 <- glm(
  formula = as.formula(paste("thirty_day_death ~ cvp_group4 +", paste(name_cov, collapse = " + "))),
  data = hfsepsis_imputed_2,
  family = binomial()
)
summary(glm_fit4)
tbl_regression(
  glm_fit4, 
  exponentiate = TRUE,
  estimate_fun = function(x) style_number(x, digits = 3),  # OR and CI rounded to 3 decimal places
  pvalue_fun = function(x) style_pvalue(x, digits = 3)  # p value rounded to 3 decimal places
)

#### IPTW for group 4 ####
library(WeightIt)
library(cobalt)
library(survey)
library(pROC)
# Covariate balance check on the original data
bal.tab(as.formula(paste("cvp_group4 ~", paste(name_cov, collapse = " + "))),
        data = hfsepsis_imputed_2,
        estimand = "ATE",
        stats = c("m","v"),
        abs = TRUE,
        thresholds = c(m = 0.1, v = 2))
# Calculate the IPTW weights
iptw_obj4_glm <- weightit(
  formula = as.formula(paste("cvp_group4 ~", paste(name_cov, collapse = " + "))),
  data = hfsepsis_imputed_2,
  estimand = "ATE",
  method = "glm"
)
# Covariate balance check after weighting
bal.tab(iptw_obj4_glm, 
        stats = c("m", "v"),
        abs = TRUE,
        thresholds = c(m = 0.1, v = 2))
# Found that multiple covariates are imbalanced

# Use the gbm method to recalculate the weights
set.seed(123)
iptw_obj4 <- weightit(
  formula = as.formula(paste("cvp_group4 ~", paste(name_cov, collapse = " + "))),
  data = hfsepsis_imputed_2,
  estimand = "ATE",
  method = "gbm", 
  stop.method = "es.mean"
)
# Covariate balance check after weighting; still not fully balanced
bal.tab(iptw_obj4, 
        stats = c("m", "v"),
        abs = TRUE,
        thresholds = c(m = 0.1, v = 2))
# Use the cbps method to calculate the weights
set.seed(123)
iptw_obj4 <- weightit(
  formula = as.formula(paste("cvp_group4 ~", paste(name_cov, collapse = " + "))),
  data = hfsepsis_imputed_2,
  estimand = "ATE",
  method = "cbps"
)
# Covariate balance check; basically balanced
bal.tab(iptw_obj4, 
        stats = c("m", "v"),
        abs = TRUE,
        thresholds = c(m = 0.1, v = 2))
summary(iptw_obj4)
# Extract the weights and add them to the original data
hfsepsis_imputed_2$iptw_w4 <- weights(iptw_obj4)
# View the ranges of weights and PS
range(hfsepsis_imputed_2$iptw_w4)
# Define the weighted survey design (id=~1 when there is no clustering)
iptw_design4 <- svydesign(
  id = ~1,
  weights = ~iptw_w4,
  data = hfsepsis_imputed_2
)
# Weighted logistic regression
iptw_fit4 <- svyglm(
  formula = as.formula("thirty_day_death ~ cvp_group4"),
  design = iptw_design4,
  family = quasibinomial(link = "logit")
)
tbl_regression(
  iptw_fit4, 
  exponentiate = TRUE,
  estimate_fun = function(x) style_number(x, digits = 3), 
  pvalue_fun = function(x) style_pvalue(x, digits = 3)
)

#### Forest plot ####
library(forestplot)
# Enter the data
fp2 <- read.table(
  header = TRUE,
  sep = "\t", 
  text = "
var	comparison	or	pvalue	or_mean	or_1	or_2
Model	Comparison	OR (95% CI)	p-value			
Multivariate-LR						
	Early vs. No	0.679 (0.572 ~ 0.804)	< 0.001	0.679	0.572	0.804
	Intermediate vs. No	1.444 (1.240 ~ 1.680)	< 0.001	1.444	1.240	1.680
	Late vs. No	3.491 (2.594 ~ 4.693)	< 0.001	3.491	2.594	4.693
IPTW						
	Early vs. No	0.881 (0.700 ~ 1.107)	0.277	0.881	0.700	1.107
	Intermediate vs. No	1.396 (1.194 ~ 1.632)	< 0.001	1.396	1.194	1.632
	Late vs. No	3.141 (2.263 ~ 4.359)	< 0.001	3.141	2.263	4.359

")
# Draw the forest plot
forestplot(
  labeltext=as.matrix(fp2[,1:4]),
  mean=fp2$or_mean,
  lower=fp2$or_1,
  upper=fp2$or_2,
  zero=1,  # position of the vertical line
  boxsize=0.2,  # box size
  graph.pos=4,  # the column in which the graph is positioned
  lwd.ci = 2,   # bolden the confidence interval lines and the box borders
  col = fpColors(line = "black", zero = "black"),   # set the colors
  ci.vertices = TRUE,  # show the vertices
  ci.vertices.height=0.1,  # vertex size
  xlim = c(0.5, 5),  # scale range
  xticks = c(0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5),  # position of the scale ticks
  txt_gp = fpTxtGp(
    ticks = gpar(cex = 0.8)  # size of the scale text
  ), 
  align = c("l", "l", "c", "c")   # add this line
) # save as EPS at 1200*500 resolution
