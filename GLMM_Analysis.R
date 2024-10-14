# --------------------------------------------------------------------------------------------
# Kajanka J Mathiaparanam
# 26th October 2020
# GLMM analysis to test the cumulative impacts of noise and light pollution on Coot abundance
# --------------------------------------------------------------------------------------------

# Assign current path
curdir <- getwd()

# Load packages
library(readxl)
library(tidyverse)
library(dplyr)
library(lme4)
library(lmerTest)
library(AICcmodavg)
library(car)
library(ggplot2)
library(stringr)


#### 1) Functions ####

# 1A) Mode function
Mode <- function(x, na.rm = TRUE) {
  if(na.rm){
    x = x[!is.na(x)]
  }
  ux <- unique(x)
  return(ux[which.max(tabulate(match(x, ux)))])
}

# 1B) Function to find waterbody ID that has an area closest to the mean waterbody ID
Rand_mode <- function(x, y) {
  return(x[which.min(abs(y - mean(y, na.rm = T)))])
}

# 1C) Functions for bootMer() and objects

# 1C-1) Return predicted values from bootstrap
Pred_func <- function(.) {
  predict(., newdata = Predictor_matrix, type = "response", re.form = NULL)
}

# 1C-2) Collapse bootstrap into median, 95% PI
Boot_sum <- function(merBoot) {
  return(
    data.frame(fit = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs = .5, na.rm = TRUE))),
               lwr = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs = .025, na.rm = TRUE))),
               upr = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs = .975, na.rm = TRUE)))
    )
  )
}


#### 2) Assign variables, filter thresholds and model attributes ####

# Read in the coot data
BLA_data <- read.csv(paste(curdir, "/BLA_Prepped_Data.csv", sep = ""))

# 2A) Assign species
Bird_Species <- "Eurasian Coot"

# 2B) Assign response, fixed and random predictors
Resp_var <- "Count"
Pred_var <- c("Poi_Noise_RR", "Poi_Lmean13_18", "Area_m2", "Poi_Temp", "Poi_Rain", "WB_GrVeg_Cat")
Pred_var_rename <- c("Noise", "Light", "Area_m2", "Temp", "Rain", "Veg")
Random_var <- "WBCat2_ID"

# 2C) Predictors to be used for interactions
Int_preds <- c("Noise", "Light")

# 2D) Assign family
Mod_family <- "poisson"

# 2E) Waterbody Area filter
Area_filt <- 10^6

# 2F) Noise Filter
Noi_thresh <- 30

# 2G) Make predictions with confidence intervals (bootstrapping)?
  # 1 - Yes
  # 0 - No
Predict_boot_YN <- 0

# 2H) Colinearity thresholds
Cor_thresh <- 0.7

# 2I) Cumulative plot matrix dimentions
Int_Pred_matrix_length <- 100


#### 3) Preparation for model ####

Model_df_cols <- c("Formula", "Family", "K", "AICc", "Delta_AICc", "AICcWt", "Intercept_Est",
                   "Intercept_Pr", "Intercept_Std_Err", "Noise_Est", "Noise_Pr", "Noise_Std_Err", "Light_Est", "Light_Pr",
                   "Light_Std_Err", "Area_Est", "Area_Pr", "Area_Std_Err", "Temp_Est", "Temp_Pr", "Temp_Std_Err", "Rain_Est",
                   "Rain_Pr", "Rain_Std_Err", "Veg1_Est", "Veg1_Pr", "Veg1_Std_Err", "Veg2_Est", "Veg2_Pr", "Veg2_Std_Err",
                   "Veg3_Est", "Veg3_Pr", "Veg3_Std_Err", "Veg4_Est", "Veg4_Pr", "Veg4_Std_Err", "Veg5_Est", "Veg5_Pr",
                   "Veg5_Std_Err", "Noise.Light_Est", "Noise.Light_Pr", "Noise.Light_Std_Err")
                   
# Empty df to save model outputs
Model_df <- data.frame(matrix(vector(), 0, length(Model_df_cols)))
names(Model_df) <- c(Model_df_cols)

# Create folder for model summary csv and plots
Plot_folder_name <- str_replace_all(paste("GLMM_", Sys.time(), sep = ""), ":", "-")
Ouput_path <- paste(curdir, "/Model_Outputs/", Plot_folder_name, sep = "")
dir.create(Ouput_path)


#### 4) Subsetting, cleaning and transformation ####

# Filter by species (and survey type), and select response and predictor columns
BLA_data_temp <- BLA_data %>%
filter(Common_Name %in% Bird_Species & Survey_Type %in% "Area") %>%
dplyr::select(all_of(c("Lat", "Lon", Resp_var, Pred_var, Random_var)))

# Remove all NA's
BLA_data_temp <- na.omit(BLA_data_temp)

# Rename columns to apply filters + plot etc. easily
names(BLA_data_temp)[4:(3 + length(Pred_var))] <- Pred_var_rename

# Noise Filter (remove surveys with less than 30dB)
BLA_data_temp <- BLA_data_temp %>% filter(Noise >= Noi_thresh)

# Apply Area filter
BLA_data_temp <- BLA_data_temp %>% filter(Area_m2 <= Area_filt)

# Convert Categorical variables to factors
Cat_var <- c("Veg", "WBCat2_ID")
for (cat_var_i in 1:length(Cat_var)) {
  BLA_data_temp[[Cat_var[cat_var_i]]] <- as.factor(BLA_data_temp[[Cat_var[cat_var_i]]])
}

# Save a copy of the data after filters - before predictor transformation for predictions and post filtering histograms
BLA_data_temp_preTrans <- BLA_data_temp

# Predictor transformation
Pred_trans <- Pred_var_rename[!(Pred_var_rename == "Veg")]
for (temp_col in 1:length(Pred_trans)) {
  BLA_data_temp[Pred_trans[temp_col]] <- (BLA_data_temp[Pred_trans[temp_col]] - mean(BLA_data_temp[[Pred_trans[temp_col]]], na.rm = T)) / sd(BLA_data_temp[[Pred_trans[temp_col]]], na.rm = T)
}


#### 5) Predictor selection based on AIC and colinearity ####

# 5A) Initial models - predictor selection based on AIC

# Predictors and intercept to be used to run initial models to choose predictors
Pred_var_sel <- c(1, Pred_var_rename)

# Add random factor
Ini_preds <- paste(Pred_var_sel, "(1|WBCat2_ID)", sep = "+")

# Pair with Response
Ini_var_comb <- expand.grid(Resp_var, Ini_preds)

# Create complete formulas
Ini_formula_vec <- sprintf("%s ~ %s", Ini_var_comb$Var1, Ini_var_comb$Var2)

# Run initial models
glm_initial <- lapply(Ini_formula_vec, function (f) {
  fitted_model <- glmer(f, data = BLA_data_temp, family = Mod_family)
  return (fitted_model)
})
names(glm_initial) <- Ini_formula_vec

# Generate AICc values
Ini_AIC <- aictab(glm_initial, Ini_formula_vec, second.ord = TRUE, sort = FALSE)

# Choose the predictors with AIC less than null model
Ini_AIC <- cbind(Pred_var_sel, Ini_AIC)
Pred_var_sel <- Pred_var_sel[Ini_AIC$AICc < Ini_AIC$AICc[1]] # Removes intercept anyway


# 5B) Colinearity filter

# Predictor correlation matrix
Cor_matrix <- abs(cor(apply(BLA_data_temp[Pred_var_rename], 2, as.numeric)))

# Remove predictors that correlate with other predictors (keep the one with lower AIC)
Pred_rename_aic <- Pred_var_rename[Pred_var_rename %in% Pred_var_sel]

if (length(Pred_rename_aic) > 1) {
  for (Cor_i in 1:(length(Pred_rename_aic) - 1)) {
    for (Cor_ii in (Cor_i + 1):length(Pred_rename_aic)) {
      if (Cor_matrix[Pred_rename_aic[Cor_i], Pred_rename_aic[Cor_ii]] > Cor_thresh) {
        Cor_pred_dual <- c(Pred_rename_aic[Cor_i], Pred_rename_aic[Cor_ii])
        cor_temp_df <- Ini_AIC %>% filter(Pred_var_sel %in% Cor_pred_dual)
        Pred_var_sel <- Pred_var_sel[!Pred_var_sel %in% cor_temp_df$Pred_var_sel[!cor_temp_df$AICc == min(cor_temp_df$AICc)]]
      }
    }
  }
}


#### 6) Formulas ####

# 6A) Generate predictor combinations

# All combinations of predictors
Pred_vars_comb <- unlist(sapply(seq_len(length(Pred_var_sel)),
                                function (i) {
                                  apply(combn(Pred_var_sel, i), 2, function (x) {paste(x, collapse = "+")})
                                }))

# Number of predictors to get combinations of = sequence of 1 to length of predictor vector
Pred_combn_i <- seq_len(length(Pred_var_sel))

# Column names for dataframe to store combinations of predictors
Pred_sel_col <- paste("Var_", Pred_combn_i, sep = "")

# Empty dataframe
Pred_sel_df <- data.frame(matrix(vector(), 0, length(Pred_var_sel)))
names(Pred_sel_df) <- Pred_sel_col

# Fill dataframe with all combinations of predictors
for (i in Pred_combn_i) {
  temp_row <- as.data.frame(t(combn(Pred_var_sel, i, simplify = T)))
  temp_col <- matrix(rep(NA, (length(Pred_var_sel) - i) * nrow(temp_row)), nrow(temp_row), length(Pred_var_sel) - i)
  temp_row_col <- cbind(temp_row, temp_col)
  names(temp_row_col) <- Pred_sel_col
  
  Pred_sel_df <- rbind(Pred_sel_df, temp_row_col)
}


# 6B) Generate interactions for each combination of predictors

# Column names for interactions df
Pred_int_col <- c(paste(Int_preds[1], Int_preds[2], sep = "_"))

# Empty dataframe
Pred_int_df <- data.frame(matrix(vector(), 0, sum(seq_len(length(Pred_int_col)))))

# Rename interactions df
names(Pred_int_df) <- Pred_int_col

# Fill dataframe with interaction terms
for (ii in 1:nrow(Pred_sel_df)) {
  temp_row <- Pred_sel_df[ii, ][!is.na(Pred_sel_df[ii, ])]
  Pred_int <- vector()
  for (i in 1:length(Int_preds)) {
    j_vec <- seq_len(length(Int_preds))
    j_vec <- j_vec[j_vec > i]
    for (j in j_vec) {
      if (any(temp_row == Int_preds[i]) & any(temp_row == Int_preds[j])) {
        Pred_int <- c(Pred_int, paste(Int_preds[i], Int_preds[j], sep = "*"))
      } else {
        Pred_int <- c(Pred_int, NA)
      }
    }
  }
  names(Pred_int) <- Pred_int_col
  Pred_int_df <- rbind(Pred_int_df, Pred_int)
  names(Pred_int_df) <- Pred_int_col
}


# 6C) Combine predictor combinations with corresponding interaction terms
Pred_var_int_comb <- vector()

for (Predvar_i in 1:length(Pred_vars_comb)) {
  
  Pred_int_df_temp <- Pred_int_df[Predvar_i, ][!is.na(Pred_int_df[Predvar_i, ])]
  
  for (Predint_i in 0:length(Pred_int_df_temp)) {
    if (Predint_i == 0) {
      Pred_comb_temp <- Pred_vars_comb[Predvar_i]
    } else {
      Pred_comb_temp <- paste(Pred_vars_comb[Predvar_i], Pred_int_df_temp[Predint_i], sep = "+")
    }
    Pred_var_int_comb <- c(Pred_var_int_comb, Pred_comb_temp)
  }
}


# 6D) Generate complete formulas

# Add the null model (intercept)
Pred_var_int_comb <- c(1, Pred_var_int_comb)

# Add random factor(s) to predictor combinations
Pred_var_int_comb <- paste(Pred_var_int_comb, "(1|WBCat2_ID)", sep = "+")

# Pair with Response
var_comb <- expand.grid(Resp_var, Pred_var_int_comb)

# Create complete formulas
formula_vec <- sprintf("%s ~ %s", var_comb$Var1, var_comb$Var2)


#### 7) Run Models ####

# Run all possible models with interactions
glm_mixed <- lapply(formula_vec, function (f) {
  fitted_model <- glmer(f, data = BLA_data_temp, family = Mod_family)
  return (fitted_model)
})
names(glm_mixed) <- formula_vec

# Generate AICc values
mod_aic <- aictab(glm_mixed, formula_vec, second.ord = TRUE, sort = FALSE)

# Extract coefficients + Generate df with summaries of all models
Pred_int_var <- c("(Intercept)", "Noise", "Light", "Area_m2", "Temp", "Rain", "Veg1", "Veg2", "Veg3", "Veg4", "Veg5",
                "Noise:Light")

for (mod_num in 1:length(glm_mixed)) {
  coef_model <- coef(summary(glm_mixed[[mod_num]]))
  est_p_vec <- data.frame(matrix(vector(), 1, 0))
  
  for (Pred_int_i in 1:length(Pred_int_var)) {
    if (any(row.names(coef_model) == Pred_int_var[Pred_int_i])) {
        est_p_vec <- cbind(est_p_vec,
                           coef_model[Pred_int_var[Pred_int_i], "Estimate"],
                           coef_model[Pred_int_var[Pred_int_i], "Pr(>|z|)"],
                           coef_model[Pred_int_var[Pred_int_i], "Std. Error"])
    } else {
      est_p_vec <- cbind(est_p_vec, as.numeric(NA), as.numeric(NA), as.numeric(NA))
    }
  }
  mod_vec <- cbind(names(glm_mixed)[mod_num], Mod_family,
                   mod_aic[mod_num, "K"], mod_aic[mod_num, "AICc"], mod_aic[mod_num, "Delta_AICc"],
                   mod_aic[mod_num, "AICcWt"], est_p_vec)
  mod_vec <- as.data.frame(mod_vec)
  names(mod_vec) <- c(Model_df_cols)
  Model_df <- rbind(Model_df, mod_vec)
}


#### 8) Save model predictions + Get limits for response scale for individual effects ####

# Extract the best model as per AIC
best_model_formula <- as.character(mod_aic$Modnames[mod_aic$AICc == min(mod_aic$AICc)])
best_model <- glm_mixed[[best_model_formula]]

# Extract the equivalent additive model - for 3D plot comparison with best model
best_model_formula_wo_int <- str_remove(best_model_formula, fixed("Noise*Light+"))
best_model_wo_int <- glm_mixed[[best_model_formula_wo_int]]

# 8A) Individual effects
for (i in 1:length(Pred_var_rename)) {

  # Predictor matrix, to which predictions should be made
  if (is.factor(BLA_data_temp[[Pred_var_rename[i]]])) {
    Predictor_matrix <- as.data.frame(unique(BLA_data_temp[Pred_var_rename[i]]))
  } else {
    Predictor_matrix <- as.data.frame(seq(min(BLA_data_temp[Pred_var_rename[i]]), max(BLA_data_temp[Pred_var_rename[i]]), length.out = 100))
  }
  names(Predictor_matrix)[1] <- Pred_var_rename[i]

  # Extract means and modes of the remaining predictors for the predictor matrix
  Pred_mean <- Pred_trans[!Pred_trans %in% Pred_var_rename[i]]
  Mean_matrix <- t(apply(BLA_data_temp[Pred_mean], 2, mean, na.rm = T))
  Pred_mode <- Pred_var_rename[Pred_var_rename %in% Cat_var]
  Pred_mode <- Pred_mode[!Pred_mode %in% Pred_var_rename[i]]
  Mode_matrix <- t(apply(BLA_data_temp[Pred_mode], 2, Mode, na.rm = T))
  Rand_matrix <- t(apply(BLA_data_temp[Random_var], 2, Rand_mode, y = BLA_data_temp$Area_m2))
  Predictor_matrix <- cbind(Predictor_matrix[Pred_var_rename[i]], Mean_matrix, Mode_matrix, Rand_matrix)
  
  # Predictions
  # Method 1 - lme4::predict.merMod (without std. errors)
  Pred_response <- predict(best_model, newdata = Predictor_matrix, type = "response", re.form = NA)
  # Method 2 - lmer::bootMer (with standard errors)
  if (Predict_boot_YN == 1) {
    Boot_out <- lme4::bootMer(best_model, Pred_func, nsim = 100, use.u = F, type = "parametric", re.form = NA)
    Boot_sum_df <- Boot_sum(Boot_out)
  }
  
  # Combine the predictor matrix and predicted responses
  if (Predict_boot_YN == 1) {
    Resp_Pred_Matrix <- cbind(Predictor_matrix, Pred_response, Boot_sum_df)
  } else {
    Resp_Pred_Matrix <- cbind(Predictor_matrix, Pred_response)
  }
  names(Resp_Pred_Matrix)[length(Predictor_matrix) + 1] <- Resp_var
  
  # Back transform predictor values to real values
  for (temp_col in 1:length(Pred_trans)) {
    Resp_Pred_Matrix[Pred_trans[temp_col]] <- (Resp_Pred_Matrix[Pred_trans[temp_col]] * sd(BLA_data_temp_preTrans[[Pred_trans[temp_col]]], na.rm = T)  + mean(BLA_data_temp_preTrans[[Pred_trans[temp_col]]], na.rm = T))
  }

  # Save the predictor-response matrix to plot predictions next
  assign(paste0("Resp_Pred_Matrix_", i), Resp_Pred_Matrix)
  
  # Save the predictor-response matrix as a csv
  dir.create(paste(Ouput_path, "/Pred_Resp_Matrices", sep = ""))
  write.csv(Resp_Pred_Matrix, file = paste(Ouput_path, "/Pred_Resp_Matrices/GLMM_Prediction_", Bird_Species, "_",
                                           Resp_var, "_", Pred_var_rename[i], ".csv", sep = ""), quote = F, row.names = F)
}


# 8B) Cumulative effects of noise and light

# Predictor matrix, to which predictions should be made
for (j in 1:2) {
  if (is.factor(BLA_data_temp[Int_preds[j]])) {
    Pred_col <- unique(BLA_data_temp[Int_preds[j]])
  } else {
    Pred_col <- seq(min(BLA_data_temp[Int_preds[j]]), max(BLA_data_temp[Int_preds[j]]), length.out = Int_Pred_matrix_length)
  }
  assign(paste0("Pred_col", j), Pred_col)
}
Int_Pred_matrix <- expand.grid(Pred_col1, Pred_col2)
names(Int_Pred_matrix) <- Int_preds

# Extract means and modes of the remaining predictors for the predictor matrix
Pred_mean <- Pred_trans[!Pred_trans %in% Int_preds]
Mean_matrix <- t(apply(BLA_data_temp[Pred_mean], 2, mean, na.rm = T))
Pred_mode <- Pred_var_rename[Pred_var_rename %in% Cat_var]
Pred_mode <- Pred_mode[!Pred_mode %in% Int_preds]
Mode_matrix <- t(apply(BLA_data_temp[Pred_mode], 2, Mode, na.rm = T))
Rand_matrix <- t(apply(BLA_data_temp[Random_var], 2, Rand_mode, y = BLA_data_temp$Area_m2))
Int_Pred_matrix <- cbind(Int_Pred_matrix, Mean_matrix, Mode_matrix, Rand_matrix)

# Predictions
# Best Model
Int_Pred_Response <- predict(best_model, newdata = Int_Pred_matrix, type = "response", re.form = NA)
# Additive equivalent of the interactive model
Int_Pred_Response_wo_Int <- predict(best_model_wo_int, newdata = Int_Pred_matrix, type = "response", re.form = NA)

# Round-up the predicted values
Int_Pred_Response <- round(Int_Pred_Response, 2)
Int_Pred_Response_wo_Int <- round(Int_Pred_Response_wo_Int, 2)

# Calculate difference between the predicted response and predicted response without the interaction term
Int_Pred_Response_diff <- Int_Pred_Response_wo_Int - Int_Pred_Response

# Combine the predictor matrix and predicted responses
Int_Resp_Pred_Matrix <- cbind(Int_Pred_matrix, Int_Pred_Response, Int_Pred_Response_wo_Int, Int_Pred_Response_diff)
Int_Resp_names <- c(Resp_var, str_c(Resp_var, "_wo_Int", sep = ""), str_c(Resp_var, "_diff", sep = ""))
names(Int_Resp_Pred_Matrix)[(length(Int_Pred_matrix) + 1) : (length(Int_Pred_matrix) + 3)] <- Int_Resp_names

# Back transform predictor values back to real values
for (temp_col in 1:length(Pred_trans)) {
  Int_Resp_Pred_Matrix[Pred_trans[temp_col]] <- (Int_Resp_Pred_Matrix[Pred_trans[temp_col]] * sd(BLA_data_temp_preTrans[[Pred_trans[temp_col]]], na.rm = T)  + mean(BLA_data_temp_preTrans[[Pred_trans[temp_col]]], na.rm = T))
}

# Select only the columns that are required for plotting (Noise, Light and Predictions Cols)
Int_Resp_Pred_Matrix <- Int_Resp_Pred_Matrix %>% select(Int_preds, Int_Resp_names)

# Save the predictor-response matrix as a csv
write.csv(Int_Resp_Pred_Matrix, file = paste(Ouput_path, "/Pred_Resp_Matrices/GLMM_Prediction_", Bird_Species, "_",
                                       Resp_var, "_Int_NL", ".csv", sep = ""), quote = F, row.names = F)


#### 9)  Save plots ####

# Save pdf
dir.create(paste(Ouput_path, "/Plots", sep = ""))
pdf(paste(Ouput_path, "/Plots/", Bird_Species, "_", Resp_var, ".pdf", sep = ""))

# Plot the predicted response curves
for (i in 1:length(Pred_var_rename)) {

  # Get the corresponding predictor-response matrix
  Resp_Pred_Matrix <- get(paste0("Resp_Pred_Matrix_", i))
  
  # Plot with confidence intervals using ggplot2
  if (Predict_boot_YN == 1) {
    print(
    ggplot(data = Resp_Pred_Matrix, aes(x = Resp_Pred_Matrix[[Pred_var_rename[i]]], y = fit, ymin = lwr, ymax = upr)) +
      geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "grey75") +
      geom_line(aes(y = fit), size = 1) +
      labs(x = Pred_var_rename[i], y = Resp_var)
    )
  }
  
  # Plot without confidence intervals using baseplot
  plot(Resp_Pred_Matrix[[Pred_var_rename[i]]], Resp_Pred_Matrix[[Resp_var]], cex = 0, xlab = Pred_var_rename[i], ylab = Resp_var)
  
  # Draw the trend for non-categorical
  if (!is.element(Pred_var_rename[i], Cat_var)) {
    Resp_spline <- smooth.spline(Resp_Pred_Matrix[[Pred_var_rename[i]]], Resp_Pred_Matrix[[Resp_var]])
    lines(Resp_spline)
  }
}

# Plot interactions (Noise and Light)
Int_plot_resp_label <- c("Count", "Count", "Difference in Count")
for (Int_plot_resp_i in 1:length(Int_Resp_names)) {
  Int_plot_resp <- Int_Resp_names[Int_plot_resp_i]
  print(
    Int_plot <- ggplot(data = Int_Resp_Pred_Matrix, aes(x = .data[[Int_preds[1]]], y = .data[[Int_preds[2]]], z = .data[[Int_plot_resp]]))
    + geom_tile(aes(fill = .data[[Int_plot_resp]]))
    + scale_fill_gradientn(colours = rainbow(7, start = 0, end = 0.4, alpha = 1, s = 0.75, v = 1), trans = "log", labels = scales::number_format(accuracy = 0.01))
    + labs(x = Int_preds[1], y = Int_preds[2], fill = Int_plot_resp_label[Int_plot_resp_i], title = Bird_Species)
    + theme_light()
    + theme(title = element_text(size = 9), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10))
    + theme(legend.position = "bottom", legend.key.height = unit(2, "mm"), legend.key.width = unit(10, "mm"))
  )
}

dev.off()


#### 10) Save dataframes with model summaries ####

# Rename model results df columns and model options df columns
names(Model_df) <- c(Model_df_cols)

# write dataframe
write.csv(Model_df, file = paste(Ouput_path, "/GLMM_Summary_Master.csv", sep = ""), quote = F, row.names = F)


#### End ####