# ------------------------------------------------------------------------------------------
# Kajanka J Mathiaparanam
# 30th October 2020
# BRT analysis to test the cumulative impacts of noise and light pollution on Coot abundance
# ------------------------------------------------------------------------------------------

curdir <- getwd()

# Load libraries
library(dplyr)
library(tidyverse)
library(gplots)
library(vegan)
library(Hmisc)
library(NADA)
library(survival)
library(gtools)
library(lattice)
library(gbm)
library(ggplot2)
library(dismo)
library(car)
library(psych)


#### 1) Functions ####

# Mode function
Mode <- function(x, na.rm = TRUE) {
  if(na.rm){
    x = x[!is.na(x)]
  }
  ux <- unique(x)
  return(ux[which.max(tabulate(match(x, ux)))])
}


#### 2) Assign variables and model properties ####

# Read in the data
BLA_data <- read.csv(paste(curdir, "/BLA_Prepped_Data_BRT.csv", sep = ""))

# Bird species
bird_species <- "Eurasian_Coot"

# Assign response and predictor variables
Resp_var <- "Mean_Count"
Pred_var <- c("Poi_Noise_RR", "Poi_Lmean13_18", "Area_m2", "Poi_Temp", "Poi_Rain", "WB_GrVeg_Cat")
Pred_var_rename <- c("Noise", "Light", "Area_m2", "Temp", "Rain", "Veg")

# Family
Mod_family <- "poisson"

# Bag fraction = 0.5, 0.75
BRT_bf <- 0.5

# Tree complexity = 1, 2, 3, 4, 5
BRT_tc <- 2

# Learning rate = 0.01, 0.005, 0.001
BRT_lr <- 0.01

# Waterbody Area filter
Area_filt <- 10^6

# Noise Filter
Noi_thresh <- 30


#### 3) Preparation for model ####

# Create folder for model summary csv and plots
Plot_folder_name <- str_replace_all(paste("BRT_", Sys.time(), sep = ""), ":", "-")
Output_path <- paste(curdir, "/Model_Outputs/", Plot_folder_name, sep = "")
dir.create(Output_path)

# Labels of the BRT summary table
BLA_Model_df_cols <- c("Learning_rate", "Tree_complexity", "Bag_fraction", "Trees_n",
                       "Lmean_influence", "Noise_influence", "WB_Area_influence", "Temperature_influence",
                       "Rainfall_influence", "Light_Noise_Interaction")

# Set seed if necessary
set.seed(47)


#### 4) Subsetting, cleaning and transformations ####

# Filter response and predictor columns
BLA_data_temp <- BLA_data %>% dplyr::select(all_of(c(Resp_var, Pred_var)))

# Rename columns to apply filters + plot etc. easily
names(BLA_data_temp)[2:(1 + length(Pred_var))] <- Pred_var_rename

# Apply Noise filter/transformation
BLA_data_temp <- BLA_data_temp %>% mutate(Noise = if_else(Noise < Noi_thresh, NA, Noise))

# Apply Area filter
BLA_data_temp <- BLA_data_temp %>% filter(Area_m2 <= Area_filt)

# Convert Categorical variables to factors
Cat_var <- c("Veg")
for (cat_var_i in 1:length(Cat_var)) {
  BLA_data_temp[[Cat_var[cat_var_i]]] <- as.factor(BLA_data_temp[[Cat_var[cat_var_i]]])
}


#### 5) Run Model and optimise ####

# Initiate Model
BLA_mod <- gbm.step(data = BLA_data_temp,
                    gbm.x = Pred_var_rename,
                    gbm.y = 1,
                    family = Mod_family,
                    tree.complexity = BRT_tc, learning.rate = BRT_lr, bag.fraction = BRT_bf)

# Rerun model with slower lr if Trees_n = NA or < 1000, or if BLA_mod - null, until lr = 10^-7
BRT_lr2 <- BRT_lr

if (!is.null(BLA_mod)) {
  find.int <- gbm.interactions(BLA_mod)
  Trees_n <- find.int$gbm.call$best.trees
} else {
  Trees_n <- NA # Dummy
}

while ((is.null(BLA_mod) | (!is.null(BLA_mod) & Trees_n < 1000)) & BRT_lr2 > 10^-7) {
  BRT_lr2 <- BRT_lr2/10

  BLA_mod <- gbm.step(data = BLA_data_temp,
                      gbm.x = Pred_var_rename,
                      gbm.y = 1,
                      family = Mod_family,
                      tree.complexity = BRT_tc, learning.rate = BRT_lr2, bag.fraction = BRT_bf)

  if (!is.null(BLA_mod)) {
    find.int <- gbm.interactions(BLA_mod)
    Trees_n <- find.int$gbm.call$best.trees
  }
}


#### 6) Model Summary ####

# Model summary
BLA_mod_summary <- summary(BLA_mod)

# Number of trees
find.int <- gbm.interactions(BLA_mod)
Trees_n <- find.int$gbm.call$best.trees

# Influences
Lmean_inf <- round(BLA_mod_summary$rel.inf[BLA_mod_summary$var == "Light"], digits = 2)
Nmax_inf <- round(BLA_mod_summary$rel.inf[BLA_mod_summary$var == "Noise"], digits = 2)
WB_Area_inf <- round(BLA_mod_summary$rel.inf[BLA_mod_summary$var == "Area_m2"], digits = 2)
Temp_inf <- round(BLA_mod_summary$rel.inf[BLA_mod_summary$var == "Temp"], digits = 2)
Rain_inf <- round(BLA_mod_summary$rel.inf[BLA_mod_summary$var == "Rain"], digits = 2)

# Interactions
Noi_Lig_int <- max(find.int$rank.list[find.int$rank.list["var1.names"] == "Light" & find.int$rank.list["var2.names"] == "Noise" |
                                        find.int$rank.list["var1.names"] == "Noise" & find.int$rank.list["var2.names"] == "Light" , "int.size"])
if (Noi_Lig_int == "-Inf" | is.na(Noi_Lig_int)) {
  Noi_Lig_int <- as.numeric("NA")
}

# Save model summary
BLA_Model_df <- as.data.frame(cbind(BRT_lr2, BRT_tc, BRT_bf, Trees_n,
                        Lmean_inf, Nmax_inf, WB_Area_inf, Temp_inf, Rain_inf, Noi_Lig_int))
names(BLA_Model_df) <- c(BLA_Model_df_cols)
BLA_Model_df <- as.data.frame(t(BLA_Model_df))
names(BLA_Model_df) <- c("Value")
write.csv(BLA_Model_df, file = paste(Output_path, "/BRT_Model_Summary.csv", sep = ""), quote = F, row.names = T)


#### 7) Model predictions ####

for (i in 1:length(Pred_var_rename)) {
  
  Resp_Pred_matrix <- plot.gbm(BLA_mod, i.var = Pred_var_rename[i], n.trees = Trees_n, return.grid = TRUE)
  
  # Create predictor matrix using the predictor response matrix above
  Pred_mean <- Pred_var_rename[!(Pred_var_rename == "Veg")]
  Pred_mean <- Pred_mean[!Pred_mean %in% Pred_var_rename[i]]
  Mean_matrix <- t(apply(BLA_data_temp[Pred_mean], 2, mean, na.rm = T))
  
  Pred_mode <- Pred_var_rename[Pred_var_rename %in% Cat_var]
  Pred_mode <- Pred_mode[!Pred_mode %in% Pred_var_rename[i]]
  Mode_matrix <- t(apply(BLA_data_temp[Pred_mode], 2, Mode, na.rm = T))
  
  Predictor_matrix <- cbind(Resp_Pred_matrix[Pred_var_rename[i]], Mean_matrix, Mode_matrix)
  
  # Get predicted responses
  Pred_response <- gbm::predict.gbm(BLA_mod, Predictor_matrix, n.trees = Trees_n, type = "response")
  
  # Combine the predictor matrix and predicted responses
  Resp_Pred_matrix <- cbind(Predictor_matrix, Pred_response)
  names(Resp_Pred_matrix)[length(Resp_Pred_matrix)] <- Resp_var
  
  # Save the predictor-response matrix as a csv and in the environment
  write.csv(Resp_Pred_matrix, file = paste(Output_path, "/BRT_Prediction_", Pred_var_rename[i], ".csv", sep = ""), quote = F, row.names = F)
  
  assign(paste0("Resp_Pred_matrix_", i), Resp_Pred_matrix)
  
  # Find the min and max of predicted responses for a common scale
  if (i == 1) {
    Resp_scale_min <- min(Pred_response, na.rm = T)
    Resp_scale_max <- max(Pred_response, na.rm = T)
  } else {
    if (Resp_scale_min > min(Pred_response, na.rm = T)) {
      Resp_scale_min <- min(Pred_response, na.rm = T)
    }
    if (Resp_scale_max < max(Pred_response, na.rm = T)) {
      Resp_scale_max <- max(Pred_response, na.rm = T)
    }
  }
}


#### 8) Save plots ####

# Save pdf
pdf(paste(Output_path, "/Plots.pdf", sep = ""))

  # 8A) Partial dependence plots
  gbm.plot(BLA_mod, plot.layout = c(2, 2), smooth = T)

  # 8B) Prediction plots
  par(mfrow = c(2, 2))
  
  for (i in 1:length(Pred_var_rename)) {
    
    # Get the corresponding predictor-response matrix
    Resp_Pred_matrix <- get(paste0("Resp_Pred_matrix_", i))

    # Plot - without points
    plot(Resp_Pred_matrix[[Pred_var_rename[i]]], Resp_Pred_matrix[[Resp_var]], cex = 0,
         ylim = c(Resp_scale_min, Resp_scale_max), xlab = Pred_var_rename[i], ylab = Resp_var)
    
    # Draw the trend for non categorical
    if (!is.element(Pred_var_rename[i], Cat_var)) {
      Resp_spline <- smooth.spline(Resp_Pred_matrix[[Pred_var_rename[i]]], Resp_Pred_matrix[[Resp_var]])
      lines(Resp_spline)
    }
  }

  # 8C) 3D plot with noise and light interaction
  par(mfrow = c(1, 1))
  
  fitted_min <- min(BLA_mod$fitted)
  fitted_max <- max(BLA_mod$fitted)
  gbm.perspec(BLA_mod, 1, 2, theta = 135, phi = 10, smooth = "average", cex.lab = 1.2)

dev.off()


#### End ####