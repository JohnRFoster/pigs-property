library(dplyr)
library(tidyr)
library(readr)
library(rsample)
library(recipes)
library(xgboost)

source("R/functions_data.R")
source("R/functions_misc.R")

config_name <- "hpc_dev"
# config_name <- "default"
config <- config::get(config = config_name)

# args <- commandArgs(trailingOnly = TRUE)
task_id <- 1

## choose response variable
responses <- c("0.5", "variance")
response <- responses[task_id]

objective <- "reg:squarederror"

iterative_dir <- config$out_dir
num_dirs <- list.files(iterative_dir)
post_dirs <- grep("posterior", num_dirs, value = TRUE)
iters <- sort(as.numeric(stringr::str_extract(post_dirs, "^\\d*(?=\\_)")), decreasing = TRUE)
last_dir <- max(iters)
read_dir <- paste0(last_dir, "_posterior")

# check if posterior exists
post <- file.path(iterative_dir, read_dir, "densitySummaries.rds")

does_not_exist <- !file.exists(post)
i <- 2
while(does_not_exist){
  read_dir <- paste0(iters[i], "_posterior")
  post <- file.path(iterative_dir, read_dir, "densitySummaries.rds")
  does_not_exist <- !file.exists(post)
  i <- i + 1
}

path <- file.path(iterative_dir, read_dir, "modelData.rds")
data <- read_rds(path) |>
  mutate(method = if_else(method == "FIREARMS", "Firearms", method),
         method = if_else(method == "FIXED WING", "Fixed wing", method),
         method = if_else(method == "HELICOPTER", "Helicopter", method),
         method = if_else(method == "SNARE", "Snare", method),
         method = if_else(method == "TRAPS", "Trap", method))

path <- file.path(iterative_dir, read_dir, "densitySummaries.rds")
density <- read_rds(path)

# want landcover for ML model
land_cover <- data |>
  select(propertyID, c_road_den, c_rugged, c_canopy) |>
  distinct()

# want ecoregion data for ML
data_repo <- config$data_repo

filename <- file.path(data_repo, 'ecoregions', 'Cnty.lower48.EcoRegions.Level2.shp')
ecoregions <- terra::vect(filename) |>
  as_tibble() |>
  rename(st_name = STATE_NAME,
         cnty_name = NAME,
         ecoregion = NA_L2NAME) |>
  select(st_name, cnty_name, ecoregion) |>
  mutate(st_name = toupper(st_name),
         cnty_name = toupper(cnty_name))

data <- data |> left_join(land_cover)
model_data <- group_join_for_ml(data, ecoregions)
model_data <- model_data |> left_join(density)


features <- c("cnty_name", "st_name", "farm_bill", "property_area_km2",
              "total_take", "return_interval", "take_density", "methods_used", "n_units",
              "c_road_den", "c_rugged", "c_canopy", "ecoregion", "effort", "effort_per")

n_strata <- 3

low_outlier <- quantile(model_data$`0.5`, c(0.025))
high_outlier <- quantile(model_data$`0.5`, c(0.975))
min_not_0 <- model_data |> filter(mean > 0) |> pull(mean) |> min()

make_c_strata <- function(df){
  df |>
    mutate(road_dens_strata = if_else(c_road_den < quantile(c_road_den, 0.25), "low", "avg"),
           road_dens_strata = if_else(c_road_den > quantile(c_road_den, 0.75), "high", road_dens_strata),
           canopy_strata = if_else(c_canopy < quantile(c_canopy, 0.25), "low", "avg"),
           canopy_strata = if_else(c_canopy > quantile(c_canopy, 0.75), "high", canopy_strata),
           rugged_strata = if_else(c_rugged < quantile(c_rugged, 0.25), "low", "avg"),
           rugged_strata = if_else(c_rugged > quantile(c_rugged, 0.75), "high", rugged_strata))
}

ml_data <- model_data |>
  select(all_of(c(response, features))) |>
  rename(y = all_of(response)) |>
  mutate(y = y ^ (1/3)) #|>
  # make_c_strata()


split <- initial_split(ml_data, prop = 0.8)
df_train <- training(split)
df_test <- testing(split)

blueprint <- my_recipe(df_train)
prepare <- prep(blueprint, training = df_train)
baked_train <- bake(prepare, new_data = df_train)
baked_test <- bake(prepare, new_data = df_test)

X <- baked_train |>
  select(-y) |>
  as.matrix()
Y <- baked_train |> pull(y)
hist(Y)
set.seed(123)

hyp_vec <- c(0, 0.01, 0.1, 1, 10, 100)

# hyperparameter grid
hyper_grid <- expand_grid(
  eta = 0.1,
  max_depth = 3,
  min_child_weight = 0.5,
  subsample = 0.5,
  colsample_bytree = 0.5,
  gamma = hyp_vec,
  lambda = hyp_vec,
  alpha = hyp_vec,
  rmse = 0,
  trees = 0
)

# hyper_grid <- hyper_grid[1:10, ]

pb <- txtProgressBar(min = 1, max = nrow(hyper_grid), style = 1)
for(i in 1:nrow(hyper_grid)){

  m <- xgb.cv(
    data = X,
    label = Y,
    nrounds = 5000,
    objective = objective,
    metrics = "rmse",
    early_stopping_rounds = 50,
    nfold = 10,
    verbose = 0,
    params = list(
      eta = hyper_grid$eta[i],
      max_depth = hyper_grid$max_depth[i],
      min_child_weight = hyper_grid$min_child_weight[i],
      subsample = hyper_grid$subsample[i],
      colsample_bytree = hyper_grid$colsample_bytree[i],
      gamma = hyper_grid$gamma[i],
      lambda = hyper_grid$lambda[i],
      alpha = hyper_grid$alpha[i]
    )
  )

  hyper_grid$rmse[i] <- min(m$evaluation_log$test_rmse_mean)
  hyper_grid$trees[i] <- m$best_iteration

  setTxtProgressBar(pb, i)

}
close(pb)

tune_grid <- hyper_grid |>
  filter(rmse == min(rmse))

params <- list(
  eta = pull(tune_grid, eta),
  max_depth = pull(tune_grid, max_depth),
  min_child_weight = pull(tune_grid, min_child_weight),
  subsample = pull(tune_grid, subsample),
  colsample_bytree = pull(tune_grid, colsample_bytree),
  gamma = pull(tune_grid, gamma),
  lambda = pull(tune_grid, lambda),
  alpha = pull(tune_grid, alpha)
)

# train final model
fit <- xgboost(
  params = params,
  data = X,
  label = Y,
  nrounds = pull(tune_grid, trees),
  objective = objective,
  verbose = 0
)

make_prediction <- function(model, new_data){

  if("y" %in% colnames(new_data)){
    new_data <- new_data |> select(-y)
  }

  pred <- predict(fit, as.matrix(new_data), interval="confidence")

  new_data |> mutate(pred = pred)
}

df_pred <- make_prediction(fit, baked_test)
df_pred$y <- baked_test$y

rmse <- sqrt(mean((baked_test$y - df_pred$pred)^2))
r2 <- cor(baked_test$y, df_pred$pred)^2
vi <- xgb.importance(model = fit)
vi$gainRelative <- vi$Gain / max(vi$Gain)


message("===============================")
message("RMSE: ", round(rmse, 2))
message("R2: ", round(r2, 2))
message("===============================")

plot(df_pred$y, df_pred$pred)
abline(0, 1)

plot(df_pred$y, df_pred$pred)
abline(0, 1)

out_list <- list(
  baked_train = baked_train,
  baked_test = df_pred,
  train_test_rmse = rmse,
  train_test_r2 = r2,
  vi = vi,
  hyper_grid = hyper_grid
)


## out of sample predictions -------------------


## MIS data ----
file <- file.path(data_repo, config$file_mis)
interval <- config$interval
dev <- config$dev

all_take <- read_csv(file, show_col_types = FALSE) |>
  filter(start.date >= lubridate::ymd("2014-01-01")) |>
  mutate(cnty_name = if_else(grepl("ST ", cnty_name), gsub("ST ", "ST. ", cnty_name), cnty_name),
         cnty_name = if_else(grepl("KERN", cnty_name), "KERN", cnty_name)) |>
  mutate(property_area_km2 = round(property.size * 0.00404686, 2)) |>
  filter(property_area_km2 >= 1.8,
         st_name != "HAWAII") |>
  mutate(effort = if_else(cmp_name %in% c("TRAPS, CAGE", "SNARE"), cmp.days, cmp.hours),
         effort_per = effort / cmp.qty,
         cmp_name = if_else(cmp_name == "TRAPS, CAGE", "TRAPS", cmp_name)) |>
  rename(method = cmp_name,
         trap_count = cmp.qty) |>
  select(-wt_work_date, -hours, -cmp.hours, -cmp.days) |>
  distinct() |>
  mutate(propertyID = paste0(agrp_prp_id, "-", alws_agrprop_id)) |>
  arrange(propertyID, start.date, end.date)

data_mis <- create_primary_periods(all_take, interval, data_repo) |>
  resolve_duplicate() |>         # resolve duplicate property areas
  order_interval() |>            # determine midpoints from start/end dates
  order_stochastic() |>          # randomly order events
  rename(statefp = st_gsa_state_cd,
         countyfp = cnty_gsa_cnty_cd) |>
  mutate(countyfp = sprintf("%03d", countyfp),
         countyfp = ifelse(cnty_name == "HUMBOLDT (E)", "013", countyfp),
         county_code = as.numeric(paste0(statefp, countyfp)),
         county_code = sprintf("%05d", county_code)) |>
  rename(primary_period = timestep)

## observation covariates ----
file <- file.path(data_repo, config$file_land)
data_obs <- get_obs_covars(file)

## join MIS with observation covariates ----
data_join <- left_join(data_mis, data_obs,
                       by = join_by(county_code))

## filter missing states ----
data_join2 <- data_join #|>
  # filter(!st_name %in% c("CALIFORNIA", "ALABAMA", "ARIZONA", "ARKANSAS"))

# targets::tar_assert_true(!any(is.na(data_join2$c_road_den)))
# targets::tar_assert_true(!any(is.na(data_join2$c_rugged)))
# targets::tar_assert_true(!any(is.na(data_join2$c_canopy)))

## join with farm bill properties ----
data_farm_bill <- read_csv(file.path(data_repo, "All_FB_Agreements_long_2024-05-30.csv"))
farm_bill_properties <- data_farm_bill |>
  rename(alws_agrprop_id = propertyID) |>
  select(-agreement_name, -property_name) |>
  mutate(farm_bill = 1)

data_join3 <- left_join(data_join2, farm_bill_properties) |>
  mutate(property = as.numeric(as.factor(propertyID)),
         county = as.numeric(as.factor(county_code)))

bayes_fit_properties <- unique(data$propertyID)

data_ml_filter <- data_join3 |>
  filter(!propertyID %in% bayes_fit_properties)

oos_data <- group_join_for_ml(data_ml_filter, ecoregions)
df_oos <- oos_data |>
  select(all_of(features)) #|>
  # make_c_strata()

baked_oos <- bake(prepare, new_data = df_oos)

oos_pred <- make_prediction(fit, baked_oos)

out_list$oos_pred <- oos_pred

ml_dir <- config$out_ml
if(!dir.exists(ml_dir)) dir.create(ml_dir, recursive = TRUE, showWarnings = FALSE)
write_rds(out_list, file.path(ml_dir, paste0(response, "_mlFits.rds")))


