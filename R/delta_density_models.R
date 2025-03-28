library(dplyr)
library(tidyr)
library(readr)
library(lubridate)
library(rsample)
library(recipes)
library(xgboost)

set.seed(123)

cutoff_date <- ymd("2023-12-31")

config_name <- "hpc_dev"
# config_name <- "default"
config <- config::get(config = config_name)

source("R/functions_data.R")

# ===================================================
# Data ingest ----
# ===================================================

data_repo <- config$data_repo

## observation covariates ----
file <- file.path(data_repo, config$file_land)
data_obs <- get_obs_covars(file) |>
  select(rural.road.density, prop.pub.land, mean.ruggedness, mean.canopy.density, county_code)


## join with farm bill properties ----
data_farm_bill <- read_csv(file.path(data_repo, "All_FB_Agreements_long_2024-05-30.csv"))
farm_bill_properties <- data_farm_bill |>
  rename(alws_agrprop_id = propertyID) |>
  select(-agreement_name, -property_name) |>
  mutate(farm_bill = 1)

if_dir <- "11_posterior"
top_dir <- config$out_dir

posterior_path <- file.path(top_dir, if_dir, "modelData.rds")
data <- read_rds(posterior_path) |>
  mutate(method = if_else(method == "FIREARMS", "Firearms", method),
         method = if_else(method == "FIXED WING", "Fixed wing", method),
         method = if_else(method == "HELICOPTER", "Helicopter", method),
         method = if_else(method == "SNARE", "Snare", method),
         method = if_else(method == "TRAPS", "Trap", method))

posterior_path <- file.path(top_dir, if_dir, "densitySummaries.rds")
density <- read_rds(posterior_path)

n_return <- data |>
  select(propertyID, primary_period, method) |>
  distinct() |>
  pivot_wider(names_from = method,
              values_from = method) |>
  unite(method,
        -c(propertyID, primary_period),
        sep = ", ",
        na.rm = TRUE) |>
  group_by(propertyID, method) |>
  mutate(return_interval = c(0, diff(primary_period))) |>
  ungroup() |>
  rename(methods_used = method)

data_grouped <- data |>
  group_by(propertyID, agrp_prp_id, start_dates, end_dates, st_name, cnty_name, farm_bill,
           alws_agrprop_id, property, primary_period, property_area_km2, county_code) |>
  summarise(total_take = sum(take),
            total_effort_per = sum(effort_per),
            n_events = n()) |>
  ungroup() |>
  mutate(take_density = total_take / property_area_km2,
         farm_bill = if_else(is.na(farm_bill), 0, farm_bill)) |>
  left_join(density) |>
  left_join(n_return) |>
  mutate(year = year(end_dates))

data_mis <- data_grouped |>
  filter(end_dates <= cutoff_date) |>
  mutate(abundance_estimate = round(`0.5` * property_area_km2)) |>
  rename(density_estimate = `0.5`) |>
  select(propertyID, end_dates, year, st_name, cnty_name, county_code, farm_bill,
         property_area_km2, total_take, take_density, density_estimate, abundance_estimate,
         n_events, methods_used)

# filename <- file.path(data_repo, 'ecoregions', 'Cnty.lower48.EcoRegions.Level2.shp')
# ecoregions <- terra::vect(filename) |>
#   as_tibble() |>
#   rename(st_name = STATE_NAME,
#          cnty_name = NAME,
#          county_code = FIPS,
#          ecoregion = NA_L2NAME) |>
#   select(st_name, county_code, ecoregion) |>
#   mutate(st_name = toupper(st_name))

first_flag <- NA

yearly_summaries <- data_mis |>
  mutate(total_take = if_else(is.na(total_take), 0, total_take),
         n_events = if_else(is.na(n_events), 0, n_events)) |>
  # left_join(ecoregions) |>
  group_by(propertyID, year, st_name, county_code, property_area_km2, farm_bill) |>
  summarise(med_density = median(density_estimate),
            sum_take = sum(total_take),
            avg_take = mean(total_take),
            sum_events = sum(n_events),
            avg_events = mean(n_events)
            ) |>
  arrange(propertyID, year) |>
  group_by(propertyID) |>
  mutate(#delta_density = c(first_flag, diff(med_density)),
         delta_take = c(first_flag, diff(sum_take)),
         delta_events = c(first_flag, diff(sum_events)),
         sum_take_density = sum_take / property_area_km2,
         avg_take_density = avg_take / property_area_km2,
         sum_events_density = sum_events / property_area_km2,
         avg_events_density = avg_events / property_area_km2) |>
  ungroup()

method_per_year <- data |>
  mutate(year = year(end_dates)) |>
  group_by(propertyID, year, method) |>
  summarise(events = n(),
            units = sum(trap_count),
            take = sum(take)) |>
  ungroup() |>
  mutate(units_per_event = units / events,
         take_per_unit = take / units,
         take_per_event = take / events,
         method = if_else(method == "Fixed wing", "fixedWing", method))

make_wide_feature <- function(df, col){
  df |>
    select(propertyID, year, method, all_of(col)) |>
    mutate(newname = paste0(method, "_", col)) |>
    select(-method) |>
    pivot_wider(names_from = newname,
                values_from = all_of(col),
                values_fill = 0)
}

events_by_method <- make_wide_feature(method_per_year, "events")
units_by_method <- make_wide_feature(method_per_year, "units")
take_by_method <- make_wide_feature(method_per_year, "take")
units_per_event_by_method <- make_wide_feature(method_per_year, "units_per_event")
take_per_unit_by_method <- make_wide_feature(method_per_year, "take_per_unit")
take_per_event_by_method <- make_wide_feature(method_per_year, "take_per_event")



data_ml <- yearly_summaries |>
  left_join(data_obs) |>
  left_join(events_by_method) |>
  left_join(units_by_method) |>
  left_join(take_by_method) |>
  left_join(units_per_event_by_method) |>
  left_join(take_per_unit_by_method) |>
  left_join(take_per_event_by_method) |>
  mutate(
    y = med_density,
    year = ymd(paste0(year, "-01-01")),
    farm_bill = as.factor(farm_bill),
    propertyID = factor(propertyID),
    st_name = factor(st_name),
    county_code = factor(county_code),
    state_year = factor(paste(st_name, year)),
    county_year = factor(paste(county_code, year)),
    property_year = factor(paste(propertyID, year))) |>
    # ecoregion = factor(ecoregion),
    # eco_year = factor(paste(ecoregion, year))) |>
  select(-med_density)

all_properties <- unique(data_ml$propertyID)

# these properties will become testing properties
single_year_properties <- data_ml |>
  group_by(propertyID) |>
  count() |>
  ungroup() |>
  filter(n == 1) |>
  pull(propertyID) |>
  unique()

# set of properties to draw train/test
multi_year_properties <- data_ml |>
  filter(!propertyID %in% single_year_properties) |>
  pull(propertyID) |>
  unique()

n_single_year <- length(single_year_properties)
n_multi_year <- length(multi_year_properties)
n_all <- length(all_properties)

assertthat::assert_that(
  n_single_year + n_multi_year == n_all
  )

n_test_single <- round(0.1 * n_single_year)
test_draws_single <- sample.int(n_single_year, n_test_single)

n_test_multi <- round(0.1 * n_multi_year)
test_draws_multi <- sample.int(n_multi_year, n_test_multi)

train_properties <- c(
  single_year_properties[-test_draws_single],
  multi_year_properties[-test_draws_multi]
)

test_properties <- c(
  single_year_properties[test_draws_single],
  multi_year_properties[test_draws_multi]
)

# training properties with one year
train1 <- data_ml |>
  filter(propertyID %in% single_year_properties[-test_draws_single])

# training properties with multiple years without their last year
train2 <- data_ml |>
  filter(propertyID %in% multi_year_properties[-test_draws_multi]) |>
  group_by(propertyID) |>
  filter(year < max(year)) |>
  ungroup()

df_train <- rbind(train1, train2) |>
  mutate(partition = "train")

# testing properties with one year
test1 <- data_ml |>
  filter(propertyID %in% single_year_properties[test_draws_single]) |>
  mutate(partition = "test_single_year")

# test properties (last year from training set) with multiple years
test2 <- data_ml |>
  filter(propertyID %in% multi_year_properties[-test_draws_multi]) |>
  group_by(propertyID) |>
  filter(year == max(year)) |>
  ungroup() |>
  mutate(partition = "test_last_year")

# held out properties because they only have one year
test3 <- data_ml |>
  filter(propertyID %in%  multi_year_properties[test_draws_multi]) |>
  mutate(partition = "test_holdout")

df_test <- bind_rows(test1, test2, test3)

p_train <- round(nrow(df_train) / nrow(data_ml), 2) * 100
p_test <- round(nrow(df_test) / nrow(data_ml), 2) * 100

message("Train test split across data: ", p_train, ":", p_test)

# move to functions script
my_recipe <- function(data){
  require(rsample)
  require(recipes)

  blueprint <- recipe(y ~ ., data = data) |>

  # recommended order from recipes package

  # Impute (does not apply)
  # Handle factor levels
    step_novel(all_nominal_predictors()) |>
    step_date(year, features = "year", keep_original_cols = FALSE) |>

  # Individual transformations for skewness and other issues
    step_YeoJohnson(all_outcomes()) |>
    step_YeoJohnson(all_numeric_predictors()) |>

  # Discretize (if needed and if you have no other choice; does not apply)
  # Create dummy variables
    step_dummy(all_nominal_predictors()) |>

  # Create interactions
    step_interact(terms = ~sum_take:sum_events) |>
    step_interact(terms = ~avg_take:avg_events) |>
    step_interact(terms = ~sum_take_density:sum_events_density) |>
    step_interact(terms = ~avg_take_density:avg_events_density) |>

    step_interact(terms = ~sum_take:rural.road.density) |>
    step_interact(terms = ~sum_take:prop.pub.land) |>
    step_interact(terms = ~sum_take:mean.ruggedness) |>
    step_interact(terms = ~sum_take:mean.canopy.density) |>

    step_interact(terms = ~sum_take_density:rural.road.density) |>
    step_interact(terms = ~sum_take_density:prop.pub.land) |>
    step_interact(terms = ~sum_take_density:mean.ruggedness) |>
    step_interact(terms = ~sum_take_density:mean.canopy.density) |>

    step_interact(terms = ~avg_take:rural.road.density) |>
    step_interact(terms = ~avg_take:prop.pub.land) |>
    step_interact(terms = ~avg_take:mean.ruggedness) |>
    step_interact(terms = ~avg_take:mean.canopy.density) |>

    step_interact(terms = ~avg_take_density:rural.road.density) |>
    step_interact(terms = ~avg_take_density:prop.pub.land) |>
    step_interact(terms = ~avg_take_density:mean.ruggedness) |>
    step_interact(terms = ~avg_take_density:mean.canopy.density) |>

    step_interact(terms = ~rural.road.density:starts_with("Helicopter")) |>
    step_interact(terms = ~prop.pub.land:starts_with("Helicopter")) |>
    step_interact(terms = ~mean.ruggedness:starts_with("Helicopter")) |>
    step_interact(terms = ~mean.canopy.density:starts_with("Helicopter")) |>

    step_interact(terms = ~rural.road.density:starts_with("fixedWing")) |>
    step_interact(terms = ~prop.pub.land:starts_with("fixedWing")) |>
    step_interact(terms = ~mean.ruggedness:starts_with("fixedWing")) |>
    step_interact(terms = ~mean.canopy.density:starts_with("fixedWing")) |>

    step_interact(terms = ~rural.road.density:starts_with("Trap")) |>
    step_interact(terms = ~prop.pub.land:starts_with("Trap")) |>
    step_interact(terms = ~mean.ruggedness:starts_with("Trap")) |>
    step_interact(terms = ~mean.canopy.density:starts_with("Trap")) |>

    step_interact(terms = ~rural.road.density:starts_with("Snare")) |>
    step_interact(terms = ~prop.pub.land:starts_with("Snare")) |>
    step_interact(terms = ~mean.ruggedness:starts_with("Snare")) |>
    step_interact(terms = ~mean.canopy.density:starts_with("Snare")) |>

    step_interact(terms = ~rural.road.density:starts_with("Firearms")) |>
    step_interact(terms = ~prop.pub.land:starts_with("Firearms")) |>
    step_interact(terms = ~mean.ruggedness:starts_with("Firearms")) |>
    step_interact(terms = ~mean.canopy.density:starts_with("Firearms")) |>

  # Normalization steps (center, scale, range, etc; does not apply)
  # Multivariate transformation (e.g. PCA, spatial sign, etc; does not apply)

  # filter
    step_nzv(all_predictors())

  return(blueprint)

}

df_blueprint <- df_train |> select(-partition)
blueprint <- my_recipe(df_blueprint)
prepare <- prep(blueprint, training = df_blueprint)
baked_train <- bake(prepare, new_data = df_blueprint)
glimpse(baked_train)

X <- baked_train |>
  select(-y) |>
  as.matrix()
Y <- baked_train |> pull(y)
hist(Y)

# hyperparameter grid
hyper_grid <- expand_grid(
  eta = config$eta,
  # eta = c(0.3, 0.1, 0.05, 0.01, 0.005),
  max_depth = 3:8,
  min_child_weight = 0.5,
  subsample = 0.5,
  colsample_bytree = 0.5,
  gamma = c(0, 1, 10, 100, 1000),
  lambda = c(0, 1e-2, 0.1, 1, 100, 1000),
  alpha = c(0, 1e-2, 0.1, 1, 100, 1000),
  rmse = 0,
  trees = 0
)

if(config_name == "default"){
  hyper_grid <- hyper_grid[1:10, ]
}

objective <- "reg:squarederror"

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

test2 <- df_test |> select(-partition)
baked_test <- bake(prepare, new_data = test2)

df_pred <- make_prediction(fit, baked_test)
df_pred$y <- baked_test$y
df_pred$partition <- df_test$partition

rmse <- sqrt(mean((baked_test$y - df_pred$pred)^2))
cc <- cor(baked_test$y, df_pred$pred)
r2 <- cc^2
vi <- xgb.importance(model = fit)
vi$gainRelative <- vi$Gain / max(vi$Gain)

message("===============================")
message("RMSE: ", round(rmse, 2))
message("COR: ", round(cc, 2))
message("R2: ", round(r2, 2))
message("===============================")

plot(df_pred$y, df_pred$pred)
abline(0, 1)

out_list <- list(
  baked_train = baked_train,
  baked_test = df_pred,
  test_test_rmse = rmse,
  test_test_r2 = r2,
  vi = vi,
  hyper_grid = hyper_grid,
  tidy_y = tidy(prepare, number = 3),
  tidy_x = tidy(prepare, number = 4),
  raw_train = df_train,
  data = bind_rows(df_train, df_test)

)

dest <- config$out_delta
filename <- file.path(dest, paste0(config$eta, "ml_YeoJohnsonALL.rds"))
write_rds(out_list, filename)


