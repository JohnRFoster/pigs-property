library(dplyr)
library(tidyr)
library(readr)
library(lubridate)
library(rsample)
library(recipes)
library(xgboost)

set.seed(123)

cutoff_date <- ymd("2023-12-31")

# config_name <- "hpc_dev"
config_name <- "default"
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

filename <- file.path(data_repo, 'ecoregions', 'Cnty.lower48.EcoRegions.Level2.shp')
ecoregions <- terra::vect(filename) |>
  as_tibble() |>
  rename(st_name = STATE_NAME,
         cnty_name = NAME,
         county_code = FIPS,
         ecoregion = NA_L2NAME) |>
  select(st_name, county_code, ecoregion) |>
  mutate(st_name = toupper(st_name))

first_flag <- NA

yearly_summaries <- data_mis |>
  mutate(total_take = if_else(is.na(total_take), 0, total_take),
         n_events = if_else(is.na(n_events), 0, n_events)) |>
  left_join(ecoregions) |>
  group_by(propertyID, year, st_name, county_code, property_area_km2, ecoregion, farm_bill) |>
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
         cum_take = cumsum(sum_take),
         cum_events = cumsum(sum_events),
         sum_take_density = sum_take / property_area_km2,
         avg_take_density = avg_take / property_area_km2,
         sum_events_density = sum_events / property_area_km2,
         avg_events_density = avg_events / property_area_km2,
         cum_take_density = cum_take / property_area_km2,
         cum_events_density = cum_events / property_area_km2) |>
  ungroup()

method_per_year <- data |>
  mutate(year = year(end_dates)) |>
  group_by(propertyID, year, method) |>
  summarise(events = n(),
            units = sum(trap_count)) |>
  ungroup() |>
  mutate(units_per_event = units / events)

events_by_method_per_year <- method_per_year |>
  select(propertyID, year, events, method) |>
  mutate(n_method_events_year = paste0("n_", method, "_events")) |>
  select(-method) |>
  pivot_wider(names_from = n_method_events_year,
              values_from = events,
              values_fill = 0)

units_by_method_per_year <- method_per_year |>
  select(propertyID, year, units, method) |>
  mutate(n_method_units_year = paste0("n_", method, "_units")) |>
  select(-method) |>
  pivot_wider(names_from = n_method_units_year,
              values_from = units,
              values_fill = 0)

units_per_event_year_by_method_per_year <- method_per_year |>
  select(propertyID, year, units_per_event, method) |>
  mutate(n_method_units_per_event_year = paste0("n_", method, "_units_per_event")) |>
  select(-method) |>
  pivot_wider(names_from = n_method_units_per_event_year,
              values_from = units_per_event,
              values_fill = 0)

data_ml <- yearly_summaries |>
  left_join(data_obs) |>
  left_join(events_by_method_per_year) |>
  left_join(units_by_method_per_year) |>
  left_join(units_per_event_year_by_method_per_year) |>
  mutate(
    y = med_density,
    year = ymd(paste0(year, "-01-01")),
    farm_bill = as.factor(farm_bill),
    propertyID = factor(propertyID),
    st_name = factor(st_name),
    county_code = factor(county_code),
    state_year = factor(paste(st_name, year)),
    county_year = factor(paste(county_code, year)),
    property_year = factor(paste(propertyID, year)),
    ecoregion = factor(ecoregion),
    eco_year = factor(paste(ecoregion, year))) |>
  select(-med_density)

properties <- unique(data_ml$propertyID)
n_props <- length(properties)
n_test_props <- 100
test_draws <- sample.int(n_props, n_test_props)

train_properties <- properties[-test_draws]
test_properties <- properties[test_draws]

df_train <- data_ml |>
  filter(propertyID %in% train_properties) |>
  group_by(propertyID) |>
  filter(year < max(year)) |>
  ungroup()

test1 <- data_ml |>
  filter(propertyID %in% train_properties) |>
  group_by(propertyID) |>
  filter(year == max(year)) |>
  ungroup()

df_test <- data_ml |>
  filter(propertyID %in% test_properties) |>
  bind_rows(test1)

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
    step_novel(eco_year, state_year, county_year, st_name,
               county_code, ecoregion, propertyID, property_year) |>
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
    step_interact(terms = ~cum_take:cum_events) |>
  # Normalization steps (center, scale, range, etc; does not apply)
  # Multivariate transformation (e.g. PCA, spatial sign, etc; does not apply)

  # filter
    step_nzv(all_predictors())

  return(blueprint)

}

blueprint <- my_recipe(df_train)
prepare <- prep(blueprint, training = df_train)
baked_train <- bake(prepare, new_data = df_train)
baked_test <- bake(prepare, new_data = df_test)

X <- baked_train |>
  select(-y) |>
  as.matrix()
Y <- baked_train |> pull(y)
hist(Y)

# hyperparameter grid
hyper_grid <- expand_grid(
  eta = c(0.3, 0.1, 0.05, 0.01, 0.005),
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

# hyper_grid <- hyper_grid[1:10, ]
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

df_pred <- make_prediction(fit, baked_test)
df_pred$y <- baked_test$y

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
  tidy_x = tidy(prepare, number = 4)

)

dest <- config$out_delta
filename <- file.path(dest, "ml_YeoJohnsonALL.rds")
write_rds(out_list, filename)


