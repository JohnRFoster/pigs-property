
return_intervals_methods_used <- function(df){
  df |>
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
}


my_recipe <- function(data){
  require(rsample)
  require(recipes)

  blueprint <- recipe(y ~ ., data = data) |>
    step_novel(methods_used, st_name, cnty_name, ecoregion) |>
    step_dummy(all_nominal_predictors()) |>
    step_nzv(all_predictors()) |>
    step_center(all_numeric_predictors()) |>
    step_scale(all_numeric_predictors())

    return(blueprint)

}

group_join_for_ml <- function(df, df_ecoregions){

  # aggregate by primary period
  df_n_return <- return_intervals_methods_used(df)

  df_landcover <- df |>
    select(propertyID, c_road_den, c_rugged, c_canopy) |>
    distinct()

  df |>
    group_by(propertyID, agrp_prp_id, start_dates, end_dates, st_name, cnty_name, farm_bill,
             alws_agrprop_id, property, primary_period, property_area_km2) |>
    summarise(total_take = sum(take),
              n_units = sum(trap_count),
              effort_per = sum(effort_per),
              effort = sum(effort)) |>
    ungroup() |>
    mutate(take_density = total_take / property_area_km2,
           farm_bill = as.character(if_else(is.na(farm_bill), 0, farm_bill)),
           cnty_name = if_else(cnty_name == "VIRGINIA BEACH CITY", "VIRGINIA BEACH", cnty_name),
           cnty_name = if_else(cnty_name == "DE WITT", "DEWITT", cnty_name)) |>
    left_join(df_n_return) |>
    left_join(df_landcover) |>
    left_join(df_ecoregions)
}

land_cover_class <- function(df){
  df |>
    mutate(canopy_class = if_else(c_canopy < -1, "low", "med"),
           canopy_class = if_else(c_canopy > 1, "high", canopy_class),
           rugged_class = if_else(c_rugged < -1, "low", "med"),
           rugged_class = if_else(c_rugged > 1, "high", rugged_class),
           road_den_class = if_else(c_road_den < -1, "low", "med"),
           road_den_class = if_else(c_road_den > 1, "high", road_den_class))
}

trap_snare_lpa <- function(log_rho, log_gamma, p_unique, effort_per, n_trap_m1){
  log(pi) +
    (2 * (log_rho + log(effort_per) -
            log(exp(log_gamma) + effort_per))) +
    log(1 + (p_unique * n_trap_m1))
}

shooting_lpa <- function(log_rho, effort_per){
  log_rho + log(effort_per)
}
