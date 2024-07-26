
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

group_join_for_ml <- function(df, df_landcover, df_ecoregions){

  # aggregate by primary period
  df_n_return <- return_intervals_methods_used(df)

  df |>
    group_by(propertyID, agrp_prp_id, start_dates, end_dates, st_name, cnty_name, farm_bill,
             alws_agrprop_id, property, primary_period, property_area_km2) |>
    summarise(total_take = sum(take)) |>
    ungroup() |>
    mutate(take_density = total_take / property_area_km2,
           farm_bill = as.character(if_else(is.na(farm_bill), 0, farm_bill)),
           cnty_name = if_else(cnty_name == "VIRGINIA BEACH CITY", "VIRGINIA BEACH", cnty_name),
           cnty_name = if_else(cnty_name == "DE WITT", "DEWITT", cnty_name)) |>
    left_join(df_n_return) |>
    left_join(df_landcover) |>
    left_join(df_ecoregions)
}
