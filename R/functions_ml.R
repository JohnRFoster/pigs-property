create_events_per_year <- function(df){
  df |>
    group_by(propertyID, agrp_prp_id, property_area_km2, year) |>
    summarise(all_take = sum(take),
              all_events = n(),
              all_units = sum(trap_count)) |>
    mutate(all_take_density = all_take / property_area_km2,
           all_events_density = all_events / property_area_km2,
           all_units_density = all_units / property_area_km2,
           avg_take = all_take / all_events,
           avg_take_density = avg_take / property_area_km2)
}

create_methods_per_year <- function(df){
  df |>
    group_by(propertyID, year, method) |>
    summarise(events = n(),
              units = sum(trap_count),
              take = sum(take)) |>
    ungroup() |>
    mutate(units_per_event = units / events,
           take_per_unit = take / units,
           take_per_event = take / events)
}

create_method_feature <- function(df, col){

  tmp <- create_methods_per_year(df)

  tmp |>
    select(propertyID, year, method, all_of(col)) |>
    mutate(newname = paste0(method, "_", col)) |>
    select(-method) |>
    pivot_wider(names_from = newname,
                values_from = all_of(col),
                values_fill = 0)
}

create_all_events_per_year <- function(df){
  events_per_year <- create_events_per_year(df)
  events_by_method <- create_method_feature(df, "events")
  units_by_method <- create_method_feature(df, "units")
  take_by_method <- create_method_feature(df, "take")
  units_per_event_by_method <- create_method_feature(df, "units_per_event")
  take_per_unit_by_method <- create_method_feature(df, "take_per_unit")
  take_per_event_by_method <- create_method_feature(df, "take_per_event")

  left_join(
    events_per_year,
    events_by_method
  ) |>
    left_join(units_by_method) |>
    left_join(take_by_method) |>
    left_join(units_per_event_by_method) |>
    left_join(take_per_unit_by_method) |>
    left_join(take_per_event_by_method)

}

create_pp_data <- function(df, density_df, first_flag = NA){
  tmp <- df |>
    group_by(propertyID, property, agrp_prp_id, year,  start_dates, end_dates, st_name, cnty_name, farm_bill,
             alws_agrprop_id, primary_period, property_area_km2, county_code, trap_count) |>
    summarise(total_take_in_pp = sum(take),
              n_events_in_pp = n(),
              n_units_in_pp = sum(trap_count)) |>
    ungroup() |>
    mutate(take_density_in_pp = total_take_in_pp / property_area_km2,
           n_events_density_in_pp = n_events_in_pp / property_area_km2,
           n_units_density  = n_units_in_pp / property_area_km2,
           farm_bill = if_else(is.na(farm_bill), 0, farm_bill)) |>
    left_join(density_df) |>
    rename(density_estimate = `0.5`) |>
    select(propertyID, end_dates, year, st_name, cnty_name, county_code, farm_bill,
           property_area_km2, total_take_in_pp, take_density_in_pp, density_estimate,
           n_events_in_pp, n_events_density_in_pp)

  tmp |>
    group_by(propertyID, year, st_name, county_code, property_area_km2, farm_bill) |>
    summarise(med_density = median(density_estimate),
              all_take = sum(total_take_in_pp),
              avg_take_in_pp = mean(total_take_in_pp),
              all_events = sum(n_events_in_pp),
              avg_events_in_pp = mean(n_events_in_pp),
              n_sampled_pp = n()
    ) |>
    arrange(propertyID, year) |>
    group_by(propertyID) |>
    mutate(delta_year = c(first_flag, diff(year)),
           delta_take = c(first_flag, diff(all_take)),
           delta_events = c(first_flag, diff(all_events)),
           avg_take_density_in_pp = avg_take_in_pp / property_area_km2,
           avg_event_density_in_pp = avg_events_in_pp / property_area_km2,
           events_per_pp = all_events / n_sampled_pp) |>
    ungroup()

}

create_ml_data <- function(df, per_year_df, land_cover_df){
  df |>
    left_join(per_year_df) |>
    left_join(land_cover_df) |>
    mutate(
      y = med_density,
      year_fac = as.factor(year),
      farm_bill = as.factor(farm_bill),
      propertyID = factor(propertyID),
      st_name = factor(st_name),
      county_code = factor(county_code),
      state_year = factor(paste(st_name, year_fac)),
      county_year = factor(paste(county_code, year_fac)),
      property_year = factor(paste(propertyID, year_fac))) |>
    select(-med_density, -agrp_prp_id)
}

get_single_year_properties <- function(df){
  df |>
    group_by(propertyID) |>
    count() |>
    ungroup() |>
    filter(n == 1) |>
    pull(propertyID) |>
    unique()
}

get_multi_year_properties <- function(df){
  df |>
    filter(!propertyID %in% single_year_properties) |>
    pull(propertyID) |>
    unique()
}

create_density_m1 <- function(df, single_property_vec, multi_property_vec){

  tmp1 <- df |>
    filter(propertyID %in% multi_property_vec) |>
    group_by(propertyID) |>
    mutate(density_m1 = c(NA, y[1:(n()-1)])) |>
    ungroup()

  tmp2 <- df |>
    filter(propertyID %in% single_property_vec) |>
    mutate(density_m1 = NA)

  bind_rows(tmp1, tmp2) |>
    mutate(rowID = 1:n())


}

