

get_next_property <- function(all_data, data_last_fit){
  require(rsample)

  # find properties not in first fit
  last_fit_properties <- unique(data_last_fit$propertyID)
  all_properties <- unique(all_data$propertyID)
  not_fit_properties <- setdiff(all_properties, last_fit_properties)
  targets::tar_assert_true(
    length(not_fit_properties) + length(last_fit_properties) == length(all_properties)
  )

  ts_length <- get_ts_length(all_data)
  n_obs <- get_n_observations(all_data, ts_length$propertyID)

  ts_sorted <- left_join(ts_length, n_obs) |>
    mutate(percent_obs_strata = n / delta) |>
    mutate(percent_obs_strata = as.numeric(as.factor(make_strata(percent_obs_strata, breaks = 10))))


  df_strata <- all_data |>
    mutate(canopy_strata = as.numeric(as.factor(make_strata(c_canopy, breaks = 10))),
           rugged_strata = as.numeric(as.factor(make_strata(c_rugged, breaks = 10))),
           road_den_strata = as.numeric(as.factor(make_strata(c_road_den, breaks = 10)))) |>
    select(propertyID, st_name, cnty_name, contains("strata")) |>
    distinct()

  # rank remaining properties
  #   proportion of timeseries that is sampled
  #     if tied, rank by length, short to long
  df_join <- left_join(ts_sorted, df_strata)

  varmax <- df_join |>
    group_by(canopy_strata, rugged_strata, road_den_strata) |>
    summarise(gmax = max(percent_obs_strata)) |>
    ungroup()

  # this is the order of properties if we didn't fit the initial batch
  property_order <- left_join(df_join, varmax) |>
    arrange(desc(gmax), desc(percent_obs_strata))

  # remaining properties
  remaining_property_order <- property_order |>
    filter(!propertyID %in% last_fit_properties)

  if(nrow(remaining_property_order) == 0){
    stop("All properties have been fit!")
  } else {
    next_property <- remaining_property_order |>
      slice(1) |>
      pull(propertyID)

    new_data <- all_data |> filter(propertyID == next_property)

    to_fit_properties <- unique(new_data$propertyID)
    not_fit_properties <- setdiff(all_properties, to_fit_properties)
    targets::tar_assert_true(
      length(not_fit_properties) + length(to_fit_properties) == length(all_properties)
    )

    return(bind_rows(data_last_fit, new_data))
  }

}
