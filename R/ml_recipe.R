my_recipe <- function(data){
  require(rsample)
  require(recipes)

  blueprint <- recipe(y ~ ., data = data) |>

    # recommended order from recipes package

    # Impute (does not apply)
    # Handle factor levels
    step_novel(all_nominal_predictors()) |>

    # Individual transformations for skewness and other issues
    step_YeoJohnson(all_outcomes()) |>
    step_YeoJohnson(all_numeric_predictors()) |>

    # Discretize (if needed and if you have no other choice; does not apply)
    # Create dummy variables
    step_dummy(all_nominal_predictors()) |>

    # Create interactions
    step_interact(terms = ~delta_density:delta_year) |>

    step_interact(terms = ~all_take:all_events) |>
    step_interact(terms = ~all_take:all_units) |>
    step_interact(terms = ~all_events:all_units) |>

    step_interact(terms = ~all_take_density:all_events_density) |>
    step_interact(terms = ~all_take_density:all_units_density) |>
    step_interact(terms = ~all_events_density:all_units_density) |>

    step_interact(terms = ~avg_take_in_pp:avg_events_in_pp) |>
    step_interact(terms = ~avg_take_density_in_pp:avg_event_density_in_pp) |>

    step_interact(terms = ~all_take:rural.road.density) |>
    step_interact(terms = ~all_take:prop.pub.land) |>
    step_interact(terms = ~all_take:mean.ruggedness) |>
    step_interact(terms = ~all_take:mean.canopy.density) |>

    step_interact(terms = ~avg_take_in_pp:rural.road.density) |>
    step_interact(terms = ~avg_take_in_pp:prop.pub.land) |>
    step_interact(terms = ~avg_take_in_pp:mean.ruggedness) |>
    step_interact(terms = ~avg_take_in_pp:mean.canopy.density) |>

    step_interact(terms = ~all_take_density:rural.road.density) |>
    step_interact(terms = ~all_take_density:prop.pub.land) |>
    step_interact(terms = ~all_take_density:mean.ruggedness) |>
    step_interact(terms = ~all_take_density:mean.canopy.density) |>

    step_interact(terms = ~avg_take_density_in_pp:rural.road.density) |>
    step_interact(terms = ~avg_take_density_in_pp:prop.pub.land) |>
    step_interact(terms = ~avg_take_density_in_pp:mean.ruggedness) |>
    step_interact(terms = ~avg_take_density_in_pp:mean.canopy.density) |>

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
