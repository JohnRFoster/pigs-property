#######################################################
### helper functions and definitions for results.qmd

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(stringr)

method_vector <-  c("Firearms", "Fixed wing", "Helicopter", "Snare", "Trap")

method_lookup <- tibble(
  method_idx = 1:5,
  method = method_vector
)

land_lookup <- tibble(
  land_idx = 1:3,
  land = c("Road density", "Ruggedness", "Canopy cover")
)



# linear
potential_area <- function(params, m, metric){

  if(m == "Firearms") id <- 1
  if(m == "Fixed wing") id <- 2
  if(m == "Helicopter") id <- 3
  if(m == "Snare") id <- 4
  if(m == "Trap") id <- 5

  log_rho <- params |>
    pull(paste0("log_rho[", id, "]"))

  beta1 <- params |>
    pull(paste0("beta1[", id, "]"))

  min_e <- effort_range |> filter(method == m) |> pull(min_e)
  max_e <- effort_range |> filter(method == m) |> pull(max_e)
  log_effort <- seq(log(min_e), log(max_e), length.out = 50)

  log_potential_area <- tibble()

  if(id <= 3){
    for(i in seq_along(log_effort)){
      log_e <- log_effort[i]

      lpa <- tibble(
        log_rho = log_rho,
        beta1 = beta1,
        log_effort_per = log_e,
        lpa = log_rho + log_e
      )

      log_potential_area <- bind_rows(log_potential_area, lpa)

    }
  } else {

    log_gamma <- params |>
      pull(paste0("log_gamma[", id - 3, "]"))
    p_unique <- params |>
      pull(paste0("p_mu[", id - 3, "]")) |>
      boot::inv.logit()
    n_trap_m1 <- if_else(id == 4, 9, 1)

    for(i in seq_along(log_effort)){
      log_e <- log_effort[i]

      lpa1 <- log(pi) +
        (2 * (log_rho + log_e -
                log(exp(log_gamma) + exp(log_e)))) +
        log(1 + (p_unique * n_trap_m1))


      lpa <- tibble(
        log_effort_per = log_e,
        beta1 = beta1,
        lpa = lpa1
      )

      log_potential_area <- bind_rows(log_potential_area, lpa)

    }
  }

  if(metric == "potential area"){
    log_potential_area |>
      mutate(pa = exp(lpa),
             effort_per = exp(log_effort_per)) |>
      group_by(effort_per) |>
      summarise(`5%` = quantile(pa, 0.05),
                `25%` = quantile(pa, 0.25),
                `50%` = quantile(pa, 0.5),
                `75%` = quantile(pa, 0.75),
                `95%` = quantile(pa, 0.95)) |>
      ungroup() |>
      mutate(method = m,
             value = "Potential area")

  } else if (metric == "capture"){
    log_potential_area |>
      mutate(pa = exp(lpa),
             effort_per = exp(log_effort_per),
             beta1 = boot::inv.logit(beta1),
             pa = beta1 * pa) |>
      group_by(effort_per) |>
      summarise(`5%` = quantile(pa, 0.05),
                `25%` = quantile(pa, 0.25),
                `50%` = quantile(pa, 0.5),
                `75%` = quantile(pa, 0.75),
                `95%` = quantile(pa, 0.95)) |>
      ungroup() |>
      mutate(method = m,
             value = "Effective area")
  }

}



get_prior <- function(mu, tau, n = 1e6, n_out = 1000){
  prior <- rnorm(n, mu, 1 / sqrt(tau))
  seq(quantile(prior, 0.05), quantile(prior, 0.95), length.out = n_out)
}





plot_post <- function(dfg, xlab, title){

  dfg |>
    ggplot() +
    aes(x = `50%`, xmin = `5%`, xmax = `95%`, y = method, color = distribution) +
    geom_linerange(position = position_dodge(width = 0.5)) +
    # geom_linerange(aes(xmin = `25%`, xmax = `75%`),
    #                linewidth = 3,
    #                position = position_dodge(width = 0.5)) +
    geom_point(position = position_dodge(width = 0.5), size = 4) +
    scale_color_manual(values = dist_colors) +
    labs(y = "Method",
         x = xlab,
         title = title,
         color = element_blank()) +
    theme_bw() +
    my_theme()

}

get_posterior <- function(df, y){
  dfp <- df |>
    select(distribution, starts_with(y)) |>
    # mutate(distribution = "Posterior") |>
    pivot_longer(cols = -distribution,
                 names_to = "node",
                 values_to = "y")

  if(grepl("log", y)) dfp <- dfp |> mutate(y = exp(y))
  if(y == "log_nu") dfp <- dfp |> mutate(y = y / 2)
  if(grepl("p_mu", y) | grepl("beta1", y)) dfp <- dfp |> mutate(y = boot::inv.logit(y))

  return(dfp)
}

join_summarise_methods <- function(df, df_prior, df_method_names){
  df |>
    bind_rows(df_prior) |>
    left_join(df_method_names) |>
    group_by(distribution, method) |>
    summarise(`5%` = quantile(y, 0.05),
              `25%` = quantile(y, 0.25),
              `50%` = quantile(y, 0.5),
              `75%` = quantile(y, 0.75),
              `95%` = quantile(y, 0.95)) |>
    ungroup() |>
    mutate(method = factor(method, levels = c("Prior", method_vector)))
}







#
# config <- config::get(config = "default")
# data_repo <- config$data_repo
#
# ## MIS data ----
# file <- file.path(data_repo, config$file_mis)
# interval <- config$interval
#
# all_take <- read_csv(file, show_col_types = FALSE) |>
#   mutate(cnty_name = if_else(grepl("ST ", cnty_name), gsub("ST ", "ST. ", cnty_name), cnty_name),
#          cnty_name = if_else(grepl("KERN", cnty_name), "KERN", cnty_name))
#
# data_mis <- all_take |>
#   mutate(property_area_km2 = property.size / 247.1) |>
#   filter(property_area_km2 >= 1.8,
#          st_name != "HAWAII") |>
#   mutate(effort = if_else(cmp_name %in% c("TRAPS, CAGE", "SNARE"), cmp.days, cmp.hours),
#          effort_per = effort / cmp.qty,
#          cmp_name = if_else(cmp_name == "TRAPS, CAGE", "TRAPS", cmp_name)) |>
#   rename(method = cmp_name,
#          trap_count = cmp.qty) |>
#   select(-wt_work_date, -hours, -cmp.hours, -cmp.days) |>
#   distinct()
#
# end_dates <- unique(sort(data_mis$end.date))
# min_date <- min(end_dates)
# max_date <- max(end_dates)
#
# start_dates <- seq(min_date, max_date, by = paste(interval, "week"))
# end_dates <- c(start_dates[-1] - 1, max_date)
#
# targets::tar_assert_identical(length(start_dates), length(end_dates))
# targets::tar_assert_true(min(data_mis$start.date) >= min_date)
# targets::tar_assert_true(max(data_mis$start.date) <= max_date)
#
# timestep_df <- tibble(start_dates, end_dates) |>
#   mutate(primary_period = 1:n())
# timestep_df$month <- month(timestep_df$end_dates)
# timestep_df$year <- year(timestep_df$end_dates)
#
# source("R/functions_data.R")
# data_timestep <- create_primary_periods(data_mis, interval) |>
#   resolve_duplicate() |>         # resolve duplicate property areas
#   take_filter() |>               # remove properties with zero pigs taken
#   dynamic_filter() |>            # filter out bad events & properties
#   condition_first_capture() |>   # condition on first positive removal event for each property
#   dynamic_filter() |>
#   rename(primary_period = timestep)
#
# timestep_df_props <- data_timestep |>
#   create_timestep_df() |>
#   left_join(timestep_df)
#
# data_dates <- left_join(data_timestep, timestep_df_props,
#                      by = join_by(agrp_prp_id, primary_period)) |>
#   mutate(primary_period = primary_period - min(primary_period) + 1) |>
#   select(agrp_prp_id, primary_period, end_dates)
#
# posterior_path <- "C:/Users/John.Foster/Downloads/iterativeFits/1_posterior/densitySummaries.rds"
# density_summary <- read_rds(posterior_path)
#
#
# all_pp <- create_all_primary_periods(data) |> select(-timestep)
#
# data_pp <- left_join(data, all_pp)
#
# data_join <- data_pp |>
#   group_by(st_name, cnty_name, property, property_area_km2, primary_period) |>
#   summarise(take = sum(take)) |>
#   ungroup() |>
#   mutate(take_density = take / property_area_km2) |>
#   left_join(density_summary) |>
#   left_join(data_dates) |>
#   distinct() |>
#   select(-take, -node, -n_id, -agrp_prp_id, -`0.025`, -`0.1`, -`0.9`, -`0.975`)
#
# stat_by_window <- function(i, dates, window, value, stat){
#   startdate <- dates[i]-window
#   enddate <- dates[i]
#   interval <- seq(startdate, enddate, 1)
#
#   tmp <- value[dates %in% interval]
#   if(stat == "average") return(mean(tmp))
#   if(stat == "sd") return(sd(tmp))
# }
#
# apply_sbw <- function(dates, window, value){
#   avg <- sapply(1:length(dates), function(x) stat_by_window(x, dates, window, value, stat = "average"))
#   avg[is.nan(avg)] <- NA
#   sd <- sapply(1:length(dates), function(x) stat_by_window(x, dates, window, value, stat = "sd"))
#   sd[is.nan(sd)] <- NA
#
#   return(data.frame(avg = avg, sd = sd))
#
# }
#
# properties <- unique(data_join$property)
#
# # standardized density index
# data_standard_index <- tibble()
# for(i in seq_along(properties)){
#   property_df <- data_join |>
#     filter(property == properties[i])
#
#   end_dates <- unique(property_df$end_dates)
#   values <- property_df$`0.5`
#
#   property_df$year0.5_avg <- apply_sbw(end_dates, 183, values)$avg
#   property_df$year0.5_sd <- apply_sbw(end_dates, 183, values)$sd
#   property_df$year1_avg <- apply_sbw(end_dates, 365, values)$avg
#   property_df$year1_sd <- apply_sbw(end_dates, 365, values)$sd
#   property_df$year5_avg <- apply_sbw(end_dates, 365*5, values)$avg
#   property_df$year5_sd <- apply_sbw(end_dates, 365*5, values)$sd
#
#   property_df <- property_df |>
#     mutate(sdi_year0.5 = (`0.5` - year0.5_avg) / year0.5_sd,
#            sdi_year1 = (`0.5` - year1_avg) / year1_sd,
#            sdi_year5 = (`0.5` - year5_avg) / year5_sd)
#
#   data_standard_index <- bind_rows(data_standard_index, property_df)
#
# }
#
# data |> group_by(property) |> summarise(effort = sum(take)) |> arrange((effort))
#
# property_filter <- 7
# property_df <- data_standard_index |> filter(property == property_filter)
#
# county <- property_df$cnty_name[1]
# state <- property_df$st_name[1]
# area <- round(property_df$property_area_km2[1], 2)
#
# color <- "darkgreen"
#
# my_theme <- function(size = 12){
#   theme(
#     title = element_text(size = size + 4),
#     axis.title = element_text(size = size + 2),
#     axis.text = element_text(size = size)
#   )
# }
#
# ###
# property_df |>
#   ggplot() +
#   aes(x = primary_period, y = `0.5`) +
#   geom_linerange(aes(ymin = `0.05`, ymax = `0.95`), color = color, alpha = 0.2, linewidth = 3) +
#   geom_linerange(aes(ymin = `0.25`, ymax = `0.75`), color = color, alpha = 0.4, linewidth = 3) +
#   geom_point(aes(shape = "Median density"), color = color, alpha = 1, size = 3) +
#   geom_point(aes(y = take_density, shape = "Total take"), size = 3) +
#   labs(x = "Primary period",
#        y = "Density (pigs / sq. km)",
#        subtitle = paste0(area, " sq. km property in ", county, " county, ", state),
#        title = "Absolute density",
#        shape = element_blank()) +
#   theme_bw() +
#   my_theme()
#
# property_df |>
#   rename(`6 month` = sdi_year0.5,
#          `1 year` = sdi_year1,
#          `5 year` = sdi_year5) |>
#   pivot_longer(cols = c(`6 month`, `1 year`, `5 year`),
#                names_to = "sdi_window",
#                values_to = "sdi") |>
#   mutate(sdi_window = factor(sdi_window, levels = c("6 month", "1 year", "5 year"))) |>
#   ggplot() +
#   aes(x = end_dates, y = sdi, color = sdi_window) +
#   geom_line() +
#   geom_point() +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   labs(x = "Date",
#        y = "Standardized Density Index",
#        subtitle = paste0(area, " sq. km property in ", county, " county, ", state),
#        title = "Standardized Density Index",
#        color = "Window") +
#   theme_bw() +
#   my_theme()
#
#
# property_df |>
#   mutate(mean_density = mean(`0.5`),
#          `0.05` = `0.05` - mean_density,
#          `0.25` = `0.25` - mean_density,
#          `0.5` = `0.5` - mean_density,
#          `0.75` = `0.75` - mean_density,
#          `0.95` = `0.95` - mean_density,
#          take_density = take_density - mean_density) |>
#   ggplot() +
#   aes(x = primary_period, y = `0.5`) +
#   geom_linerange(aes(ymin = `0.05`, ymax = `0.95`), color = color, alpha = 0.2, linewidth = 3) +
#   geom_linerange(aes(ymin = `0.25`, ymax = `0.75`), color = color, alpha = 0.4, linewidth = 3) +
#   geom_point(aes(shape = "Median density"), color = color, alpha = 1, size = 3) +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   labs(x = "Primary period",
#        y = "Density anomaly (pigs / sq. km)",
#        subtitle = paste0(area, " sq. km property in ", county, " county, ", state),
#        title = "Density relative to historical average",
#        shape = element_blank()) +
#   theme_bw() +
#   my_theme()
#
# property_df |>
#   mutate(mean_density = mean(`0.5`),
#          `0.05` = `0.05` - mean_density,
#          `0.25` = `0.25` - mean_density,
#          `0.5` = `0.5` - mean_density,
#          `0.75` = `0.75` - mean_density,
#          `0.95` = `0.95` - mean_density,
#          rel_change = `0.5` / mean_density * 100) |>
#   ggplot() +
#   aes(x = primary_period, y = `0.5`) +
#   geom_linerange(aes(ymin = `0.05`, ymax = `0.95`, color = rel_change), alpha = 0.2, linewidth = 3) +
#   geom_linerange(aes(ymin = `0.25`, ymax = `0.75`, color = rel_change), alpha = 0.4, linewidth = 3) +
#   geom_point(aes(shape = "Median density", color = rel_change), alpha = 1, size = 3) +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   labs(x = "Primary period",
#        y = "Density anomaly (pigs / sq. km)",
#        subtitle = paste0(area, " sq. km property in ", county, " county, ", state),
#        title = "Density relative to historical average",
#        shape = element_blank(),
#        color = "Percent difference from\nhistorical average") +
#   scale_color_gradient(low = "blue", high = "red") +
#   theme_bw() +
#   my_theme()
#
# property_df |>
#   mutate(max_density = max(`0.95`),
#          `0.05` = `0.05` / max_density,
#          `0.25` = `0.25` / max_density,
#          `0.5` = `0.5` / max_density,
#          `0.75` = `0.75` / max_density,
#          `0.95` = `0.95` / max_density) |>
#   ggplot() +
#   aes(x = primary_period, y = `0.5`) +
#   geom_linerange(aes(ymin = `0.05`, ymax = `0.95`), color = color, alpha = 0.2, linewidth = 3) +
#   geom_linerange(aes(ymin = `0.25`, ymax = `0.75`), color = color, alpha = 0.4, linewidth = 3) +
#   geom_point(aes(shape = "Median density"), color = color, alpha = 1, size = 3) +
#   labs(x = "Primary period",
#        y = "Relative density",
#        subtitle = paste0(area, " sq. km property in ", county, " county, ", state),
#        title = "Relative density",
#        shape = element_blank()) +
#   theme_bw() +
#   my_theme()
#
# property_df |>
#   mutate(`0.5` = c(0, diff(`0.5`))) |>
#   ggplot() +
#   aes(x = primary_period, y = `0.5`) +
#   geom_point(aes(shape = "Median density"), color = color, alpha = 1, size = 3) +
#   geom_line(color = color, alpha = 1) +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   labs(x = "Primary period",
#        y = "Density delta (pigs / sq. km)",
#        subtitle = paste0(area, " sq. km property in ", county, " county, ", state),
#        title = "Absolute change in density",
#        shape = element_blank()) +
#   theme_bw() +
#   my_theme()
#
# property_df |>
#   mutate(`0.5` = c(0, diff(`0.5`)) / `0.5`) |>
#   ggplot() +
#   aes(x = primary_period, y = `0.5`) +
#   geom_point(aes(shape = "Median density"), color = color, alpha = 1, size = 3) +
#   geom_line(color = color, alpha = 1) +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   labs(x = "Primary period",
#        y = "Percent delta",
#        subtitle = paste0(area, " sq. km property in ", county, " county, ", state),
#        title = "Percent change in density",
#        shape = element_blank()) +
#   theme_bw() +
#   my_theme()



