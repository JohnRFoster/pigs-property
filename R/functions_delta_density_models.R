# delta density model functions

model_message <- function(m, dest){
  print(k.check(m))
  print(summary(m))
  write_rds(m, dest)
  message("\n\n")
}

# Groups that could matter
# state
# county
# year
# state-year (see above)
# county-year (see above)

# create models by group/smoother and variable criteria

# T models - yearly totals

# M models - yearly means

# D models - yearly change

m0 <- function(d, m, f, dest, dul){
  mod <- gam(delta_density ~ 1, data = d, method = m, family = f, drop.unused.levels = dul)
  model_message(mod, dest)
}

m1 <- function(d, m, f, dest, dul){
  mod <- gam(delta_density ~ state_year, data = d, method = m, family = f, drop.unused.levels = dul)
  model_message(mod, dest)
}

# G model - A single common smoother for all observations; only has a Global smoother

G_T <- function(d, m, f, k_sum_take, k_state_year, dest, dul){
  mod <- gam(delta_density ~
        s(sum_take, k = k_sum_take, bs = "tp") +
        s(sum_events, bs = "tp") +
        s(property_area_km2, bs = "tp") +
        s(c_road_den, bs = "tp") +
        s(c_rugged, bs = "tp") +
        s(c_canopy, bs = "tp") +
        s(c_prop.pub.land, bs = "tp") +
        s(state_year, bs = "re", k = k_state_year),
      data = d,
      method = m,
      family = f, drop.unused.levels = dul)
  model_message(mod, dest)
}

G_M <- function(d, m, f, k_sum_take, k_state_year, dest, dul){
  mod <- gam(delta_density ~
        s(avg_take, k = k_sum_take, bs = "tp") +
        s(avg_events, bs = "tp") +
        s(property_area_km2, bs = "tp") +
        s(c_road_den, bs = "tp") +
        s(c_rugged, bs = "tp") +
        s(c_canopy, bs = "tp") +
        s(c_prop.pub.land, bs = "tp") +
        s(state_year, bs = "re", k = k_state_year),
      data = d,
      method = m,
      family = f, drop.unused.levels = dul)
  model_message(mod, dest)
}

G_D <- function(d, m, f, k_sum_take, k_state_year, dest, dul){
  mod <- gam(delta_density ~
        s(delta_take, k = k_sum_take, bs = "tp") +
        s(delta_events, bs = "tp") +
        s(property_area_km2, bs = "tp") +
        s(c_road_den, bs = "tp") +
        s(c_rugged, bs = "tp") +
        s(c_canopy, bs = "tp") +
        s(c_prop.pub.land, bs = "tp") +
        s(state_year, bs = "re", k = k_state_year),
      data = d,
      method = m,
      family = f, drop.unused.levels = dul)
  model_message(mod, dest)
}

# GS model - A global smoother plus group-level smoothers that have the *same* wiggliness
#   - Global smoother with individual effects that have a *Shared* penalty

GS_T <- function(d, m, f, dest, dul){
  mod <- gam(delta_density ~
        s(sum_take, bs = "tp", m = 2) +
        s(sum_events, bs = "tp", m = 2) +
        s(property_area_km2, bs = "tp", m = 2) +
        s(c_road_den, bs = "tp", m = 2) +
        s(c_rugged, bs = "tp", m = 2) +
        s(c_canopy, bs = "tp", m = 2) +
        s(c_prop.pub.land, bs = "tp", m = 2) +
        s(sum_take, state_year, bs = "fs", m = 2) +
        s(sum_events, state_year, bs = "fs", m = 2),
      method = m,
      family = f,
      data = d, drop.unused.levels = dul)
  model_message(mod, dest)
}

GS_M <- function(d, m, f, dest, dul){
  mod <- gam(delta_density ~
        s(avg_take, bs = "tp", m = 2) +
        s(avg_events, bs = "tp", m = 2) +
        s(property_area_km2, bs = "tp", m = 2) +
        s(c_road_den, bs = "tp", m = 2) +
        s(c_rugged, bs = "tp", m = 2) +
        s(c_canopy, bs = "tp", m = 2) +
        s(c_prop.pub.land, bs = "tp", m = 2) +
        s(avg_take, state_year, bs = "fs", m = 2) +
        s(avg_events, state_year, bs = "fs", m = 2),
      method = m,
      family = f,
      data = d, drop.unused.levels = dul)
  model_message(mod, dest)
}

GS_D <- function(d, m, f, dest, dul){
  mod <- gam(delta_density ~
        s(delta_take, bs = "tp", m = 2) +
        s(delta_events, bs = "tp", m = 2) +
        s(property_area_km2, bs = "tp", m = 2) +
        s(c_road_den, bs = "tp", m = 2) +
        s(c_rugged, bs = "tp", m = 2) +
        s(c_canopy, bs = "tp", m = 2) +
        s(c_prop.pub.land, bs = "tp", m = 2) +
        s(delta_take, state_year, bs = "fs", m = 2) +
        s(delta_events, state_year, bs = "fs", m = 2),
      method = m,
      family = f,
      data = d, drop.unused.levels = dul)
  model_message(mod, dest)
}

# GI model - A global smoother plus group-level smoothers with *differing* wiggliness
#   - Global smoother with individual effects that have *Individual* penalties
#   - This is useful if different groups differ substantially in how wiggly they are.

GI_T <- function(d, m, f, k_state_year, dest, dul){
  mod <- gam(delta_density ~
        s(sum_take, bs = "tp", m = 2) +
        s(sum_events, bs = "tp", m = 2) +
        s(property_area_km2, bs = "tp", m = 2) +
        s(c_road_den, bs = "tp", m = 2) +
        s(c_rugged, bs = "tp", m = 2) +
        s(c_canopy, bs = "tp", m = 2) +
        s(c_prop.pub.land, bs = "tp", m = 2) +
        s(sum_take, by = state_year, bs = "tp", m = 1) +
        s(sum_events, by = state_year, bs = "tp", m = 1) +
        s(state_year, bs = "re", k = k_state_year),
      method = m,
      family = f,
      data = d, drop.unused.levels = dul)
  model_message(mod, dest)
}

GI_M <- function(d, m, f, k_state_year, dest, dul){
  mod <- gam(delta_density ~
        s(avg_take, bs = "tp", m = 2) +
        s(avg_events, bs = "tp", m = 2) +
        s(property_area_km2, bs = "tp", m = 2) +
        s(c_road_den, bs = "tp", m = 2) +
        s(c_rugged, bs = "tp", m = 2) +
        s(c_canopy, bs = "tp", m = 2) +
        s(c_prop.pub.land, bs = "tp", m = 2) +
        s(avg_take, by = state_year, bs = "tp", m = 1) +
        s(avg_events, by = state_year, bs = "tp", m = 1) +
        s(state_year, bs = "re", k = k_state_year),
      method = m,
      family = f,
      data = d, drop.unused.levels = dul)
  model_message(mod, dest)
}


GI_D <- function(d, m, f, k_state_year, dest, dul){
  mod <- gam(delta_density ~
        s(delta_take, bs = "tp", m = 2) +
        s(delta_events, bs = "tp", m = 2) +
        s(property_area_km2, bs = "tp", m = 2) +
        s(c_road_den, bs = "tp", m = 2) +
        s(c_rugged, bs = "tp", m = 2) +
        s(c_canopy, bs = "tp", m = 2) +
        s(c_prop.pub.land, bs = "tp", m = 2) +
        s(delta_take, by = state_year, bs = "tp", m = 1) +
        s(delta_events, by = state_year, bs = "tp", m = 1) +
        s(state_year, bs = "re", k = k_state_year),
      method = m,
      family = f,
      data = d, drop.unused.levels = dul)
  model_message(mod, dest)
}

# S model - Group-specific smoothers without a global smoother
#   - all smoothers having the *same* wiggliness

S_T <- function(d, m, f, dest, dul){
  mod <- gam(delta_density ~
        s(sum_take, state_year, bs = "fs", m = 2) +
        s(sum_events, state_year, bs = "fs", m = 2) +
        s(property_area_km2, bs = "tp", m = 2) +
        s(c_road_den, bs = "tp", m = 2) +
        s(c_rugged, bs = "tp", m = 2) +
        s(c_canopy, bs = "tp", m = 2) +
        s(c_prop.pub.land, bs = "tp", m = 2),
      method = m,
      family = f,
      data = d, drop.unused.levels = dul)
  model_message(mod, dest)
}

S_M <- function(d, m, f, dest, dul){
  mod <- gam(delta_density ~
        s(avg_take, state_year, bs = "fs", m = 2) +
        s(avg_events, state_year, bs = "fs", m = 2) +
        s(property_area_km2, bs = "tp", m = 2) +
        s(c_road_den, bs = "tp", m = 2) +
        s(c_rugged, bs = "tp", m = 2) +
        s(c_canopy, bs = "tp", m = 2) +
        s(c_prop.pub.land, bs = "tp", m = 2),
      method = m,
      family = f,
      data = d, drop.unused.levels = dul)
  model_message(mod, dest)
}

S_D <- function(d, m, f, dest, dul){
  mod <- gam(delta_density ~
        s(delta_take, state_year, bs = "fs", m = 2) +
        s(delta_events, state_year, bs = "fs", m = 2) +
        s(property_area_km2, bs = "tp", m = 2) +
        s(c_road_den, bs = "tp", m = 2) +
        s(c_rugged, bs = "tp", m = 2) +
        s(c_canopy, bs = "tp", m = 2) +
        s(c_prop.pub.land, bs = "tp", m = 2),
      method = m,
      family = f,
      data = d, drop.unused.levels = dul)
  model_message(mod, dest)
}


# I model - Group-specific smoothers with different wiggliness
#   - all smoothers having the *different* wiggliness

I_T <- function(d, m, f, k_state_year, dest, dul){
  mod <- gam(delta_density ~
        s(property_area_km2, bs = "tp", m = 2) +
        s(c_road_den, bs = "tp", m = 2) +
        s(c_rugged, bs = "tp", m = 2) +
        s(c_canopy, bs = "tp", m = 2) +
        s(c_prop.pub.land, bs = "tp", m = 2) +
        s(sum_take, by = state_year, bs = "tp", m = 2) +
        s(sum_events, by = state_year, bs = "tp", m = 2) +
        s(state_year, bs = "re", k = k_state_year),
      method = m,
      family = f,
      data = d, drop.unused.levels = dul)
  model_message(mod, dest)
}

I_M <- function(d, m, f, k_state_year, dest, dul){
  mod <- gam(delta_density ~
        s(property_area_km2, bs = "tp", m = 2) +
        s(c_road_den, bs = "tp", m = 2) +
        s(c_rugged, bs = "tp", m = 2) +
        s(c_canopy, bs = "tp", m = 2) +
        s(c_prop.pub.land, bs = "tp", m = 2) +
        s(avg_take, by = state_year, bs = "fs", m = 2) +
        s(avg_events, by = state_year, bs = "fs", m = 2) +
        s(state_year, bs = "re", k = k_state_year),
      method = m,
      family = f,
      data = d, drop.unused.levels = dul)
  model_message(mod, dest)
}


I_D <- function(d, m, f, k_state_year, dest, dul){
  mod <- gam(delta_density ~
        s(property_area_km2, bs = "tp", m = 2) +
        s(c_road_den, bs = "tp", m = 2) +
        s(c_rugged, bs = "tp", m = 2) +
        s(c_canopy, bs = "tp", m = 2) +
        s(c_prop.pub.land, bs = "tp", m = 2) +
        s(delta_take, by = state_year, bs = "fs", m = 2) +
        s(delta_events, by = state_year, bs = "fs", m = 2) +
        s(state_year, bs = "re", k = k_state_year),
      method = m,
      family = f,
      data = d, drop.unused.levels = dul)
  model_message(mod, dest)
}






