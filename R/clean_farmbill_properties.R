library(dplyr)
library(readr)

data_farm_bill <- read_csv("../pigs-statistical/data/All_FB_Agreements_052124.csv")

# there are lots of empty rows and columns, drop them
filter_nas <- data_farm_bill |>
  filter(if_any(everything(), ~ !is.na(.))) |>
  select(where(~ !all(is.na(.)))) |>
  distinct() |>
  rename(propertyID = `PRP ID`)

dat <- filter_nas |>
  select(STATE, `AGR ID`, AGREEMENT, propertyID) |>
  mutate(property_name = NA)

row <- filter_nas |>
  filter(propertyID == 390782)

extra_properties <- tibble()
pb <- txtProgressBar(max = nrow(filter_nas), style = 3)
for(i in 1:nrow(filter_nas)){
  row <- filter_nas |> slice(i)
  row_data <- dat |> slice(i)

  skip <- FALSE
  for(j in 6:ncol(row)){

    if(skip){
      skip <- FALSE
      next
    }

    cell <- row[j]

    if(grepl("\\d{5}", cell) | grepl("\\d{6}", cell)){
      skip <- TRUE
      dat_j <- row[,c(j, j + 1)]
      colnames(dat_j) <- c("propertyID", "property_name")

      if(!grepl("\\.", dat_j$propertyID)){
        row_data$propertyID <- dat_j$propertyID
        row_data$property_name <- dat_j$property_name
        extra_properties <- bind_rows(extra_properties, row_data)
      }
    }
  }
  setTxtProgressBar(pb, i)
}
close(pb)


fb <- bind_rows(dat, extra_properties) |>
  mutate(property_name = if_else(is.na(property_name), AGREEMENT, property_name)) |>
  rename(agreementID = `AGR ID`,
         agreement_name = AGREEMENT)

write_csv(fb, "../pigs-statistical/data/All_FB_Agreements_long_2024-05-30.csv")
