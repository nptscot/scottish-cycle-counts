## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  message=FALSE,
  warning=FALSE
)


## -----------------------------------------------------------------------------
library(tidyverse)


## -----------------------------------------------------------------------------
zipped_data = list.files(pattern = ".zip")
zipped_data


## ---- echo=FALSE, eval=FALSE--------------------------------------------------
#> tail(zipped_data, 1)


## -----------------------------------------------------------------------------
unzip(zipped_data, exdir = "data-raw")
files_csv = list.files("data-raw", pattern = ".csv", full.names = TRUE)
files_csv


## -----------------------------------------------------------------------------
# counts = arrow::open_dataset("data-raw")
counts = map_dfr(files_csv, read_csv)
dim(counts)
counts


## -----------------------------------------------------------------------------
counts_monthly = counts |>
  mutate(
    year = year(endTime),
    month = month(endTime, label = TRUE),
    # date rounded to nearest month in 2020-01-01 format:
    date = lubridate::floor_date(endTime, unit = "month")
  ) |>
  group_by(date, area) |>
  summarise(
    count = sum(count)
  )
# Add column with names for most common areas:
# Most common areas are:
area_counts = counts_monthly |>
  group_by(area) |>
  summarise(
    count = sum(count)
  ) |>
  arrange(desc(count))
top_5_areas = head(area_counts, n = 5)
# Add column that is area name if in top 5, else "Other":
counts_monthly_top = counts_monthly |>
  mutate(
    Area = ifelse(
      area %in% top_5_areas$area,
      area,
      "Other"
    )
  ) |>
  group_by(date, Area) |>
  summarise(
    count = sum(count)
  )


## -----------------------------------------------------------------------------
counts_monthly_top |>
  ggplot(aes(x = date, y = count, colour = Area)) +
  geom_line() +
  # Add log y-axis:
  scale_y_log10()


## ---- echo=FALSE--------------------------------------------------------------
# Convert README.Rmd to counts.R:  
knitr::purl("README.Rmd", "counts.R")

