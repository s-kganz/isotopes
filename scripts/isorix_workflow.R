library(mgcv)
library(MuMIn)
library(tidyverse)
library(IsoriX)
library(lubridate)
library(janitor)
library(DHARMa)
library(terra)

# Load data ----
trees <- read_csv("data_in/Tree_ds_sep16.csv") %>%
  clean_names() %>%
  # Drop observations on west coast
  filter(longitude > -100) %>%
  # Parse datetimes
  mutate(date_collected = parse_date_time(date_collected, orders=c("dmy", "mdy")))

# Convert to isorix format ----
trees_iso <- trees %>%
  transmute(
    source_ID = factor(paste("site", latitude, longitude, elevation_dem, sep = "_")),
    lat=latitude,
    long=longitude,
    elev=elevation_dem,
    year=year(date_collected),
    month=month(date_collected),
    source_value=d18o
  ) %>%
  # shuffle for split
  sample_frac()

# Split into data we will use for the isoscape and data we will use to
# test geographic assignment.
trees_isoscape <- head(trees_iso, 50) %>% prepsources()
trees_assign   <- tail(trees_iso, nrow(trees_iso)-50)

# Fit model ----
tree_fit <- trees_isoscape %>%
  isofit(
    mean_model_fix=list(elev=TRUE, lat=TRUE)
  )

# Get R2 and RMSE
y_hat <- predict(tree_fit$mean_fit)
y     <- tree_fit$mean_fit$y

r2 <- cor(y, y_hat)^2
rmse <- sqrt(mean((y - y_hat)^2))

# TODO model diagnostic tests. These do not look great.
sims <- simulateResiduals(tree_fit$mean_fit)
plot(sims)

# Build isoscape ----
# TODO create inference polygon from mean distance b/t points
getelev(
  file="data_in/elevation_world.tif",
  z=4,
  long_min=min(trees$longitude),
  long_max=max(trees$longitude),
  lat_min=min(trees$latitude),
  lat_max=max(trees$latitude),
  overwrite=TRUE
)

elev <- rast("data_in/elevation_world.tif")
plot(elev)

elev_prep <- prepraster(
  raster = elev,
  isofit = tree_fit,
  aggregation_factor = 1
)

plot(elev_prep)

tree_isoscape <- isoscape(
  raster = elev_prep,
  isofit = tree_fit
)

plot(tree_isoscape)

# Calibration model ----
# This uses the sample we held out from the isoscape
trees_assign <- trees_assign %>%
  rename(sample_value=source_value,
         site_ID=source_ID)
tree_calib <- calibfit(trees_assign, isofit=tree_fit)

plot(tree_calib)

# Infer origins ----
tree_locs <- isofind(
  data = trees_assign %>% mutate(sample_ID=paste("sample", row_number(), sep="_")),
  isoscape = tree_isoscape,
  calibfit = tree_calib
)

geom_vect <- vect(trees_assign, geom=c("long", "lat"), crs="+proj=longlat") %>%
  project(crs(tree_locs$sample$stat$sample_1))

plot_sample <- function(sample_num) {
  plot(tree_locs, who=sample_num,
       sources = list(draw = FALSE),
       calibs = list(draw = FALSE),
       assigns = list(draw = FALSE)
  )
  #points(geom_vect[sample_num], col="red", cex=5)
}

plot_sample(2)
