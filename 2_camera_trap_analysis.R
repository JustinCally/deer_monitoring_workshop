# Analysing Camera Data ---------------------------------------------------------------------


#### Step 1: Setup ####
# install.packages("tidyverse")
library(tidyverse)
# install.packages("camtrapR")
library(camtrapR)
# note you will need to install exiftoolr and make camtrapR recognise it
# eg:
# camtrapR::addToPath()
# install.packages("mgcv")
library(mgcv)
# install.packages("mapview")
library(mapview)
# install.packages("Distance")
library(Distance)
# install.packages("dsm")
library(dsm)
# install.packages("kableExtra")
library(kableExtra)
# install.packages("terra")
library(terra)

#### Step 2: Extract metadata from images ####
cam_records_example <- recordTable(inDir = "example_data/camera_trap/raw_photos",
                                   IDfrom = "metadata",
                                   metadataSpeciesTag = "Species")

#### Step 3: Summarising and processing data ####
camtrap_records_sambar <- read_csv("example_data/camera_trap/camtrap_records_sambar.csv") %>%
  mutate(SiteID = as.character(SiteID),
         DateTime = as.POSIXct(paste(Date, Time)))

camera_deployments <- read_csv("example_data/camera_trap/camera_deployments.csv")

# get the total encounters per day/week for each camera
camtrap_records_delta <- camtrap_records_sambar %>%
  group_by(SiteID) %>%
  arrange(DateTime) %>%
  mutate(delta = DateTime - lag(DateTime)) %>% # get lag between photos
  filter(delta > lubridate::minutes(10))  # filter out records with 10 min of previous record

# per site calculate encounters
camtrap_records_delta_summary <- camtrap_records_delta %>%
  summarise(n_encounters = n()) %>%
  ungroup()

# Format per week using deployment data
camera_trap_rai <- camera_deployments %>%
  left_join(camtrap_records_delta_summary) %>%
  mutate(n_encounters = coalesce(n_encounters, 0), # if no encounters make it zero
         StartTime = DateTimeDeploy, # when camera is deployed
         EndTime = coalesce(Problem1_from, DateTimeRetrieve),
         DeploymentTime = as.numeric(EndTime - StartTime)/7,
         RAI = n_encounters/DeploymentTime) # when camera is retrieved/stops working

# Summary of RAI: lots of zeros for statewide surveys
summary(camera_trap_rai$RAI)

#### Step 4: Model RAI ####
# Simple X-Y spline to model relative abundance
fit_rai <- gam(RAI ~ s(scale(Longitude), scale(Latitude)),
               family = nb(link = "log"),
               data = camera_trap_rai)

# get summary and plot spline effect
summary(fit_rai)
plot(fit_rai, scheme=2, hcolors=heat.colors(999, rev =T), rug=F)

# generate predictions
rai_predictions <- predict(fit_rai, se.fit = T, type = "response") %>%
  as.data.frame()

# combine the predictions with our site data and make a little interactive map
map_data_rai <- camera_trap_rai %>%
  select(SiteID, Latitude, Longitude, RAI) %>%
  bind_cols(rai_predictions)

str(map_data_rai)

mapview(x = map_data_rai,
        xcol = "Longitude",
        ycol = "Latitude",
        zcol = "fit",
        grid = FALSE,
        crs = 4283)

#### Step 5: Camera trap distance sampling: distance function ####
# Tidy format for the records and deployments
snapshot_interval <- 2 # this basically filters down the data to every 2nd second
theta <- 40 * pi / 180 # camera angle in radians

# Based on the camera field of view and the deployment duration you can calculate a surey effort
camera_deployments_effort <- camera_deployments %>%
  dplyr::mutate(Tk = as.numeric(DateTimeRetrieve - DateTimeDeploy, units = "secs"), #seconds
                Tk_prob = dplyr::coalesce(as.numeric(as.POSIXct(Problem1_to,
                                                                format = "%Y-%m-%d %H:%M:%OS") -
                                                       as.POSIXct(Problem1_from, format = "%Y-%m-%d %H:%M:%OS"),
                                                     units = "secs"), 0), # How much time was it when cameras were not functioning
                Tk_adj = Tk-Tk_prob, # Remove the time cameras were not functioning
                Tkt = Tk_adj / snapshot_interval, # snapshot moments: every 2nd second
                Effort = (Tk_adj*theta)/(snapshot_interval * pi * 2)) %>% # This is basically the survey effort/area that is expended
  dplyr::arrange(SiteID) # Make sure sites are ordered

# combine the records with site information
dcount<- camtrap_records_sambar %>%
  dplyr::select(Species, SiteID, Distance, size, Date, Time, DateTime) %>%
  dplyr::full_join(camera_deployments_effort %>%
                     dplyr::select(SiteID, Tkt, Effort), by="SiteID") %>%
  dplyr::mutate(object=1:nrow(.)) %>% # ID for each observation
  dplyr::mutate(size = if_else(is.na(size),0L, size)) %>% # if size is NA than size = 0 (size is group size)
  dplyr::arrange(SiteID) # Make sure sites are ordered

# get the midpoints of each distance bin
summarised_count <- dcount %>%
  mutate(distance = case_when(Distance == "0 - 2.5" ~ 1.25,
                              Distance == "2.5 - 5" ~ 3.75,
                              Distance == "5 - 7.5" ~ 6.25,
                              Distance == "7.5 - 10" ~ 8.75,
                              Distance == "10+" ~ 11.25),
         Sample.Label = SiteID)

# setup parameters for distance sampling model
mybreaks <- c(0,2.5,5,7.5,10,12.5)
trunc.list <- list(left=0, right=12.5) # no left truncation but truncate at 12.5m
conversion <- convert_units("metre", "metre", "metre")

# Different models to test
# half-normal
hn0 <- ds(summarised_count, transect = "point", key="hn", adjustment = NULL,
          convert_units = conversion, cutpoints = mybreaks, truncation = trunc.list)

plot(hn0)

# hazard-rate model
hr0 <- ds(summarised_count, transect = "point", key="hr", adjustment = NULL,
          convert_units = conversion, cutpoints = mybreaks, truncation = 12.5)

plot(hr0)

# half-normal with adjustment
hn1 <- ds(summarised_count, transect = "point", key="hn", adjustment = "cos", order = 2,
          convert_units = conversion, cutpoints = mybreaks, truncation = trunc.list)

plot(hn1)

# hazard with adjustment
hr1 <- ds(summarised_count, transect = "point", key="hr", adjustment = "cos", order = 2,
          convert_units = conversion, cutpoints = mybreaks, truncation = trunc.list)

plot(hr1)

# visualise the model selction table
kableExtra::kbl(summarize_ds_models(hn0, hr0, hn1, hr1, output = "latex") %>%
                  mutate(Model = stringr::str_extract(Model, pattern = '(?<=\\{)[^\\}]+')),
                digits=3,
                caption="Model selection table for deer detection", format = "html") %>%
  kable_styling(bootstrap_options = "striped")

# plot the detection function
par(mfrow=c(1,2))
plot(hn1, showpoints=FALSE, lwd=2,xlab="Distance from camera (m)",pl.col="blue",
     pl.den=20,pdf=TRUE, main="(a) Observed distances")
plot(hn1, showpoints=FALSE, lwd=2,ylim=c(0,1),xlab="Distance from camera (m)",pl.col="blue",
     pl.den=20,pdf=FALSE, main ="(b) Detection function")

dev.off()

# see the goodness of fit
ddf.gof(hn1$ddf)

#### Step 6: Camera trap distance sampling: availability function ####
trigger.events <- dcount %>%
  filter(size != 0) %>% ### remove sites with zero counts
  dplyr::mutate(Sample.Label = as.character(SiteID)) # Distance requires a column of Sample.Label

# Make time column in radians
trigger.events$rtime <- activity::gettime(x = trigger.events$DateTime,
                                          tryFormats = "%Y-%m-%d %H:%M:%OS")

# Get availability
act_result <- activity::fitact(trigger.events$rtime,
                               sample="model",
                               reps=50, # bootstrap reps
                               bw = 30) # bandwidth
# Plot over 24 hours
plot(act_result)

# availability rate
avail <- list(creation=data.frame(rate = act_result@act[1],
                                  SE   = act_result@act[2]))

#### Step 7: Camera trap distance sampling: density surface model ####
# Create data for the density-surface model
segdata <- left_join(camera_deployments, summarised_count) %>% # join site info to the summarised counts
  mutate(Effort = Effort*avail$creation$rate) %>% # Multiple the Effort (Area x time deployed) by the availability for a new 'Effort'
  select(Sample.Label, Effort, X = Longitude, Y = Latitude) %>%
  distinct() %>% # Avoid duplication of site info
  sf::st_as_sf(coords = c("X", "Y"), crs = 4283) %>% # convert to projected geographics: epsg: 3111
  sf::st_transform(3111) %>%
  bind_cols(sf::st_coordinates(.)) %>%
  sf::st_drop_geometry()

# Make a tweedie distribution model of abundance
mod_tw <- dsm(count~s(X, Y, k = 50), # count is dependent on a spline of X and Y
              ddf.obj=hn1, # detecton function
              segment.data=segdata, # site data
              observation.data=hn1$ddf$data, # observation data (same as in the ds function)
              family=tw(),
              transect="point")

summary(mod_tw)
qq.gam(mod_tw, pch=16, cex=0.5, main="QQ-Plot",rep=100)
rqgam_check(mod_tw)

# Get standardised (1km2 predictions for each site)
preddata <- segdata %>% mutate(off.set=1e6)
pp <- preddata %>% mutate(dens = predict(mod_tw, preddata) %>% round(2) %>% as.numeric())

# Get the mean and sd of densities
mean(pp$dens)
sd(pp$dens)

# visualise density
mapview(x = pp,
        xcol = "X",
        ycol = "Y",
        zcol = "dens",
        grid = FALSE,
        crs = 3111)

# Now predict to an unsampled area
# Raster is of the proportion of each 1km2 grid covered by public land
offset.rast <- rast("example_data/camera_trap/offset.tif")

# Convert raster to data.frame and make off.set in m2
statewide_pred_data <- as.data.frame(offset.rast, xy = T) %>%
  transmute(X = x, Y = y, off.set = off.set*1e6)

statewide_pred <- statewide_pred_data %>%
  mutate(dens = predict(mod_tw, statewide_pred_data) %>% round(2) %>% as.numeric())

# total sum for state
sum(statewide_pred$dens)

# plot the predictions
prediction_raster <- offset.rast
values(prediction_raster)[!is.na(values(prediction_raster))] <- statewide_pred$dens
names(prediction_raster) <- "density_km2"

plot(prediction_raster)

# variance estimation
# Variance estimate assuming independence between components
sambar.var <- dsm_var_gam(mod_tw, statewide_pred_data, off.set=statewide_pred_data$off.set)
sambar.var

# Bootstrap variance
# calculate the variance by 100 moving block bootstraps
sambar.boot <- dsm_var_movblk(dsm.object = mod_tw,
                              pred.data = statewide_pred_data,
                              n.boot = 100,
                              block.size = 1,
                              samp.unit.name = "Sample.Label",
                              off.set = statewide_pred_data$off.set,
                              bar = TRUE,
                              ds.uncertainty = FALSE)
sambar.boot
