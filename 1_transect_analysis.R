# Analysing Transect Data ---------------------------------------------------------------------

#### Step 1: Setup ####
# install.packages("tidyverse")
library(tidyverse)
# install.packages("unmarked")
library(unmarked)
# install.packages("mapview")
library(mapview)


#### Step 2: Read in data ####
sambar_detection <- read_csv("example_data/transects/sambar_detection.csv")
site_data <- read_csv("example_data/transects/site_data.csv")

#### Step 3: Format detection data ####
# convert the data into observation covariates and presence-absence

# order the SiteIDs
sambar_detection_ordered <- sambar_detection %>%
  arrange(SiteID) %>%
  group_by(SiteID) %>% # add a number to each row for a site (used later)
  mutate(survey_n = 1:n()) %>%
  ungroup()

# get a presence-absence for the site
pa <- sambar_detection_ordered %>%
  group_by(SiteID) %>%
  summarise(pa = max(Presence)) %>%
  pull(pa)

# get site data that is in sambar_detection_ordered
site_covariates <- site_data %>%
  filter(SiteID %in% sambar_detection_ordered$SiteID)

# check same length
length(pa) == nrow(siteCovs)

# Make observation covariates for different survey methods.
det_covariates <- list()
det_covariates$Method <- pivot_wider(sambar_detection_ordered,
                                   id_cols = "SiteID",
                                   values_from = "Survey",
                                   names_from = "survey_n") %>%
  select(-SiteID)

# Make detection histories
y_hist <- pivot_wider(sambar_detection_ordered,
                      id_cols = "SiteID",
                      values_from = "Presence",
                      names_from = "survey_n") %>%
  select(-SiteID)

# Combine data into the one object
umf <- unmarkedFrameOccu(y = y_hist,
                     siteCovs = site_covariates,
                     obsCovs = det_covariates)

#### Step 4: Run analysis ####

# Run the model with the occu function
# BIO04 = temperature seasonality
# NPP = Net primary productivity
# SLOPE = slope of terrain
model_fit <- occu(formula = ~ Method ~ scale(BIO04) + scale(NPP) + scale(SLOPE),
                  data = umf,
                  linkPsi = "logit")

# get model summary
summary(model_fit)

# Plot the detection model effects
plotEffects(model_fit, "det", covariate = "Method")

# Plot the effect of BIO01, NPP and ForestEdge
par(mfrow = c(1, 3))
plotEffects(model_fit, "state", covariate = "BIO04")
plotEffects(model_fit, "state", covariate = "NPP")
plotEffects(model_fit, "state", covariate = "SLOPE")

#### Step 5: Generate predictions ####
# you can back-predict to the sites to get average occupancy
# state = occupancy, det = detection
site_predictions <- predict(model_fit, type = "state")

# combine the predictions with our site data and make a little interactive map
map_data <- site_covariates %>%
  select(SiteID, BIO04, NPP, SLOPE, Latitude, Longitude) %>%
  bind_cols(site_predictions)

str(map_data)

mapview(x = map_data,
        xcol = "Longitude",
        ycol = "Latitude",
        zcol = "Predicted",
        grid = FALSE,
        crs = 4283)
