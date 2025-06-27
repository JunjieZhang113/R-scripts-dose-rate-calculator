# R script for environmental dose rate calculation in luminescence and ESR dating

# Open the doserate_rProject.Rproj file inside 'doserate_rProject' folder.
# Then, in the 'Files' panel (right-bottom), double click the doserate_main.R file to open it.
# In this way, the working directory will be automatically set to the location of the 'doserate_rProject' folder.
# If not, the working directory can also be set manually, by setwd(".../doserate_rProject"). '...'is the location of the 'doserate_rProject' folder.


rm(list = ls()) # to clear the workspace
getwd()         # to check whether the current working directory. It should be ".../doserate_rProject". If it is not, use setwd()

Dr <- read.csv("Template_input.csv", header = TRUE) # load samples

print(Dr) # to check whether the data has been successfully loaded or not
          # Be careful, if a space ' ' is entered into empty rows inside the csv template by an accident, 
          # R will also take that row as a sample (with 'NA' values) and will give error message when the calculation goes to that 'sample'.
         
output_filename <- "Doserate_output.csv"  # give the name of the output csv file

conversion <- read.csv("conversion_data/Liritzis_2013_conversion.csv", header = TRUE) # load the conversion factors, 3 options: Guerin_2011_conversion, Liritzis_2013_conversion, Cresswell_2018_conversion

print(conversion) # to check whether it is successfully loaded or not


#################################################################################
### The code below are for the calculations. Select all (ctrl+A), and run.#######
#################################################################################

# load the functions
source("functions/calc_cosmicray.R")
source("functions/calc_fsp_internal_coarse.R")
source("functions/calc_fsp_internal_fine.R")
source("functions/calc_dr_alpha_coarse.R")
source("functions/calc_dr_alpha_fine.R")
source("functions/calc_dr_beta_coarse.R")
source("functions/calc_dr_beta_fine.R")
source("functions/calc_dr_gamma.R")

# to calculate the samples one by one
for (i in 1:nrow(Dr)) {
  # mineral and grain size classification
  mineral <- Dr[i, "mineral"]
  grainsize <- Dr[i, "grainsize"]

  # etch depth for coarse grains
  etch.depth <- Dr[i, "etch_depth_um"] * 2 ## to get the reduction in grain diameter

  # alpha efficiency
  alpha.value <- Dr[i, "a_value"]
  alpha.value.err <- Dr[i, "a_value_err"]

  # radon (222Rn) loss. Values from 0 to 1. e.g., 0 means no Rn loss, 0.2 means 20% of the radon is escaped.
  Radon.loss <- Dr[i, "Rn_loss"]

  # geography information
  geo.lat <- Dr[i, "latitude"] # latitude
  geo.lon <- Dr[i, "longitude"] # longitude
  altitude <- Dr[i, "altitude_m"] # altitude，unit m
  depth <- Dr[i, "depth_m"] # the sampling depth from the surface，unit m

  # grain size
  grain.min <- Dr[i, "grain_min"] # the minimum grain size，unit um
  grain.max <- Dr[i, "grain_max"] # the maximum grain size，unit um

  # water
  water <- Dr[i, "water_percent"] / 100
  water.error <- Dr[i, "water_err"] / 100

  # U, Th, K concentrations
  U <- Dr[i, "U_ppm"]
  U.error <- Dr[i, "U_err"]
  Th <- Dr[i, "Th_ppm"]
  Th.error <- Dr[i, "Th_err"]
  K <- Dr[i, "K_percent"]
  K.error <- Dr[i, "K_err"]

  if (!(mineral %in% c("fsp", "qz"))) {
    stop("The input of 'mineral' should be either 'fsp' or 'qz' in the template file.")
  }
  if (!(grainsize %in% c("coarse", "fine"))) {
    stop("The input of 'grainsize' should be either 'coarse' or 'fine' in the template file.")
  }


  ############## cosmic-ray dose rate ###########
  cosmic <- calc_cosmicray()
  Dr[i, "cosmicray"] <- cosmic[1]
  Dr[i, "cosmicray_err"] <- cosmic[2]

  ############# internal dose rate ##############
  if (mineral == "qz") {
    internal <- c(0, 0) # for quartz, internal dose rate is set as 0.
    Dr[i, "internal"] <- internal[1]
    Dr[i, "internal_err"] <- internal[2]
  }

  if (mineral == "fsp") {
    if (grainsize == "coarse") {
      internal <- calc_fsp_internal_coarse()
      Dr[i, "internal"] <- internal[1]
      Dr[i, "internal_err"] <- internal[2]
    }

    if (grainsize == "fine") {
      internal <- calc_fsp_internal_fine()
      Dr[i, "internal"] <- internal[1]
      Dr[i, "internal_err"] <- internal[2]
    }
  }

  ################ external dose rate  ###################
  # external alpha dose rate
  if (grainsize == "coarse") {
    #  if no etching, calculate the alpha dose rate
    if (etch.depth == 0) {
      alpha <- calc_dr_alpha_coarse()
    }
    # if etched, just assume alpha dose rate to be 0
    if (etch.depth != 0) {
      alpha <- c(0, 0)
    }

    Dr[i, "external_alpha"] <- alpha[1]
    Dr[i, "external_alpha_err"] <- alpha[2]
  }

  if (grainsize == "fine") {
    alpha <- calc_dr_alpha_fine()
    Dr[i, "external_alpha"] <- alpha[1]
    Dr[i, "external_alpha_err"] <- alpha[2]
  }

  # external beta dose rate
  if (grainsize == "coarse") {
    beta <- calc_dr_beta_coarse()
    Dr[i, "external_beta"] <- beta[1]
    Dr[i, "external_beta_err"] <- beta[2]
  }

  if (grainsize == "fine") {
    beta <- calc_dr_beta_fine()
    Dr[i, "external_beta"] <- beta[1]
    Dr[i, "external_beta_err"] <- beta[2]
  }

  # external gamma dose rate
  gamma <- calc_dr_gamma()
  Dr[i, "external_gamma"] <- gamma[1]
  Dr[i, "external_gamma_err"] <- gamma[2]

  # total external dose rate
  Dr.external <- alpha[1] + beta[1] + gamma[1]
  Dr.external.error <- sqrt(alpha[2]^2 + beta[2]^2 + gamma[2]^2)
  Dr[i, "external_total"] <- Dr.external
  Dr[i, "external_total_err"] <- Dr.external.error

  # total  dose rate
  Dr.total <- Dr.external + internal[1] + cosmic[1]
  Dr.total.error <- sqrt(Dr.external.error^2 + internal[2]^2 + cosmic[2]^2)
  Dr[i, "doserate_Gyperka"] <- Dr.total
  Dr[i, "doserate_err"] <- Dr.total.error

  if (i == nrow(Dr)) {
    write.csv(Dr, output_filename, row.names = FALSE)
    cat(
      "\n####################### Well done! ##########################\n",
      "Dose rate calculation of", nrow(Dr), "samples is finished.\n",
      "The results are saved in the csv file:", output_filename, "\n",
      "The dose rates of the first five samples are shown below:\n\n "
    )
    print(Dr[1:5, c("Sample", "mineral", "grain_min", "grain_max", "doserate_Gyperka", "doserate_err")])
  }
}
