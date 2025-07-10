# function of alpha attenuation factor, suitable for coarse grains between 20 um and 1000 um, data from Brenann (1991) 

calc_alpha_attenuation <- function() {
  grain <- (grain.min + grain.max) / 2

  f_attenuation_alphaU <- function(x) {
    y <- 0.8624 * exp(-x / 23.477) + 0.2193 * exp(-x / 160.12) + 0.01983 # double exponential function fitted by origin
  }

  attenuation.alphaU <- f_attenuation_alphaU(grain)
  attenuation.alphaU.err <- 0.5 * (f_attenuation_alphaU(grain.min) - f_attenuation_alphaU(grain.max))

  f_attenuation_alphaTh <- function(x) {
    y <- 0.768 * exp(-x / 28.532) + 0.2741 * exp(-x / 122.47) + 0.03927 # double exponential function fitted by origin
  }

  attenuation.alphaTh <- f_attenuation_alphaTh(grain)
  attenuation.alphaTh.err <- 0.5 * (f_attenuation_alphaTh(grain.min) - f_attenuation_alphaTh(grain.max))

  c(attenuation.alphaU, attenuation.alphaU.err, attenuation.alphaTh, attenuation.alphaTh.err) # return 4 parameters
}
