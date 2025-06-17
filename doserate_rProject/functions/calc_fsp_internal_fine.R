# function for internal dose rate of 4-11 um

calc_fsp_internal_fine <- function() {
  
  # internal K and Rb inside K-feldspar, can be modified
  K.internal <- 12.5        ##Internal K of KF is set as 12.5+-0.5%, from Huntely and Baril (1997) and Zhao and Li (2005).
  K.internal.err <- 0.5     ## unit percent
  Rb.internal <- 400        ##Internal Rb of KF is 400+-100 ppm, from Huntley and Hancock (2001).
  Rb.internal.err <- 100    ##unit ppm
  
  #conversion parameters for K
  conversion.beta.K <- conversion[1,'K_beta']
  conversion.beta.K.err <- conversion[2,'K_beta']
  
  #beta absorption factor specifically for 4-11 um, from Guerin et al (2012)
  absorb.K <- 0.003 
  absorb.K.err <- 0.001
  
  #absorb.Rb of 4-11 um is from Readhead(2002). It include the conversion factor. Unit is uGy/a/(ppm Rb)
  absorb.Rb <- 0.022 
  absorb.Rb.err <- 0.007

  # to calculate dose rate
  internal <- absorb.K*conversion.beta.K*K.internal + absorb.Rb*Rb.internal/1000 
  internal.error <- sqrt((absorb.K*conversion.beta.K*K.internal)^2*((absorb.K.err/absorb.K)^2 + (conversion.beta.K.err/conversion.beta.K)^2 + (K.internal.err/K.internal)^2) + (absorb.Rb*Rb.internal/1000)^2*((absorb.Rb.err/absorb.Rb)^2 + (Rb.internal.err/Rb.internal)^2))
  #cat('Internal dose rate=',internal,'+-',internal.error,'Gy','\n\n')

  c(internal, internal.error) # The function will return these two values

}


