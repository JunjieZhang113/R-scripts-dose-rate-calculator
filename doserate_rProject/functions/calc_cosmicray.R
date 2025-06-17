# function for cosmic ray dose rate calculation

calc_cosmicray <- function() {
  
  density <- 2.0 #Density is set as 2.0 g cm-3
  
  geomag.lat <- 180/pi*asin(0.203*cos(geo.lat*pi/180)*cos((geo.lon-291)*pi/180)+0.979*sin(geo.lat*pi/180))
  geomag.lat.abs<-abs(geomag.lat)
  
  #To set the three parameters Fi，J，H in the function D=D0*[Fi+J*exp(-h/H)]
  #From Prescott and Hutton，1994 (Fig.2).Do not change these values
  if(geomag.lat.abs<2.5)                         {Fi <- 0.40; J <- 0.52; H <- 4.40}
  if(2.5<=geomag.lat.abs & geomag.lat.abs<7.5)   {Fi <- 0.39; J <- 0.54; H <- 4.38}
  if(7.5<=geomag.lat.abs & geomag.lat.abs<12.5)  {Fi <- 0.38; J <- 0.55; H <- 4.36}
  if(12.5<=geomag.lat.abs & geomag.lat.abs<17.5) {Fi <- 0.37; J <- 0.57; H <- 4.33}
  if(17.5<=geomag.lat.abs & geomag.lat.abs<22.5) {Fi <- 0.35; J <- 0.60; H <- 4.30}
  if(22.5<=geomag.lat.abs & geomag.lat.abs<27.5) {Fi <- 0.32; J <- 0.63; H <- 4.25}
  if(27.5<=geomag.lat.abs & geomag.lat.abs<32.5) {Fi <- 0.28; J <- 0.69; H <- 4.18}
  if(32.5<=geomag.lat.abs & geomag.lat.abs<37.5) {Fi <- 0.25; J <- 0.75; H <- 4.12}
  if(37.5<=geomag.lat.abs)                       {Fi <- 0.24; J <- 0.76; H <- 4.10}
  
  #calculate the cosmic ray
  x <- density*depth   # unit 100 g cm-2
  if(x<1.7) {D0 <- 0.115*exp(-x*100/77.86)+0.172}      #fitted from data of Prescott and Hutton (1988), containing soft and hard components
  else {D0 <- exp(-0.00055*x)*6072/(((x+11.6)^1.68+75)*(x+212))} #Prescott and Hutton (1994), only  hard component is significant
  
  cosmic.ray <- D0*(Fi+J*exp(altitude/1000/H)) 
  cosmic.error <- 0.1*cosmic.ray #usually it is taken as 10 %
  
  c(cosmic.ray, cosmic.error) #return cosmic dose rate and error
  
}

