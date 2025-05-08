## Dose rate calculation for 4-11 um K-Feldspar (polymineral fine grains)

## Input data into the csv file - Template 4-11um Feldspar.csv, before running the codes.

## The units of U, Th is ppm, the unit of K (not K2O) is %.
## For the values of the latitude and longitude, it is positive values when it is north and east, otherwise it is negative.

##------------parameters used for calculation---------------###
##Conversion factors are from Liritzis et al. (2013); 
##Alpha-dose attenuation factor is fitted from Brenann (1991);
##Beta absorption factor is from Guerin et al. (2012);
##Absorbed dose fraction of Rb is from Readhead (2002);
##Cosmic ray calculation from Prescott and Hutton (1994);

rm(list = ls())
setwd("D:/R scripts for dose rate calculation/R scripts and templates") 
## change the working directory to the folder location on your computer

####### Prepare the CSV template, select all the script (ctrl+A), and run  ################
###########################################################################################

Dr<-read.csv("Template 4-11um Feldspar.csv",header = TRUE)

result_filename<-"Dr results feldspar 4-11 um.csv" 
#result_filename can be changed as you wish

########################  parameters ##############
Radon.loss<-0      #For radon (222Rn) loss. Values from 0 to 1. For example, 0 means no Rn loss, and 0.2 means 20% of the radon is escaped.The degree of Rn loss can be deduced from 210Pb/(214Pb & Bi214) (De Corte, 2006) if gamma spectrometry is used for U, Th, K measurements.

########################  alpha efficiency value ##############
alpha.factor<-0.09     #It is the a-value (i.e., k3.7 value) from Schmidt et al. 2018. It will be converted to effective k-value in the code below. keff = 0.80*K3.7 for U, and 0.86*k3.7 for Th (Zimmerman 1971, page 313 in Aitken 1985)
alpha.factor.err<-0.02 #Error of alpha factor, here it is set as 0.02.
########################

## internal K and Rb inside K-feldspar####
K.internal<-12.5        ##percent. For K-feldspar, internal K is set as 12.5+-0.5%. Huntely and Baril (1997). Zhao and Li (2005).
K.internal.err<-0.5     ##percent
Rb.internal<-400        ##ppm, internal Rb of K-feldspar is 400+-100 ppm. Huntley and Hancock (2001).
Rb.internal.err<-100    ##ppm
############################################################

####################### don't change the following codes ###################

###U unit 1 ppm
conversion.alpha.U.preRadon<-1.265
conversion.alpha.U.preRadon.err<-0.011
conversion.alpha.U.afterRadon<-2.793-1.265
conversion.alpha.U.afterRadon.err<-0.011

conversion.beta.U.preRadon<-0.0602
conversion.beta.U.preRadon.err<-0.0002
conversion.beta.U.afterRadon<-0.1459-0.0602
conversion.beta.U.afterRadon.err<-0.0004

conversion.gamma.U.preRadon<-0.0044
conversion.gamma.U.preRadon.err<-0.0000
conversion.gamma.U.afterRadon<-0.1118-0.0044
conversion.gamma.U.afterRadon.err<-0.0002

###Th unit 1ppm
conversion.alpha.Th<-0.7375
conversion.alpha.Th.err<-0.0026

conversion.beta.Th<-0.0275
conversion.beta.Th.err<-0.0009

conversion.gamma.Th<-0.0481
conversion.gamma.Th.err<-0.0002

### K unit 1%
conversion.beta.K<-0.8011
conversion.beta.K.err<-0.0073

conversion.gamma.K<-0.2498
conversion.gamma.K.err<-0.0048


##Rb  for 50 ppm Rb (may be used)
conversion.beta.Rb<-0.0185
conversion.beta.Rb.err<-0.0004

for(i in 1:nrow(Dr))  ##recycling from the first row to the last row.
    {
  geo.lat<-Dr[i,2]        #latitude
  geo.lon<-Dr[i,3]        #longitude
  altitude<-Dr[i,4]        #altitude，unit m
  depth<-Dr[i,5]           #the sampling depth from the surface，unit m
  grain.min<-Dr[i,6]       #the minimum grain size，unit um
  grain.max<-Dr[i,7]       #the maximum grain size，unit um
  
  water=Dr[i,8]/100     
  U=Dr[i,10]
  Th=Dr[i,12]
  K=Dr[i,14] 
  
  water.error=Dr[i,9]/100
  U.error=Dr[i,11]
  U.error=U*sqrt((U.error/U)^2+0.05^2)  #add 5 % system error
  Th.error=Dr[i,13]
  Th.error=Th*sqrt((Th.error/Th)^2+0.05^2)
  K.error=Dr[i,15]
  K.error=K*sqrt((K.error/K)^2+0.05^2)
  ######################################### cosmic ray dose rate calculation #############################
  
  geomag.lat=180/pi*asin(0.203*cos(geo.lat*pi/180)*cos((geo.lon-291)*pi/180)+0.979*sin(geo.lat*pi/180))
  geomag.lat.abs<-abs(geomag.lat)
  #To set the three parameters Fi,J,H in the function D=D0*[Fi+J*exp(-h/H)]
  #from Prescott and Hutton，1994 (Fig.2).Do not change these values.
  if(geomag.lat.abs<2.5)                     {Fi=0.40; J=0.52; H=4.40}
  if(2.5<=geomag.lat.abs & geomag.lat.abs<7.5)   {Fi=0.39; J=0.54; H=4.38}
  if(7.5<=geomag.lat.abs & geomag.lat.abs<12.5)  {Fi=0.38; J=0.55; H=4.36}
  if(12.5<=geomag.lat.abs & geomag.lat.abs<17.5) {Fi=0.37; J=0.57; H=4.33}
  if(17.5<=geomag.lat.abs & geomag.lat.abs<22.5) {Fi=0.35; J=0.60; H=4.30}
  if(22.5<=geomag.lat.abs & geomag.lat.abs<27.5) {Fi=0.32; J=0.63; H=4.25}
  if(27.5<=geomag.lat.abs & geomag.lat.abs<32.5) {Fi=0.28; J=0.69; H=4.18}
  if(32.5<=geomag.lat.abs & geomag.lat.abs<37.5) {Fi=0.25; J=0.75; H=4.12}
  if(37.5<=geomag.lat.abs)                       {Fi=0.24; J=0.76; H=4.10}
  
  #calculate the cosmic ray
  x=2.0*depth   ##Density is set as 2.0 g/cm3, unit 100 g/cm-2
  if(x<1.7) {D0=0.115*exp(-x*100/77.86)+0.172}      #fitted from data of Prescott and Hutton (1988), containing soft and hard components
  else {D0=exp(-0.00055*x)*6072/(((x+11.6)^1.68+75)*(x+212))} #Prescott and Hutton (1994), only  hard component is significant
  
  cosmic.ray=D0*(Fi+J*exp(altitude/1000/H)) 
  Dr[i,16]=cosmic.ray
  cosmic.error=0.1*Dr[i,16] #usually it is taken as 10 %
  Dr[i,17]=cosmic.error
  
  
  #### alpha attenuation factor for 4-11 um ################################# 
  # the alpha-dose attenuation factor is fitted from Brenann (1991). 
  
  attenuation.alphaU<-0.84
  attenuation.alphaU.err<-0.075
  attenuation.alphaTh<-0.86
  attenuation.alphaTh.err<-0.06

  ######### beta absorption factor for 4-11 um #################################   
  #the beta absorption factor is from Guerin (2012), cite from DARC
  absorb.K=0.003
  absorb.K.err=0.001
  
  absorb.U=0.015
  absorb.U.err=0.007
  
  absorb.Th=0.019
  absorb.Th.err=0.009
  
  #absorb.Rb is from Readhead(2002), unit �Gy/a/(ppm Rb), for internal dose
  absorb.Rb=0.022
  absorb.Rb.err=0.007
  
  
  sampleID = as.character(Dr[i,1])
  cat('*******************************Sample ID: ',sampleID, '*******************************************************************','\n\n')
  
   ############################################################ internal dose rate ###########################################
  Dr.internal=absorb.K*conversion.beta.K*K.internal+absorb.Rb*Rb.internal/1000 
  Dr.internal.error=sqrt((absorb.K*conversion.beta.K*K.internal)^2*((absorb.K.err/absorb.K)^2+(conversion.beta.K.err/conversion.beta.K)^2+(K.internal.err/K.internal)^2)+(absorb.Rb*Rb.internal/1000)^2*((absorb.Rb.err/absorb.Rb)^2+(Rb.internal.err/Rb.internal)^2))
  cat('Internal dose rate=',Dr.internal,'+-',Dr.internal.error,'Gy','\n\n')
  
  Dr[i,18]=Dr.internal
  Dr[i,19]=Dr.internal.error
  
  ########################################################## external dose rate ############################################

    ########### alpha dose rate ###########
    
    # Note that 0.80 for U and 0.86 for Th are the conversion coefficients from k3.7 (a-value) to keff.
  
    alpha.U=0.80*U*attenuation.alphaU*alpha.factor*(conversion.alpha.U.preRadon+(1-Radon.loss)*conversion.alpha.U.afterRadon)
    alpha.U.error=alpha.U*sqrt((U.error/U)^2+(attenuation.alphaU.err/attenuation.alphaU)^2+(alpha.factor.err/alpha.factor)^2+(conversion.alpha.U.preRadon.err^2+(1-Radon.loss)^2*conversion.alpha.U.afterRadon.err^2)/(conversion.alpha.U.preRadon+(1-Radon.loss)*conversion.alpha.U.afterRadon)^2)
    cat('Individual alpha dose rate before water:','\n',' alpha.U=',alpha.U,'+-',alpha.U.error,'Gy','\n')
    
    alpha.Th=0.86*Th*attenuation.alphaTh*alpha.factor*conversion.alpha.Th
    alpha.Th.error=alpha.Th*sqrt((Th.error/Th)^2+(attenuation.alphaTh.err/attenuation.alphaTh)^2+(alpha.factor.err/alpha.factor)^2+(conversion.alpha.Th.err/conversion.alpha.Th)^2)  
    cat('  alpha.Th=',alpha.Th,'+-',alpha.Th.error,'Gy','\n')
    
    alpha=(alpha.U+alpha.Th)/(1+1.50*water)
    alpha.error=sqrt((1/(1+1.5*water))^2*alpha.U.error^2+(1/(1+1.5*water))^2*alpha.Th.error^2+((alpha.U+alpha.Th)/(1+1.5*water)^2)^2*1.5^2*water.error^2)
    cat('Alpha dose rate after water=',alpha,'+-',alpha.error,'Gy','\n\n')
    
    Dr[i,20]=alpha
    Dr[i,21]=alpha.error
    
    ########## beta dose rate ##############
    
    beta.K=K*(1-absorb.K)*conversion.beta.K;
    beta.K.error=beta.K*sqrt((K.error/K)^2+absorb.K.err^2/(1-absorb.K)^2+(conversion.beta.K.err/conversion.beta.K)^2)
    cat('Individual beta dose rate before water:','\n',' beta.K=',beta.K,'+-',beta.K.error,'Gy','\n')
    
    beta.U=U*(1-absorb.U)*(conversion.beta.U.preRadon+(1-Radon.loss)*conversion.beta.U.afterRadon);
    beta.U.error=beta.U*sqrt((U.error/U)^2+absorb.U.err^2/(1-absorb.U)^2+(conversion.beta.U.preRadon.err^2+(1-Radon.loss)^2*conversion.beta.U.afterRadon.err^2)/(conversion.beta.U.preRadon+(1-Radon.loss)*conversion.beta.U.afterRadon)^2)
    cat('  beta.U=',beta.U,'+-',beta.U.error,'Gy','\n')
    
    beta.Th=Th*(1-absorb.Th)*conversion.beta.Th;
    beta.Th.error=beta.Th*sqrt((Th.error/Th)^2+absorb.Th.err^2/(1-absorb.Th)^2+(conversion.beta.Th.err/conversion.beta.Th)^2)
    cat('  beta.Th=',beta.Th,'+-',beta.Th.error,'Gy','\n')
    
    beta=(beta.K+beta.U+beta.Th)/(1+1.25*water);
    beta.error=sqrt((1/(1+1.25*water))^2*beta.K.error^2+(1/(1+1.25*water))^2*beta.U.error^2+(1/(1+1.25*water))^2*beta.Th.error^2+((beta.K+beta.U+beta.Th)/(1+1.25*water)^2)^2*1.25^2*water.error^2)

    cat('Beta dose rate after water=',beta,'+-',beta.error,'Gy','\n\n')
    
    Dr[i,22]=beta
    Dr[i,23]=beta.error

    
    ########## gamma dose rate ##############    
    
    gamma.K=K*conversion.gamma.K
    gamma.K.error=gamma.K*sqrt((K.error/K)^2+(conversion.gamma.K.err/conversion.gamma.K)^2)
    cat('Individual gamma dose rate before water:','\n',
        ' gamma.K=',gamma.K,'+-',gamma.K.error,'Gy','\n')
    
    gamma.U=U*(conversion.gamma.U.preRadon+(1-Radon.loss)*conversion.gamma.U.afterRadon)
    gamma.U.error=gamma.U*sqrt((U.error/U)^2+(conversion.gamma.U.preRadon.err^2+(1-Radon.loss)^2*conversion.gamma.U.afterRadon.err^2)/(conversion.gamma.U.preRadon+(1-Radon.loss)*conversion.gamma.U.afterRadon)^2)
    cat('  gamma.U=',gamma.U,'+-',gamma.U.error,'Gy','\n')
    
    gamma.Th=Th*conversion.gamma.Th
    gamma.Th.error=gamma.Th*sqrt((Th.error/Th)^2+(conversion.gamma.Th.err/conversion.gamma.Th)^2)
    cat('  gamma.Th=',gamma.Th,'+-',gamma.Th.error,'Gy','\n')
    
    gamma=(gamma.K+gamma.U+gamma.Th)/(1+1.14*water);
    gamma.error=sqrt((1/(1+1.14*water))^2*gamma.K.error^2+(1/(1+1.14*water))^2*gamma.U.error^2+(1/(1+1.14*water))^2*gamma.Th.error^2+((gamma.K+gamma.U+gamma.Th)/(1+1.14*water)^2)^2*1.14^2*water.error^2)
    
    Dr[i,24]=gamma
    Dr[i,25]=gamma.error
    cat('Gamma dose rate after water=',gamma,'+-',gamma.error,'Gy','\n\n')
    
    ######################################################################
    
    Dr.external=alpha+beta+gamma
    Dr.external.error=sqrt(alpha.error^2+beta.error^2+gamma.error^2)
    Dr[i,26]=Dr.external
    Dr[i,27]=Dr.external.error
    
    Dr.total=Dr.internal+Dr.external+cosmic.ray;
    Dr.total.error=sqrt(Dr.internal.error^2+Dr.external.error^2+cosmic.error^2)
    Dr[i,28]<-Dr.total
    Dr[i,29]<-Dr.total.error
   
    ######################################
    
    if(i==nrow(Dr)) 
      {
      filename<-result_filename
      write.csv(Dr,filename,row.names=FALSE)
      cat('\n','*********************Well done!**************************************************',"\n",'The results are written in the file:',filename,"\n",
          "The dose rates and errors are also displayed below.","\n\n")
     print(Dr[,c(1,16,17,18,19,20,21,22,23,24,25,26,27,28,29)])
      }
}





