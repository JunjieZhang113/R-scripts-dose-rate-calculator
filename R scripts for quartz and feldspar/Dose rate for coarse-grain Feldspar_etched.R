##Dose rate calculation for K-Feldspar coarse grains
##surface layer etched by HF. alpha dose rate is not considered.

##Applicable for grain size in the range between 20 and 1000 um,
##because alpha attenuation factor and beta absorption factor are fitted within this range

##Input data into the csv file - Template coarse-grain Feldspar.csv

##The units of U,Th is ppm, the unit of K (not K2O) is %.
## For the values of the latitude and longitude, it is positive values when it is north and east, otherwise it is negative.

##------------parameters used for calculation---------------###
##Conversion factors are from Liritzis et al. (2013); 
##Beta absorption factor is from Guerin et al. (2012);
##Absorbed dose fraction of Rb is from Readhead (2002);
##Cosmic ray calculation from Prescott and Hutton (1994);

rm(list = ls())
setwd("D:/R scripts for dose rate calculation/R scripts and templates")  ## change the working directory to the folder location on your computer


####### Prepare the CSV template, select all the script (ctrl+A), and run  ################

Dr<-read.csv("Template coarse-grain Feldspar.csv",header = TRUE)
#head(Dr) ###check the 'read in' is successful

result_filename<-"Dr results feldspar coarse grain - etched 20perRnloss.csv" 
#result_filename can be changed as you wish

########################  parameters #######################
Radon.loss<-0.20      #For radon (222Rn) loss. Values from 0 to 1. For example, 0 means no Rn loss, and 0.2 means 20% of the radon is escaped.The degree of Rn loss can be deduced from 210Pb/(214Pb & Bi214) (De Corte, 2006), if gamma spectrometry is used for U, Th, K measurements.

## internal K and Rb inside K-feldspar####
K.internal<-12.5        ##percent. For K-feldspar, internal K is set as 12.5+-0.5%. Huntely and Baril (1997). Zhao and Li (2005).
K.internal.err<-0.5     ##percent
Rb.internal<-400        ##ppm, internal Rb of K-feldspar is 400+-100 ppm. Huntley and Hancock (2001).
Rb.internal.err<-100    ##ppm
############################################################

etch.depth<-10      #the thickness of the layer of K-feldspar removed by HF etching (unit um) 
########################

###conversion factors are from Liritzis 2013 (Table 1-3).
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
  Th.error=Dr[i,13]
  K.error=Dr[i,15]
  
  ### add 5 % error in quadrature for U, Th, K measured by gamma spectrometry in LIAG. Not applicable for others.
  #U.error=U*sqrt((U.error/U)^2+0.05^2) 
  #Th.error=Th*sqrt((Th.error/Th)^2+0.05^2)
  #K.error=K*sqrt((K.error/K)^2+0.05^2)
  
  ######################################### cosmic ray dose rate calculation #############################
  
  geomag.lat=180/pi*asin(0.203*cos(geo.lat*pi/180)*cos((geo.lon-291)*pi/180)+0.979*sin(geo.lat*pi/180))
  geomag.lat.abs<-abs(geomag.lat)
  #To set the three parameters Fi，J，H in the function D=D0*[Fi+J*exp(-h/H)]
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
  
  

  ########################################################## beta absorption factor #################################   
  #calculate the beta absorption factor, function fitted from data in Guerin (2012), suitable for 20-1000um
  grain<-(grain.min+grain.max)/2
  absorb.K=-6.86E-4+3.66E-4*grain
  absorb.K.err=(3.66E-4*grain.max-3.66E-4*grain.min)/2
  
  f_absorb_U <- function(x) {
    y <- 0.01505+0.6178*(1-exp(-x/1273.14))+0.03534*(1-exp(-x/44.40)) #fitted in Origin
    return(y)
  }
  
  absorb.U<-f_absorb_U(grain)
  absorb.U.err<-0.5*(f_absorb_U(grain.max)-f_absorb_U(grain.min))
  
  f_absorb_Th <- function(x) {
    y <- 0.02003+0.5690*(1-exp(-x/1128.41))+0.10357*(1-exp(-x/95.99)) #fitted in Origin
    return(y)
  }
  
  absorb.Th<-f_absorb_Th(grain)
  absorb.Th.err<-0.5*(f_absorb_Th(grain.max)-f_absorb_Th(grain.min))
  
  #absorb.Rb is fitted from Readhead(2002), unit µGy/a/(ppm Rb), for internal dose, it include the conversion factor
  absorb.Rb=0.00735+0.30394*(1-exp(-grain/125.6)) 
  absorb.Rb.err=(0.30394*(1-exp(-grain.max/125.6))-0.30394*(1-exp(-grain.min/125.6)))/2
 
  ##combined etching factor: the degree of increase of the beta absorption factor after etching, etch depth from 0-40 um, data from Brennan (2003)in DARC
  etch.factor<-1.00+0.00954*etch.depth-1.618E-4*etch.depth^2+1.425E-6*etch.depth^3
  
 
  #output the parameters
   sampleID = as.character(Dr[i,1])
   cat('*********************  Sample ID: ',sampleID, '*********************************************************','\n', 
       'Grain size:',grain.min,'-',grain.max,'um.  Mean grain size:',grain,'um.', '\n',
       "absorb.K=",absorb.K,'+-',absorb.K.err,"\n",
       "absorb.U=",absorb.U,'+-',absorb.U.err,"\n",
       "absorb.Th=",absorb.Th, '+-',absorb.Th.err,"\n",
       "absorb.Rb=",absorb.Rb, '+-',absorb.Rb.err,"\n",
       'etching depth=',etch.depth,'um',"\n",
       'etch.factor=',etch.factor,"\n",

        "\n")
   
  
  ############################################################ internal dose rate ###########################################
  Dr.internal=absorb.K*etch.factor*conversion.beta.K*K.internal+absorb.Rb*etch.factor*Rb.internal/1000 
  Dr.internal.error=sqrt((absorb.K*etch.factor*conversion.beta.K*K.internal)^2*((absorb.K.err/absorb.K)^2+(conversion.beta.K.err/conversion.beta.K)^2+(K.internal.err/K.internal)^2)+(absorb.Rb*etch.factor*Rb.internal/1000)^2*((absorb.Rb.err/absorb.Rb)^2+(Rb.internal.err/Rb.internal)^2))
  
  cat('Internal dose rate=',Dr.internal,'+-',Dr.internal.error,'Gy','\n\n')
  
  Dr[i,18]=Dr.internal
  Dr[i,19]=Dr.internal.error
  
  ########################################################## external dose rate ############################################

    ########### alpha dose rate ###########
    # since the grains are etched, we just set alpha dose rate to be 0.
  
    alpha=0
    alpha.error=0
  
    Dr[i,20]=alpha
    Dr[i,21]=alpha.error
    
    ########## beta dose rate ##############
    
    beta.K=K*(1-etch.factor*absorb.K)*conversion.beta.K;
    beta.K.error=beta.K*sqrt((K.error/K)^2+etch.factor^2*absorb.K.err^2/(1-etch.factor*absorb.K)^2+(conversion.beta.K.err/conversion.beta.K)^2)
    cat('Individual external beta dose rate before water:','\n',' beta.K=',beta.K,'+-',beta.K.error,'Gy','\n')
    
    beta.U=U*(1-etch.factor*absorb.U)*(conversion.beta.U.preRadon+(1-Radon.loss)*conversion.beta.U.afterRadon);
    beta.U.error=beta.U*sqrt((U.error/U)^2+etch.factor^2*absorb.U.err^2/(1-etch.factor*absorb.U)^2+(conversion.beta.U.preRadon.err^2+(1-Radon.loss)^2*conversion.beta.U.afterRadon.err^2)/(conversion.beta.U.preRadon+(1-Radon.loss)*conversion.beta.U.afterRadon)^2)
    cat('  beta.U=',beta.U,'+-',beta.U.error,'Gy','\n')
    
    beta.Th=Th*(1-etch.factor*absorb.Th)*conversion.beta.Th
    beta.Th.error=beta.Th*sqrt((Th.error/Th)^2+etch.factor^2*absorb.Th.err^2/(1-etch.factor*absorb.Th)^2+(conversion.beta.Th.err/conversion.beta.Th)^2)
    cat('  beta.Th=',beta.Th,'+-',beta.Th.error,'Gy','\n')
    
    beta=(beta.K+beta.U+beta.Th)/(1+1.25*water);
    beta.error=sqrt((1/(1+1.25*water))^2*beta.K.error^2+(1/(1+1.25*water))^2*beta.U.error^2+(1/(1+1.25*water))^2*beta.Th.error^2+((beta.K+beta.U+beta.Th)/(1+1.25*water)^2)^2*1.25^2*water.error^2)
    
    cat('External beta dose rate after water=',beta,'+-',beta.error,'Gy','\n\n')
    
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
    cat('Gamma dose rate after water=',gamma,'+-',gamma.error,'Gy','\n\n')
    
    Dr[i,24]=gamma
    Dr[i,25]=gamma.error
    
    
    ########################################### total dose rate ################################################################
    
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
        write.csv(Dr,result_filename,row.names=FALSE)
      cat('\n',"Well done!","\n",'The results are written in the file:',result_filename,"\n",
          "The dose rates and errors are also displayed below.","\n\n")
      print(Dr[,c(1,16,17,18,19,20,21,22,23,24,25,26,27,28,29)])
      }
}





