###This code was originally authored by Devin W. Goodsman, postdoctoral researcher at LANL#####

##############################################################################################
# This software has been authored by an employee or employees of Los Alamos National Security,
# LLC, operator of the Los Alamos National Laboratory (LANL) under contract No. 
# DE-AC52-06NA25396 with the U.S. Department of Energy. The U.S. government has rights to use,
# reproduce, and distribute this software. The public may copy, distribute, prepare derivative
# works and publicly display this software without charge, provided that this Notice and
# any statement of authorship are reproduced in all copies. Neither the government nor LANS
# makes any warranty, express or implied, or assumes any liability or responsibility for the
# use of this software. If software is modified to produce derivative works, such modified
# software should be clearly marked, so as not to confuse it with the version available from
# LANL
##############################################################################################

# This version of the IMAPfunction is designed to simulate mountain pine beetle phenology 
# given a vector of daily temperatures that is 443 days long. This version has been extensively 
# de-bugged and should be reliable.

# This version simulates phenology and egg, pupae, teneral adult, and 
# adult mortality only (no larval mortality).

# I call it IMAPfunctionNewJ. J is for Jasper.

# Here I take my phenology and mortality model and turn it into a function that can be applied in a loop
# This version differs from the original (Imapfunction.R) in that intead of taking as input mean air 
# temperatures, we compute adjusted mean phloem temperatures from the maximum and minimum air temperatures
# as described in Bolstad et al.

# Here I compute the analytic solution for the expectation of 
# the mountain pine beetle phenology.

######################################################################
# The rate equations are given in Regniere 2012
# I assume (as they did) that there is log-normal variability
# in development rates.
######################################################################

# The ImapFunc function takes as input
# a vector of minimum temperatures (Tmin), and a vector of maximum
# temperatures (Tmax).

# The ImapFuncNewJ function returns as output a matrix
# with rows that comprise
# vectors of the number of individuals in each of the 
# nine life stages (9 rows and 443 columns). 
# Note also that for this version, Tmin and Tmax
# must be exactly 443 days in length.

ImapFuncNewJ = function(Tmin, Tmax, StartT){
  #########################################################################################
  # Components of the phenology model
  #########################################################################################
  
  ## We need to compute the mean phloem temperature according to Bolstad, Bentz and Logan.
  ## We compute mean phloem temperature by averaging maximum and minimum phloem temperature.
  
  # Here I use the average temperature differential (6.6)
  Temp = 0.5*(Tmax + Tmin) + 0.9 + 6.6*(Tmax - Tmin)/(2*24.4)
  Tmax2 = Tmax + 6.6*(Tmax - Tmin)/24.4
  Tmin2 = Tmin + 1.8
  
  # The Regniere function for temperature dependent beetle development
  RegniereFunc = function(TC, TB, DeltaB, TM, DeltaM, omega, psi){
    
    Y = 0.0
    
    if(TC >= TB & TC <= TM){
      Y = psi*(exp(omega*(TC - TB)) - (TM - TC)/(TM - TB)*exp(-omega*(TC - TB)/DeltaB)
               - (TC - TB)/(TM - TB)*exp(omega*(TM - TB) - (TM - TC)/DeltaM))
    }
    
    if(Y <= 0.0) Y = 0.0
    
    return(Y)
  }
  
  # the lognormal probability density function
  LnormDist = function(x, mulog, sigmalog){
    
    y = 1/(x*sigmalog*sqrt(2*pi))*exp(-1/(2*sigmalog^2)*(log(x)-mulog)^2)
    
    return(y)
  }
  
  # This function does a convolution of two variables of the same length
  # but with padding on the right hand side
  ConvolveFunc = function(x1, y1, padsize){
    # padding
    #x2 = c(rep(0,padsize), x1, rep(0,padsize))
    #y2 = c(rep(0,padsize), y1, rep(0,padsize))
    x2 = c(x1, rep(0,padsize))
    y2 = c(y1, rep(0,padsize))
    
    # fast Fourier transforming
    fx2 = fft(x2)
    fy2 = fft(y2)
    
    # Doing the convolution and inverse transforming
    Convolution1 = Re(fft(fx2*fy2, inverse = T)/length(x2))
    
    # Now clipping
    #Convolution2 = Convolution1[(padsize+1):(length(x1)+padsize)]
    Convolution2 = Convolution1[1:length(x1)]
    
    # This funky step adds in all of the individuals that have developed beyond the next
    # threshold, but adds them in the last possible spot so that they do not effect
    # the development distribution behind them.
    Convolution2[length(x1)] = Convolution2[length(x1)] + 
      sum(Convolution1[(length(x1)+1):length(x2)])
    
    return(Convolution2)
  }
  
  # setting up the domain in physiological age
  avec = seq(1e-20, 2, length.out = 2^8) # domain for the larval stage
  da = avec[3] - avec[2]
  
  # Figuring out where in the domain avec = 1
  # (upper breakpoint for egg stage)
  # which.min(abs(avec - 1))
  # avec[128]
  # avec[129]
  
  # Defining all of the parameters ------
  
  # Defining the time step
  # We will assume that the reaction rate is measured per day
  # like the time step
  deltat = 1.0 # units are days
  
  # parameters for oviposition rate
  # (from Regniere et al 2012)
  sigma0 = 0.2458         # controls rate variability
  TB0 = 4.6341            # base temperature in degrees C
  DeltaB0 = 0.1
  TM0 = 27.7587           # max temperature in degree C
  DeltaM0 = 3.0759
  omega0 = 0.3684
  psi0 = 0.005199
  
  # parameters for development rate for eggs
  # (from Regniere et al 2012)
  sigma1 = 0.1799         # controls rate variability
  TB1 = 7.0               # base temperature in degrees C
  DeltaB1 = 0.019297569
  TM1 = 30.0928           # max temperature in degree C
  DeltaM1 = 4.4175
  omega1 = 0.2563
  psi1 = 0.02317
  
  # parameters for development rate for L1 larvae
  # (from Regniere et al 2012)
  sigma2 = 0.2911
  TB2 = 3.5559
  DeltaB2 = 0.1
  TM2 = 29.2647
  DeltaM2 = 3.8227
  omega2 = 0.2398
  psi2 = 0.01082
  
  # parameters for development rate for L2 larvae
  # (from Regniere et al 2012)
  sigma3 = 0.3799
  TB3 = 6.9598
  DeltaB3 = 0.097087379
  TM3 = 28.9047
  DeltaM3 = 3.0374
  omega3 = 0.3714
  psi3 = 0.01072
  
  # parameters for development rate for L3 larvae
  # (from Regniere et al 2012)
  sigma4 = 0.3868
  TB4 = 6.8462
  DeltaB4 = 0.1
  TM4 = 28.7013
  DeltaM4 = 2.5359
  omega4 = 0.4399
  psi4 = 0.003892
  
  # parameters for development rate for L4 larvae
  # (from Regniere et al 2012)
  sigma5 = 0.3932
  TB5 = 16.2464
  DeltaB5 = 0.039052889
  TM5 = 28.0
  DeltaM5 = 4.5504
  omega5 = 0.2593
  psi5 = 0.05034
  
  # parameters for development rate for pupae
  # (from Regniere et al 2012)
  sigma6 = 0.2998
  TB6 = 5.63
  DeltaB6 = 0.10989011
  TM6 = 28.55
  DeltaM6 = 2.86
  omega6 = 0.1532
  psi6 = 0.02054
  
  # parameters for development rate for teneral adults
  # (from Regniere et al 2012)
  sigma7 = 0.5284
  TB7 = 4.24
  DeltaB7 = 0.099967011
  TM7 = 35.0
  DeltaM7 = 7.1479
  omega7 = 0.1463
  psi7 = 0.01173
  
  # a container for predictions
  Fec = rep(0,length(Tmin))
  Fec[1:StartT] = rep(82.0, StartT)
  
  # containers for the other predictions
  Eggs = rep(0,length(Tmin))
  L1 = rep(0,length(Tmin))
  L2 = rep(0,length(Tmin))
  L3 = rep(0,length(Tmin))
  L4 = rep(0,length(Tmin))
  Pupae = rep(0,length(Tmin))
  TenAd = rep(0,length(Tmin))
  Adults = rep(0,length(Tmin))
  Flown = rep(0,length(Tmin))
  
  ##########################################
  # Initializing all of the stages ------
  ##########################################
  # Initializing the old eggs
  OldE = rep(0,length(avec))
  OE = rep(0,length(avec))
  NewEggstm1 = 0.0
  
  # Initializing the old L1 larvae
  OldL1 = rep(0,length(avec))
  OL1 = rep(0,length(avec))
  NewL1tm1 = 0.0
  
  # Initializing the old L2 larvae
  OldL2 = rep(0,length(avec))
  OL2 = rep(0,length(avec))
  NewL2tm1 = 0.0
  
  # Initializing the old L3 larvae
  OldL3 = rep(0,length(avec))
  OL3 = rep(0,length(avec))
  NewL3tm1 = 0.0
  
  # Initializing the old L4 larvae
  OldL4 = rep(0,length(avec))
  OL4 = rep(0,length(avec))
  NewL4tm1 = 0.0
  
  # Initializing the old pupae
  OldP = rep(0,length(avec))
  OP = rep(0,length(avec))
  NewPtm1 = 0.0
  
  # Initializing the old teneral adults
  OldT = rep(0,length(avec))
  OT = rep(0,length(avec))
  NewTtm1 = 0.0
  
  #################################################################################
  # The integrated IMAP (Insect Mortality and Phenology) model ------
  #################################################################################
  
  # Doing the iteration
  for(i in (StartT+1):length(Temp)){
    
    # Making an adjustment to minimum air temperature to account for the buffering
    # effect of tree bark as described by Bolstad, Bentz, and Logan (1997)
    Tmn = Tmin[i] + 1.8
    
    #############################################################
    # The two step process for oviposition
    #############################################################
    
    # Step 1: computing the development rate.
    # Development for oviposition (Eggs increase)
    Y1 = RegniereFunc(TC = Temp[i], TB0, DeltaB0, TM0, DeltaM0, omega0, psi0)
    
    # Aplying winter mortality to egg laying adults
    if(Tmn <= -18.0){
      Fec[i-1] = 0.0
    }
    
    # Step 2: Simulating oviposition 
    # (Fec represents the number of eggs remaining)
    Fec[i] = Fec[i-1]*exp(-Y1)
    
    ############################################################
    # The six step process for egg development
    ############################################################
    
    # Step 1: New eggs individuals that developed into eggs in this 
    # time step are individuals that did not stay in the fecundity stage.
    NewEggs = Fec[i-1]*(1 - exp(-Y1))
    
    # Step 2: computing the development rate.
    # Develpment for eggs
    Y2 = RegniereFunc(TC = Temp[i], TB1, DeltaB1, TM1, DeltaM1, omega1, psi1)
    
    # Step 3: To compute the number of new L1 larvae in the next stage,
    # we need to know how many old Eggs there were 
    # in the previous step. 
    OldE = OE
    
    # I use the -18 threshold to kill all eggs as
    # described in Regniere 2015
    if(Tmn <= -18.0){
      OldE = rep(0,length(avec))
      NewEggstm1 = 0.0
    }
    
    if(Y2 > 0.0){
      # Step 4 computing the Fourier transform of the aging kernel
      mu1 = log(Y2*deltat)
      G1 = LnormDist(x = avec, mulog = mu1, sigmalog = sigma1)
      G1 = G1/sum(G1)     # normalizing
      
      # Step 5: Doing the convolution
      # Inverse Fourier transforming
      OE = ConvolveFunc(x1 = OldE, y1 = G1, padsize = length(avec)) + NewEggstm1*G1
      NewEggstm1 = NewEggs
      
      # Step 6: Computing how eggs there are
      Eggs[i] = sum(na.omit(OE[1:128])) + NewEggs
      
    }else{
      OE = OldE
      NewEggstm1 = NewEggstm1 + NewEggs
      Eggs[i] = sum(na.omit(OE[1:128])) + NewEggstm1
    }
    
    ############################################################
    # The six step process for L1 
    ############################################################
    
    # Step 1: New L1 individuals that developed into L1 larvae in this 
    # time step are Egg individuals that exceeded the Egg breakpoint
    OldLarv1 = sum(OldE[129:length(avec)], na.rm = T)
    NewL1 = sum(OE[129:length(avec)], na.rm = T) - OldLarv1
    
    # Just to make sure that silly things don't happen
    if(NewL1 < 0.0) NewL1 = 0.0
    
    # Step 2: computing the development rate for L1 larvae
    Z2 = RegniereFunc(TC = Temp[i], TB2, DeltaB2, TM2, DeltaM2, omega2, psi2)
    
    # Step 3: To compute the number of new L2 larvae in the next stage, 
    # we need to know how many old L1 larvae there were 
    # in the previous step. 
    OldL1 = OL1
    
    if(Z2 > 0.0){
      # Step 4: computing the Fourier transform of the aging kernel
      mu2 = log(Z2*deltat)
      G2 = LnormDist(x = avec, mulog = mu2, sigmalog = sigma2)
      G2 = G2/sum(G2)     # normalizing
      
      # Step 5: Doing the convolution 
      # Inverse Fourier tranforming
      OL1 = ConvolveFunc(x1 = OldL1, y1 = G2, padsize = length(avec)) + NewL1tm1*G2
      NewL1tm1 = NewL1
      
      # Step 6: Computing how many L1 larvae there are
      L1[i] = sum(na.omit(OL1[1:128])) + NewL1
    }else{
      OL1 = OldL1
      NewL1tm1 = NewL1tm1 + NewL1
      
      L1[i] = sum(na.omit(OL1[1:128])) + NewL1tm1
    }
    
    ############################################################
    # The six step process for L2
    ############################################################
    
    # Step 1: New L2 individuals that developed into L2 larvae in this 
    # time step are L1 individuals that exceeded the L1 breakpoint
    OldLarv2 = sum(OldL1[129:length(avec)], na.rm = T)
    NewL2 = sum(OL1[129:length(avec)], na.rm = T) - OldLarv2
    
    # Just to make sure that silly things don't happen
    if(NewL2 < 0.0) NewL2 = 0.0
    
    # Step 2: computing the development rate for L2 larvae
    Z3 = RegniereFunc(TC = Temp[i], TB3, DeltaB3, TM3, DeltaM3, omega3, psi3)
    
    # Step 3: To compute the number of new L3 larvae in the next stage, 
    # we need to know how many Oldold and old L2 larvae there were 
    # in the previous step.
    OldL2 = OL2
    
    # Step 4: computing the Fourier transform of the aging kernel
    if(Z3 > 0.0){
      mu3 = log(Z3*deltat)
      G3 = LnormDist(x = avec, mulog = mu3, sigmalog = sigma3)
      G3 = G3/sum(G3)     # normalizing
      
      # step 5: Doing the covolution
      # Inverse Fourier transforming
      OL2 = ConvolveFunc(x1 = OldL2, y1 = G3, padsize = length(avec)) + NewL2tm1*G3
      NewL2tm1 = NewL2
      
      # Step 6: Computing how many L2 larvae there are
      L2[i] = sum(na.omit(OL2[1:128])) + NewL2
    }else{
      OL2 = OldL2
      NewL2tm1 = NewL2tm1 + NewL2
      L2[i] = sum(na.omit(OL2[1:128])) + NewL2tm1
    }
    
    ############################################################
    # The six step process for L3
    ############################################################
    
    # Step 1: New L3 individuals that developed into L3 larvae in this 
    # time step are L2 individuals that exceeded the L2 breakpoint 
    OldLarv3 = sum(OldL2[129:length(avec)], na.rm = T)
    NewL3 = sum(OL2[129:length(avec)], na.rm = T) - OldLarv3
    
    # Just to make sure that silly things don't happen
    if(NewL3 < 0.0) NewL3 = 0.0
    
    # Step 2: computing the development rate for L3 larvae
    Z4 = RegniereFunc(TC = Temp[i], TB4, DeltaB4, TM4, DeltaM4, omega4, psi4)
    
    # Step 3: To compute the number of new L4 larvae in the next stage, 
    # we need to know how many Oldold and old L3 larvae there were 
    # in the previous step.
    OldL3 = OL3
    
    if(Z4 > 0.0){
      # Step 4: computing the Fourier transform of the aging kernel
      mu4 = log(Z4*deltat)
      G4 = LnormDist(x = avec, mulog = mu4, sigmalog = sigma4)
      G4 = G4/sum(G4)     # normalizing
      
      # Step 5: Doing the convolution
      # Inverse Fourier transforming
      OL3 = ConvolveFunc(x1 = OldL3, y1 = G4, padsize = length(avec)) + NewL3tm1*G4
      NewL3tm1 = NewL3
      
      # Step 6: Computing how many L3 larvae there are
      L3[i] = sum(na.omit(OL3[1:128])) + NewL3
    }else{
      OL3 = OldL3
      NewL3tm1 = NewL3tm1 + NewL3
      L3[i] = sum(na.omit(OL3[1:128])) + NewL3tm1
    }
    
    ############################################################
    # The six step process for L4
    ############################################################
    
    # Step 1: New L4 individuals that developed into L4 larvae in this 
    # time step are L3 individuals that exceeded the L3 breakpoint
    OldLarv4 = sum(OldL3[129:length(avec)], na.rm = T)
    NewL4 = sum(OL3[129:length(avec)], na.rm = T) - OldLarv4
    
    # Just to make sure that silly things don't happen
    if(NewL4 < 0.0) NewL4 = 0.0
    
    # Step 2: computing the development rate for L4 larvae
    Z5 = RegniereFunc(TC = Temp[i], TB5, DeltaB5, TM5, DeltaM5, omega5, psi5)
    
    # Step 3: To compute the number of new Pupae in the next stage, 
    # we need to know how many Oldold and old L4 larvae there were 
    # in the previous step. 
    OldL4 = OL4
    
    if(Z5 > 0.0){
      # Step 4: computing the Fourier transform of the aging kernel
      mu5 = log(Z5*deltat)
      G5 = LnormDist(x = avec, mulog = mu5, sigmalog = sigma5)
      G5 = G5/sum(G5)     # normalizing
      
      # Step 5: Doing the convolution
      # Inverse Fourier transforming
      OL4 = ConvolveFunc(x1 = OldL4, y1 = G5, padsize = length(avec)) + NewL4tm1*G5
      NewL4tm1 = NewL4
      
      # Step 6: Computing how many L4 larvae there are
      L4[i] = sum(na.omit(OL4[1:128])) + NewL4
    }else{
      OL4 = OldL4
      NewL4tm1 = NewL4tm1 + NewL4
      L4[i] = sum(na.omit(OL4[1:128])) + NewL4tm1
    }
    
    ############################################################
    # The six step process for Pupae
    ############################################################
    
    # Step 1: New pupae that developed into pupae in this 
    # time step are L4 individuals that exceeded the L4 breakpoint 
    OldPup = sum(OldL4[129:length(avec)], na.rm = T)
    NewP =sum(OL4[129:length(avec)], na.rm = T) - OldPup
    
    # Just to make sure that silly things don't happen
    if(NewP < 0.0) NewP = 0.0
    
    # Step 2: computing the development rate for pupae
    Z6 = RegniereFunc(TC = Temp[i], TB6, DeltaB6, TM6, DeltaM6, omega6, psi6)
    
    # Step 3: To compute the number of new teneral adults in the next stage, 
    # we need to know how many Oldold and old Pupae there were 
    # in the previous step. 
    OldP = OP
    
    # I use the -18 threshold to kill all pupae as
    # described in Regniere 2015
    if(Tmn <= -18.0){
      OldP = rep(0,length(avec))
      NewPtm1 = 0.0
    }
    
    # Step 4: computing the Fourier transform of the aging kernel
    if(Z6 > 0.0){
      mu6 = log(Z6*deltat)
      G6 = LnormDist(x = avec, mulog = mu6, sigmalog = sigma6)
      G6 = G6/sum(G6)     # normalizing
      
      # Step 5: Doing the convolution
      # Inverse Fourier transforming
      OP = ConvolveFunc(x1 = OldP, y1 = G6, padsize = length(avec)) + NewPtm1*G6
      NewPtm1 = NewP
      
      # Step 6: Tallying up pupae
      Pupae[i] = sum(na.omit(OP[1:128])) + NewP
    }else{
      OP = OldP
      NewPtm1 = NewPtm1 + NewP
      Pupae[i] = sum(na.omit(OP[1:128])) + NewPtm1
    }
    
    ############################################################
    # The six step process for teneral adults
    ############################################################
    
    # Step 1: New teneral adults that developed into tenerals in this 
    # time step are pupae that exceeded the pupal breakpoint 
    OldT = sum(OldP[129:length(avec)], na.rm = T)
    NewT = sum(OP[129:length(avec)], na.rm = T) - OldT
    
    # Just to make sure that silly things don't happen
    if(NewT < 0.0) NewT = 0.0
    
    # Step 2: computing the development rate for teneral adults
    Z7 = RegniereFunc(TC = Temp[i], TB7, DeltaB7, TM7, DeltaM7, omega7, psi7)
    
    # Step 3: To compute the number of new adults (Step 1.8), 
    # we need to know the number of Oldold and old Teneral adults 
    # in the previous step. 
    OldT = OT 
    
    # I use the -18 threshold to kill all teneral adults as
    # described in Regniere 2015
    if(Tmn <= -18.0){
      OldT = rep(0,length(avec))
      NewTtm1 = 0.0
    }
    
    # Step 4: computing the Fourier transform of the aging kernel
    if(Z7 > 0.0){
      mu7 = log(Z7*deltat)
      G7 = LnormDist(x = avec, mulog = mu7, sigmalog = sigma7)
      G7 = G7/sum(G7)     # normalizing
      
      # Step 5: Doing the convolution
      # Inverse Fourier transforming
      OT = ConvolveFunc(x1 = OldT, y1 = G7, padsize = length(avec)) + NewTtm1*G7
      NewTtm1 = NewT
      
      # Step 6: Tallying up teneral adults
      TenAd[i] = sum(na.omit(OT[1:128])) + NewT
    }else{
      OT = OldT
      NewTtm1 = NewTtm1 + NewT
      TenAd[i] = sum(na.omit(OT[1:128])) + NewTtm1
    }
    
    ############################################################
    # The two step process for adults
    ############################################################
    
    # Step 1: New adults that developed into adults in this 
    # time step are teneral individuals that exceeded the teneral 
    # breakpoint in this time step
    OldAds = sum(OldT[129:length(avec)], na.rm = T)
    NewA = sum(OT[129:length(avec)], na.rm = T) - OldAds
    
    # Just to make sure that silly things don't happen
    if(NewA < 0.0) NewA = 0.0
    
    # To kill adult beetles, I use the -18 temperature used
    # by Regniere et al 2015.
    if(Tmn <= -18.0) Adults[i-1] = 0.0
    
    # Counting up adults
    Adults[i] = Adults[i-1] + NewA
  }
  
  ###############################################
  ##Now all of the phenology has been simulated##
  ###############################################
  LifeCycleMat = matrix(c(Fec, Eggs, L1, L2, L3, L4, Pupae, TenAd, Adults), 
                        nrow = 9, ncol = 443, byrow = TRUE)
  
  return(LifeCycleMat)
}

# Now I will test the Imap function
# Reading in some temperature data (in degrees C)
Tmin = read.csv("C:\\Users\\Grace.Veenstra\\Documents\\GitHub\\Goodsman_IPM\\ece33590-sup-0004-supinfo1\\JasperDailyMin2.csv", h = F)   # Min temp  
Tmax = read.csv("C:\\Users\\Grace.Veenstra\\Documents\\GitHub\\Goodsman_IPM\\ece33590-sup-0004-supinfo1\\JasperDailyMax2.csv", h = F)   # Max temp

Tmin = as.numeric(Tmin$V1)
Tmax = as.numeric(Tmax$V1)

ptm <- proc.time()
Out1 = ImapFuncNewJ(Tmin = Tmin, Tmax = Tmax, StartT = 15)
proc.time() - ptm

Times = 1:443
plot(Out1[1,] ~ Times, type = 'l')
lines(Out1[2,] ~ Times, lty = 2)
lines(Out1[3,] ~ Times, lty = 3)
lines(Out1[4,] ~ Times, lty = 1, col = 'blue')
lines(Out1[5,] ~ Times, lty = 1, col = 'orange')
lines(Out1[6,] ~ Times, lty = 1, col = 'green')
lines(Out1[7,] ~ Times, lty = 2, col = 'purple')
lines(Out1[8,] ~ Times, lty = 2, col = 'grey')
lines(Out1[9,] ~ Times, lty = 1, col = 'red')
colSums(Out1)
# Yes; it works.


plot(Eggs, type='l', xlim=c(0,30))
