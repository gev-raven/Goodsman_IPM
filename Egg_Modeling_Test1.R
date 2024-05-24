# Uses code modified from Goodsman 2018 R code "IMAPfunctionNewJ2v2.R". This
# code attempts to simulate chinook salmon egg development and
# emergence timing.

# Development rate equation from Abbie and Beacham and Murray Model
# Temperature data for Upper Middle Fork Salmon River at Marsh Creek.

# Model assumes no egg mortality.

#####################################################################

# The ImapFunc function takes temperature as input (Temp)

# The ImapFuncNewJ function returns as output a matrix
# with rows that comprise vectors of life stages
# and columns as individual days.


####################################################################
###### IPM Model Start ----- 
####################################################################

ImapFuncNewJ = function(Tmin, Tmax, StartT){

  ## Components of Model -----
  ############################

  # temperature in degrees Celsius
  
  Temp

  ## Beacham & Murray Egg Development Model
  ## Predicts emergence timing of Chinook Salmon in Salmon River
  # Notes: Constants from 1990 paper, spawn date assumed August 1st
    
  #### NEEDS MODIFICATION FOR "BearValley..."
    
  development.func = function(Tp, a, b, c) {

    # develop.time is the development time (in days)

    develop.time <- exp(log(a)+log((Tp-c)^b))
    
    daily.develop <- (1/develop.time) #daily development rate
        #### BearValleyElkCreekTemperaturedaily$dailydevelopment<-(1/develop.time) #daily development rate

    
    #define development period

    timetoemerge <- NULL
    startspawn <- 2454680
    endspawn <- startspawn + 40
  
    for (spawndate in startspawn:endspawn) {
    
      DevelopmentPeriod <- subset(BearValleyElkCreekTemperaturedaily, 
      Ydate >= spawndate & Ydate <= (startspawn+366))
    
      DevelopmentPeriod$totaldevelopment <- cumsum(DevelopmentPeriod$dailydevelopment)
    
      y <- min(which((cumsum(DevelopmentPeriod$dailydevelopment)) >= 1))
    
      timetoemerge <- rbind(timetoemerge,y)
    }

  }
  
  
  ## the lognormal probability density function
  
  LnormDist = function(x, mulog, sigmalog){
    
    y = 1/(x*sigmalog*sqrt(2*pi))*exp(-1/(2*sigmalog^2)*(log(x)-mulog)^2)
    
    return(y)
  }

  
  ## Defining the parameters -----
  ################################
  
  # Parameters for egg development

    Tp <- BearValleyElkCreekTemperaturedaily$Temperature
    a_egg <- exp(10.404)
    b_egg <-- 2.043
    c_egg <-- 7.575

    
  ####################################################################
  ##### Integrated Model -----
  ####################################################################
    
    # Doing the iteration
    for(i in (StartT+1):length(temp)){
    
      ### Process of Egg Development -----
      ####################################
      
      # Step 2: Computing the development rate for eggs
      Y2 = development.func(Tp = , a_egg, b_egg, c_egg)
      
      if(Y2 > 0.0){
        # Step 4 computing the Fourier transform of the aging kernel
        mu1 = log(Y2*deltat)
        G1 = LnormDist(x = avec, mulog = mu1, sigmalog = sigma1) #probability density function
        G1 = G1/sum(G1)     # normalizing
        
        # Step 5: Doing the convolution
        # Inverse Fourier transforming
        OE = ConvolveFunc(x1 = OldE, y1 = G1, padsize = length(avec)) + NewEggstm1*G1
        NewEggstm1 = NewEggs
        
        # Step 6: Computing how many eggs there are
        Eggs[i] = sum(na.omit(OE[1:128])) + NewEggs
        
      }else{
        OE = OldE
        NewEggstm1 = NewEggstm1 + NewEggs
        Eggs[i] = sum(na.omit(OE[1:128])) + NewEggstm1
      }
      
      ### Others ~
        
    }
  
    ### Setup the IPM Matrix Output
    
    LifeCycleMat = matrix(c(Eggs),
                          nrow = 1, ncol = 443, byrow = TRUE)

    return(LifeCycleMat)
}   
  
####################################################################
##### IPM Model Test -----
####################################################################


temp_df = read.csv("", header=TRUE) #takes csv file
temp_df_cols <- select(temp_df, HUC_10, Mean) #filters columns
filter(temp_df_cols, HUC_10 == "1706020503") # filters stream identifier to Marsh
temp <- temp_df_cols["Mean"]


