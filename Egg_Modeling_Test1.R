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
###### Packages Set Up ----- 
####################################################################

install.packages("lubridate")
library(lubridate)
install.packages("anytime")
library(anytime)
install.packages("tidyverse")
library(tidyverse)
install.packages("dbplyr")
library(dbplyr)
install.packages("dplyr")
library(dplyr)
install.packages("dtplyr")
library(dtplyr)
install.packages("insol") #insol package no longer available
library(insol)


####################################################################
###### Temperature Set Up ----- 
####################################################################

temp_df <- read.csv("C:\\Users\\Grace.Veenstra\\Documents\\GitHub\\Goodsman_IPM\\Egg Modeling Test\\MarshCreek_Temp_2013.csv",
                    header=TRUE) #takes csv file and creates dataframe
temp_df <- temp_df[c("Date","HUC_10","Mean")] #removes 'min' 'max' and 'HUC_8' columns from frame
temp_df <- temp_df[which(temp_df$HUC_10 %in% "1706020503"),] 
                  #filters to only include Marsh Creek ("1706020503" is creek identifier)
temp_df$Ydate <- yday(temp_df$Date) #convert date to day of year (value 1 to 365)
temp_df$Ydate
temp_df$Jdate <- JD(temp_df$Date) ## Doesn't work. JD not available
temp_df$Jdate

temp <- temp_df$Mean #assigns variable temp to the mean temperatures
temp <- as.numeric(temp)


####################################################################
###### IPM Model Start ----- 
####################################################################

ImapFuncNewJ = function(Tmin, Tmax, StartT){

  ## Components of Model -----
  ############################

  # temperature in degrees Celsius
  
  temp

  ## the Beacham & Murray Egg Development Model
  ## Predicts emergence timing of Chinook Salmon in Salmon River
  # Notes: Constants from 1990 paper, spawn date assumed August 1st
    
  #### NEEDS MODIFICATION FOR "BearValley..."
    
  development_func = function(Tp, a, b, c) {

    # develop.time is the development time (in days)
    # Tp is mean daily temp;

    develop_time <- exp(log(a)+log((Tp-c)^b))
    
    daily_develop <- (1/develop_time) #daily development rate
        frame <- temp_df
        frame$DailyDevelopment < daily_develop
        #### BearValleyElkCreekTemperaturedaily$dailydevelopment<-(1/develop.time) #daily development rate

    
    #define development period

    timetoemerge <- NULL
    startspawn <- 2454680 #set to august 1st, change format?
    endspawn <- startspawn + 40
  
    for (spawndate in startspawn:endspawn) {
    
      DevelopmentPeriod <- subset(frame, #development period made a subset
      Ydate >= spawndate & Ydate <= (startspawn+366)) #year set to count up to a year from spawning
    
      DevelopmentPeriod$TotalDevelopment <- cumsum(DevelopmentPeriod$DailyDevelopment) #total development = sum(dailydevelopment)
    
      y <- min(which((cumsum(DevelopmentPeriod$dailydevelopment)) >= 1)) #once y hits 1 is emergence
    
      timetoemerge <- rbind(timetoemerge,y)
    }

  }

  # the lognormal probability density function
  LnormDist = function(x, mulog, sigmalog){
    
    y = 1/(x*sigmalog*sqrt(2*pi))*exp(-1/(2*sigmalog^2)*(log(x)-mulog)^2)
    
    return(y)
  }
  

  ## Defining the parameters -----
  ################################
  
  # defining time step
  deltat = 1.0 #units are days
  
  #initial number of eggs laid
  imax = 100
  
  # Parameters for egg development rate
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

    }
  
    ### Setup the IPM Matrix Output
    
    LifeCycleMat = matrix(c(Eggs),
                          nrow = 1, ncol = 365, byrow = TRUE)

    return(LifeCycleMat)
}   
  
####################################################################
##### IPM Model Test -----
####################################################################



ptm <- proc.time()
Freq_Egg = ImapFuncNewJ(Tmin = Tmin, StartT = 15)
proc.time() - ptm

Times = 1:360
plot(Out1[1,] ~ Times, type = 'l')
lines(Out1[2,] ~ Times, lty = 2)

colSums(Out1)

