
# This code creates a distribution of egg emergence days for Chinook salmon
# based on the spawn date, with variance applied to the spawn date.
# The development rate equation comes from Abbie, and from the Beacham and Murray Model

# Temperature data is for Upper Middle Fork Salmon River at Marsh Creek.
# Model assumes no egg mortality.


##############################################
###### Packages Set Up ----- 
##############################################

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


##############################################

## Reading in the temperature data (in degrees C)
temp_df <- read.csv("C:\\Users\\Grace.Veenstra\\Documents\\GitHub\\Goodsman_IPM\\Egg Modeling Test\\MarshCreek_Temp_2015-2016.csv",
                    header=TRUE) #takes csv file and creates dataframe
temp_df <- temp_df[c("Date","HUC_10","Mean")] 
                    #removes 'min' 'max' and 'HUC_8' columns from frame
temp_df <- temp_df[which(temp_df$HUC_10 %in% "1706020503"),] 
                    #filters to only get Marsh Creek (Creek Identifier: "1706020503")
temp_df$Date <- parse_date_time(temp_df$Date, orders=c('mdy','ymd')) 
                    #reads different date formats
temp_df$Date <- format.Date(temp_df$Date, format="%Y-%m-%d")
                    #sets the date format of Date list
temp_df$Ydate <- yday(temp_df$Date) 
                    #convert date to 'day of year' (value 1 to 365)
temp_df$Jdate <- 1:nrow(temp_df) 
                    #sets multi-year 'julian' values (1-730+)

temp <- temp_df$Mean #assigns variable to the mean temp in data
temp <- as.numeric(temp)


############################################
###### Development model -----
############################################

development_func = function(Tp, a, b, c, startspawn) {
  
  # develop_time is the development time (in days)
  # Tp is daily mean temperature
  
  develop_df <- temp_df
  develop_time <- exp(log(a)+log((Tp-c)^b))
  develop_df$DailyDevelopment <- (1/develop_time) #daily development rate
  
  #define development period
  
  develop_df <- subset(develop_df, 
                       develop_df$Jdate >= startspawn & develop_df$Jdate <= (startspawn + 366))
  ## restricts develop_df to the period of development
  
  develop_df$TotalDevelopment <- cumsum(develop_df$DailyDevelopment)
  ## total development = sum of daily development
  
  y <- develop_df[min(which((develop_df$TotalDevelopment) >= 1)), "Jdate"] 
  ## once development = 1, it is emergence time and we extract emergence date
  
  return(y)
}

## Parameters

Tp <- temp
a <- exp(10.404)
b <- -2.043
c <- -7.575

## Spawn Start Range
spawn_start <- 213 #aug 1; mean spawn start date
spawn_var <- 5 #spawn start variance is 5 days
spawndate_lower <- spawn_start - spawn_var #lower limit of spawn start range
spawndate_upper <- spawn_start + spawn_var #upper limit of spawn start range
spawnstart.range <- spawndate_lower:spawndate_upper

i = 1
develop_df <- temp_df
emergence <- NULL


## Iteration over Spawn Range

for(i in 1:length(spawnstart.range)) { #runs loop for the range of spawning start dates
  
  startspawn <- spawnstart.range[i] #start of spawning is somewhere in the spawn_start range
  endspawn <- startspawn + 40 #define spawning window as 40 days
  
  for(spawndate in startspawn:endspawn) { #runs loop over the spawn window
    
    emergence[i] <- development_func(Tp = Tp, a = a, b = b, c = c, startspawn = startspawn)
    #collects the emergence days into a vector
    
  }
}


## Graphing

emergence

mu1 = mean(emergence)
sigma1 = sd(emergence)
emerge_freq <- dnorm(x = develop_df$Jdate, mean = mu1, sd = sigma1)
plot(emerge_freq, xlim=c(min(emergence), max(emergence)))
  
  ## so with a start date of aug 1 (213 julian), mean emergence is march 5 (430 julian)
  ## this correct?


plot(develop_df$TotalDevelopment ~ develop_df$Jdate, type='l', ylim=c(0,1))
