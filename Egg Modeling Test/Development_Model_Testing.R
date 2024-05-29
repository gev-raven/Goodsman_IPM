# Workshopping the development function

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

i = spawndate_lower
emergence <- NULL


for(i in spawndate_lower:spawndate_upper) { 
  
  s <- i
  startspawn <- s
  endspawn <- startspawn + 40
 
  development_func = function(Tp, a, b, c, startspawn, endspawn) {
    
    # develop_time is the development time (in days)
    # Tp is daily mean temperature
    
    develop_df <- temp_df
    develop_time <- exp(log(a)+log((Tp-c)^b))
    develop_df$DailyDevelopment <- (1/develop_time) #daily development rate
    
    #define development period
    
    for (s in startspawn:endspawn) {
        
      develop_df <- subset(develop_df, 
                           develop_df$Jdate >= startspawn & develop_df$Jdate <= (startspawn + 366))
        ## restricts develop_df to the period of development
      
      develop_df$TotalDevelopment <- cumsum(develop_df$DailyDevelopment)
        ## total development = sum of daily development
      
      y <- develop_df[min(which((develop_df$TotalDevelopment) >= 1)), "Jdate"] 
        ## once development = 1, it is emergence time and we extract emergence date
      
      emergence[i] = y
      
    }

  }
  
  i <- i + 1 #will count up for each iteration in the spawn_start range
  
}
  
development_func()
emergence

plot(develop_df$TotalDevelopment ~ develop_df$Jdate, type='l', ylim=c(0,1))
    
plot(emergence)


### the development function for temperature dependent egg development


plot(DailyDevelopment ~ times, type='l')
plot(TotalDevelopment ~ times, type='l', ylim=c(0,1))

