# Workshopping the development function


##############################################

# Reading in the temperature data (in degrees C)
temp_df <- read.csv("C:\\Users\\Grace.Veenstra\\Documents\\GitHub\\Goodsman_IPM\\Egg Modeling Test\\MarshCreek_Temp_2015-2016.csv",
                    header=TRUE) #takes csv file and creates dataframe
temp_df <- temp_df[c("Date","HUC_10","Mean")] #removes 'min' 'max' and 'HUC_8' columns from frame
temp_df <- temp_df[which(temp_df$HUC_10 %in% "1706020503"),] #filters to only include Marsh Creek ("1706020503" is creek identifier)
temp_df$Date <- format.Date(temp_df$Date, format="%Y-%m-%d")
temp_df$Ydate <- yday(temp_df$Date) #convert date to 'day of year' (value 1 to 365)
temp_df$Jdate <- temp_df$Date

#temp_df$Jdate <- as.numeric(as.Date(temp_df$Jdate, origin="1993-01-01")) #converts date to 'days since 1993-1-1'
### still trying to figure out ordinal

temp <- temp_df$Mean #assigns variable to the mean temp in data
temp <- as.numeric(temp)


############################################
# Development model -----
############################################

Tp <- temp
a <- exp(10.404)
b <- -2.043
c <- -7.575
timetoemerge <- NULL

spawn_start <- 213 #julian
spawn_var <- 5 #spawn variance is 5 days
spawndate_lower <- spawn_start - spawn_var
spawndate_upper <- spawn_start + spawn_var

for(i in spawndate_lower:spawndate_upper) {
  
  i = spawndate_lower
  s <- i
  startspawn <- s
  endspawn <- startspawn + 40
 
  development_func = function(Tp, a, b, c, startspawn, endspawn) {
      
  # develop_time is the development time (in days)
  # Tp is mean daily temperature
      
  develop_time <- exp(log(a)+log((Tp-c)^b))
  DailyDevelopment <- (1/develop_time) #daily development rate
      
  #define development period
      
  for (s in startspawn:endspawn) {
        
    DevelopmentPeriod <- subset(temp_df$Ydate, 
                                temp_df$Ydate >= spawndate & temp_df$Ydate <= (startspawn + 366))
      ## DevelopmentPeriod is a vector with *just* ordinal dates
        
    TotalDevelopment <- cumsum(DailyDevelopment)
      ## total development = sum(dailydevelopment)
        
    y <- min(which((cumsum(DailyDevelopment)) >= 1)) #once y = 1, it is emergence time
      #in integration, warning: no non missing arguments
        
    timetoemerge <- rbind(timetoemerge, y)
        
        return(TotalDevelopment)
        return(DevelopmentPeriod)
        return(y)
      }
      
    }
    
  }
  
  i = i + 1
    
}


### the development function for temperature dependent egg development


times = 1:length(temp)
plot(DailyDevelopment ~ times, type='l')
plot(TotalDevelopment ~ times, type='l', ylim=c(0,1))

