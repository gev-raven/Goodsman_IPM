# Workshopping the development function


# Reading in the temperature data (in degrees C)
temp_df <- read.csv("C:\\Users\\Grace.Veenstra\\Documents\\GitHub\\Goodsman_IPM\\Egg Modeling Test\\MarshCreek_Temp_2013.csv",
                    header=TRUE) #takes csv file and creates dataframe
temp_df <- temp_df[c("Date","HUC_10","Mean")] #removes 'min' 'max' and 'HUC_8' columns from frame
temp_df <- temp_df[which(temp_df$HUC_10 %in% "1706020503"),] #filters to only include Marsh Creek ("1706020503" is creek identifier)
temp_df$Ydate <- yday(temp_df$Date) #convert date to 'day of year' (value 1 to 365)
temp_df$Jdate <- temp_df$Date
temp_df$Jdate <- as.numeric(as.Date(temp_df$Jdate, origin="1993-01-01")) #converts date to 'days since 1993-1-1'

temp <- temp_df$Mean #assigns variable to the mean temp in data
temp <- as.numeric(temp)


############################################
# Development model -----
############################################

### Parameters ----

# defining time step
deltat = 1.0 #units are days

#initial number of eggs laid
imax = 100

# Parameters for egg development rate
Tp <- temp
a_egg <- exp(10.404)
b_egg <-- 2.043
c_egg <-- 7.575


### the development function for temperature dependent egg development
development_func = function(Tp, a, b, c) {
  
  # develop_time is the development time (in days)
  # Tp is mean daily temp;
  
  develop_time <- exp(log(a)+log((Tp-c)^b))
  
  develop_df <- temp_df
  
  develop_df$DailyDevelopment <- (1/develop_time) #daily development rate
  
  
  #define development period
  
  timetoemerge <- NULL
  startspawn <- 2456505 #aug 1, 2013 #aug 1 2008, ordinal date : 2454680
  endspawn <- startspawn + 40
  
  for (spawndate in startspawn:endspawn) {
    
    develop_df$DevelopmentPeriod <- subset(develop_df$Jdate, 
                                           develop_df$Jdate <= spawndate & temp_df$Jdate <= (startspawn + 366))
    ## DevelopmentPeriod is a vector with *just* ordinal dates
    
    develop_df$TotalDevelopment <- cumsum(develop_df$DailyDevelopment) 
    ## total development = sum(dailydevelopment)
    
    emerge <- min(which((cumsum(develop_df$DailyDevelopment)) >= 1)) #once y = 1, it is emergence time
    #in integration, warning: no non missing arguments
    
    timetoemerge <- rbind(timetoemerge,emerge)
  }
  
}


############################################
# Integrated model ------
############################################

# The start time
StartT = 15

# Containers for predictions
Eggs = rep(0,length(temp))


ptm <- proc.time()

# Doing the iteration
for(i in (StartT+1):length(temp)){
  
  # Computing the median development rate for each life stage in this time step
  median1 = development_func(Tp = temp[i], a_egg, b_egg, c_egg)# for eggs
  
  # The mu parameter of the lognormal distribution is given by ln(R[T]*deltat)
  mu1 = log(median1*deltat) # for eggs
  
  # In each time step we set the counters back to zero
  EggCount = 0.0
  
  # Now we slot in the observed number of individuals in each life stage
  # at this time step (corresponding to i)
  Eggs[i] = EggCount
  
}
proc.time() - ptm

times = 1:length(temp)

# Now I write the data to csv (20 times)
IBM = data.frame(times, Eggs)
write.csv(IBM, "C:\\Users\\Grace.Veenstra\\Documents\\GitHub\\Goodsman_IPM\\Egg Modeling Test\\PhenSimR.csv")

graphics.off()
plot(Eggs ~ times, type = 'l', ylim = c(0,100))

plot(develop_df$DailyDevelopment ~ times, type = 'l', ylim=c(0,100))


