

## TEMP ----

Tmin = read.csv("C:/Users/Grace.Veenstra/Downloads/ece33590-sup-0004-supinfo/JasperDailyMin2.csv", h = F)   # Min temp  
Tmax = read.csv("C:/Users/Grace.Veenstra/Downloads/ece33590-sup-0004-supinfo/JasperDailyMax2.csv", h = F)   # Max temp

Tmin = as.numeric(Tmin$V1)
Tmax = as.numeric(Tmax$V1)

Temp = 0.5*(Tmax + Tmin) + 0.9 + 6.6*(Tmax - Tmin)/(2*24.4)
Tmax2 = Tmax + 6.6*(Tmax - Tmin)/24.4
Tmin2 = Tmin + 1.8

## Function ----

RegniereFunc = function(TC, TB, DeltaB, TM, DeltaM, omega, psi){
  
  Y = 0.0
  
  if(TC >= TB & TC <= TM){
    Y = psi*(exp(omega*(TC - TB)) - (TM - TC)/(TM - TB)*exp(-omega*(TC - TB)/DeltaB)
             - (TC - TB)/(TM - TB)*exp(omega*(TM - TB) - (TM - TC)/DeltaM))
  }
  
  if(Y <= 0.0) Y = 0.0
  
  return(Y)
}

##Parameters ----

# Defining the time step
# We will assume that the reaction rate is measured per day
# like the time step
deltat = 1.0 # units are days

# containers for the other predictions
Eggs = rep(0,length(Tmin))

# parameters for development rate for eggs
# (from Regniere et al 2012)
sigma1 = 0.1799         # controls rate variability
TB1 = 7.0               # base temperature in degrees C
DeltaB1 = 0.019297569
TM1 = 30.0928           # max temperature in degree C
DeltaM1 = 4.4175
omega1 = 0.2563
psi1 = 0.02317

Y=NULL

##temp??
t.mean <- mean(Temp)
t.sd <- sd(Temp)
n.days <- length(Tmin)
t.dist <- rnorm(n = n.days, mean = t.mean, sd = t.sd)
t.dist2 <- round(t.dist)
t.dist2
hist(t.dist2)

for(i in 1:length(Tmin)){
  
  # Making an adjustment to minimum air temperature to account for the buffering
  # effect of tree bark as described by Bolstad, Bentz, and Logan (1997)
  Tmn = Tmin[i] + 1.8

  Y[i] = RegniereFunc(TC = Temp[i], TB1, DeltaB1, TM1, DeltaM1, omega1, psi1)

}

Y

times <- 1:443

plot(Y ~ times)
plot(Temp ~ times)
plot(Y ~ Temp)
plot(Y ~ num.day)

num.day <- as.vector(table(t.dist2)) 
num.day

test <- rep(Y, num.day) 
hist(test)

qqnorm(log(test))
qqline(log(test))
