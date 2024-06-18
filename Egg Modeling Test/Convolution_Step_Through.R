develop_rate = 0.08
mu1 = log(develop_rate*deltat)
egg.dist = LnormDist(x = avec, mulog = mu1, sigmalog = sigma1) #probability density distribution of the development rate
egg.dist = egg.dist/sum(egg.dist) # normalizing
#Egg.B = ConvolveFunc(x1 = PrevEgg, y1 = egg.dist, padsize = length(avec)) + NewEggsT1*egg.dist
#plot(avec[1:10],egg.dist[1:10],main=paste("develop_rate",develop_rate))
#lines(avec[1:10],egg.dist[1:10],main=paste("develop_rate",develop_rate))


Fec = rep(0,t)
Eggs = rep(0, t)
develop <- rep(0,t)

# Initializing the previous time step eggs
PrevEgg = rep(0,length(avec))
Egg.B = rep(0,length(avec))
NewEggsT1 = 10.0


plot(avec[1:10],Egg.B[1:10],main=paste("develop_rate",develop_rate),ylim=c(0,10))
lines(avec[1:10],egg.dist[1:10],col="grey")
  Egg.B = ConvolveFunc(x1 = PrevEgg, y1 = egg.dist, padsize = length(avec)) #+ NewEggsT1*egg.dist
lines(avec[1:10],Egg.B[1:10],col="red")
  PrevEgg = Egg.B
  Egg.B = ConvolveFunc(x1 = PrevEgg, y1 = egg.dist, padsize = length(avec)) + NewEggsT1*egg.dist
lines(avec[1:10],Egg.B[1:10],col="green")
  PrevEgg = Egg.B
  Egg.B = ConvolveFunc(x1 = PrevEgg, y1 = egg.dist, padsize = length(avec)) #+ NewEggsT1*egg.dist
lines(avec[1:10],Egg.B[1:10],col="blue")


Egg.B[1:10]
PrevEgg = Egg.B


####

sigma2 = 0.1
growth_step = 0.008
mu2 = log(growth_step*deltat)
juv.growth.dist = LnormDist(x = avec2, mulog = mu2, sigmalog = sigma2)
juv.growth.dist = juv.growth.dist/sum(juv.growth.dist)
#Egg.B = ConvolveFunc(x1 = PrevEgg, y1 = egg.dist, padsize = length(avec)) + NewEggsT1*egg.dist
#plot(avec[1:10],egg.dist[1:10],main=paste("develop_rate",develop_rate))
#lines(avec[1:10],egg.dist[1:10],main=paste("develop_rate",develop_rate))


Fec = rep(0,t)
Eggs = rep(0, t)
develop <- rep(0,t)

# Initializing the previous time step
PreJuvWeight = rep(0,length(avec2))
JuvWeight = rep(0,length(avec2))
NewJuvT1 = 10

juv.growth.dist

i=1
plot(avec2[1:20],JuvWeight[1:20],main=paste("growth_rate",growth_step),ylim=c(0,5))
lines(avec2[1:20],juv.growth.dist[1:20],col="grey")

JuvWeight = ConvolveFunc(x1 = PreJuvWeight, y1 = juv.growth.dist, padsize = length(avec2)) #+ NewEggsT1*egg.dist
#JuvWeight = abs((JuvWeight))^(3)
lines(avec2[1:20],JuvWeight[1:20],col="red")
i=i+1

PreJuvWeight = (JuvWeight)#^(e1)
JuvWeight = ConvolveFunc(x1 = PreJuvWeight, y1 = juv.growth.dist, padsize = length(avec2)) 
#JuvWeight = abs((JuvWeight))^(3)
JuvWeight = JuvWeight + NewJuvT1*juv.growth.dist
lines(avec2[1:20],JuvWeight[1:20],col="green")
i=i+1

PreJuvWeight = (JuvWeight)#^(e1)
JuvWeight = ConvolveFunc(x1 = PreJuvWeight, y1 = juv.growth.dist, padsize = length(avec2)) #+ NewEggsT1*egg.dist
#JuvWeight = abs((JuvWeight))^(3)
lines(avec2[1:20],JuvWeight[1:20],col="blue")
i=i+1

  








############

i=500
omega <- JuvGrowthFunc(Tp=temp[i], a1, b1, c1, d1, e1)
omega <- omega*(1/e1)*100

sigma2 = 0.1
growth_rate = omega
mu2 = log(growth_rate*deltat)
juv.growth.dist = LnormDist(x = avec2, mulog = mu2, sigmalog = sigma2)
juv.growth.dist = juv.growth.dist/sum(juv.growth.dist)

# Initializing the previous time step
PreJuvWeight = rep(0,length(avec2))
JuvWeight = rep(0,length(avec2))
NewJuvT1 = 10


juv.dist2 <- juv.growth.dist * (1/100) * PreJuvWeight^(1-e1)

#juv.growth.dist
#juv.growth.step.dist = juv.growth.dist*(1/100)*PreJuvWeight^(1-e1)


plot(avec2[1:2000],JuvWeight[1:2000],main=paste("growth_rate",growth_step),ylim=c(0,5))
lines(avec2[1:2000],juv.growth.dist[1:2000],col="grey")
lines(avec2[1:2000],juv.dist2[1:2000],col="grey")

JuvWeight = ConvolveFunc(x1 = PreJuvWeight, y1 = juv.dist2, padsize = length(avec2)) #+ NewEggsT1*egg.dist
#JuvWeight = abs((JuvWeight))^(1/e1)
#JuvWeight = (JuvWeight)^(1/e1)
lines(avec2[1:2000],JuvWeight[1:2000],col="red")

#PreJuvWeight = (JuvWeight)^(e1)
PreJuvWeight = (JuvWeight)
juv.dist2 <- juv.growth.dist * (1/100) * PreJuvWeight^(1-e1)
lines(avec2[1:2000],juv.dist2[1:2000],col="grey")
JuvWeight = ConvolveFunc(x1 = PreJuvWeight, y1 = juv.dist2, padsize = length(avec2)) 
#JuvWeight = abs((JuvWeight))^(1/e1) + NewJuvT1*juv.growth.dist
JuvWeight = JuvWeight + NewJuvT1*juv.growth.dist
lines(avec2[1:2000],JuvWeight[1:2000],col="green")

#PreJuvWeight = (JuvWeight)^(e1)
PreJuvWeight = (JuvWeight)
juv.dist2 <- juv.growth.dist * (1/100) * PreJuvWeight^(1-e1)
JuvWeight = ConvolveFunc(x1 = PreJuvWeight, y1 = juv.dist2, padsize = length(avec2)) #+ NewEggsT1*egg.dist
#JuvWeight = abs((JuvWeight))^(1/e1)
#JuvWeight = (JuvWeight)^(1/e1) #Doing this produces the NaN values...
lines(avec2[1:2000],JuvWeight[1:2000],col="blue")


JuvLength <- ((JuvWeight/0.00001432)^(1/2.953)) #distribution of juv lengths (in mm) at this time step
  


####

i=500
Tp=temp[i]

emergence_weight = 35^2.953*0.00001432 #emergence mass (size 35 mm)
omega <- (a1*(Tp-b1)*(1-exp(c1*(Tp-d1)))) #growth rate
omega <- if_else(omega > 0, omega, 0) #growth rate can not be 'negative'

juv_weight <- rep(0,(date-j))

for(k in 1:(date - j + 1)) { #(k in j:date)
  if(k == 1) {
    juv_weight[k] <- (emergence_weight^(e1) + omega[k]*(1/100)*e1)^(1/e1) 
    
    #juv weight in grams
  }else{
    juv_weight[k] <- (juv_weight[k-1]^(e1) + omega[k]*(1/100)*e1)^(1/e1) 
    #juv weight in grams
  }
}
