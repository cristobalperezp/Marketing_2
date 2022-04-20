#--------------------------------------------------------------------------------------------------
# Author		:	Marcel Goic
# Usage			:	source("discrete duration models - v01.r", print.eval=TRUE)
# Notes     : For cycles can be replaced by matrix multiplication to improve efficiency 
# Description	:	Estimate Simple Regression Model
#						- v.0.0 (01/10/2018). shifted geometric and shifted geometric with two segments
#						- v.0.1 (02/10/2018). add the beta-geometric 
#--------------------------------------------------------------------------------------------------

rm(list=ls())				  # clear the list of objects
graphics.off()				# clear the list of graphs
options(digits = 5)		# number of digits to display

# read the data
mydata <- read.table("retentiondata.txt", sep="\t", header=TRUE)
# compute how many customers are lost every year
mydata$lost = c(NA, -diff(mydata$ncust))                      

calibration = 2:8;  nc = length(calibration)
validation  = 9:12; nv = length(validation)

caldata = mydata[calibration, ]
valdata = mydata[validation, ]

#loglikelihod for the shifted geometric model --------------------------------
ll.sg <-function(theta){
  # transform variable to impose is in the (0,1) range
  mytheta = exp(theta)/(1+exp(theta))
  # define auxiliary arrays
  pr = array(NA, dim=nc+1)
  ll = array(NA, dim=nc+1)
  #detection likelihood
  for(i in 1:nc){
    pr[i] = mytheta*(1-mytheta)^(caldata$year[i]-1)
    ll[i] = caldata$lost[i]*log(pr[i])
  }
  #survival likelihood
  pr[nc+1] = 1-sum(pr[1:nc])
  ll[nc+1] = caldata$ncust[nc]*log(pr[nc+1])
  return(-sum(ll))
}

theta.start = 0.666
cat("Maximizing likelihood shifted geometric ... \n")
mle.sg = optim(par=theta.start, fn=ll.sg, hessian=TRUE, method="BFGS", control = list(maxit=30000, trace=TRUE, REPORT=10))
mytheta = exp(mle.sg$par)/(1+exp(mle.sg$par))
print(mytheta)

#loglikelihod for the shifted geometric model with two segments --------------
ll.sg2s <-function(theta){
  # transform variable to impose is in the (0,1) range
  mytheta = exp(theta)/(1+exp(theta))
      # mytheta[1]: churn probability segment 1
      # mytheta[2]: churn probability segment 2
      # mytheta[3]: probability belonging segment 1
  # define auxiliary arrays
  pr = array(NA, dim=c(2,nc+1))
  ll = array(NA, dim=nc+1)
  #detection likelihood
  for(i in 1:nc){
    pr[1,i] = mytheta[1]*(1-mytheta[1])^(caldata$year[i]-1)
    pr[2,i] = mytheta[2]*(1-mytheta[2])^(caldata$year[i]-1)
    ll[i] = caldata$lost[i]*log(mytheta[3]*pr[1,i] + (1-mytheta[3])*pr[2,i])
  }
  #survival likelihood
  pr[1,nc+1] = 1-sum(pr[1,1:nc])
  pr[2,nc+1] = 1-sum(pr[2,1:nc])
  ll[nc+1] = caldata$ncust[nc]*log(mytheta[3]*pr[1,nc+1]+(1-mytheta[3])*pr[2,nc+1])
  return(-sum(ll))
}

theta.start = rep(0.1,3)
cat("Maximizing likelihood shifted geometric with two segments ...\n")
mle.sg2s = optim(par=theta.start, fn=ll.sg2s, hessian=TRUE, method="BFGS", control = list(maxit=30000, trace=TRUE, REPORT=10))
mytheta <- exp(mle.sg2s$par)/(1+exp(mle.sg2s$par))
print(mytheta)

#loglikelihod for beta geometric model ------------------------------------------------------------
ll.bg <-function(theta){
  # transform variable to impose they are positive
  myalpha = exp(theta[1])
  mybeta  = exp(theta[2])
  # define auxiliary arrays
  pr = array(NA, dim=nc+1)
  ll = array(NA, dim=nc+1)
  #detection likelihood
  for(i in 1:nc){
    pr[i] = beta(myalpha+1,mybeta+caldata$year[i]-1) / beta(myalpha,mybeta)
    ll[i] = caldata$lost[i]*log(pr[i])
  }
  #survival likelihood
  pr[nc+1] = 1-sum(pr[1:nc])
  ll[nc+1] = caldata$ncust[nc]*log(pr[nc+1])
  return(-sum(ll))
}

theta.start = rep(0.666,2)
cat("Maximizing likelihood shifted beta geometric ... \n")
mle.bg = optim(par=theta.start, fn=ll.bg, hessian=TRUE, method="BFGS", control = list(maxit=30000, trace=TRUE, REPORT=10))
myalpha = exp(mle.bg$par[1])
mybeta  = exp(mle.bg$par[2])
print(c(myalpha, mybeta))