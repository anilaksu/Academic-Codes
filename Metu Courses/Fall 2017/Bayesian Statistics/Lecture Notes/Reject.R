####################################################
#	                                             #
#            Rejection sampling plots              #
#                  by Anil Aksu                    #    
#     						         #
####################################################	


require(SMPracticals)
## rejection sampling

## the range of sampling
x=seq(-10,10,length=101)
## this function gets numbers from console
posterior=0.6*dnorm(x, mean = 0, sd = 4, log = FALSE)+0.4*dnorm(x, mean = 6, sd = 2, log = FALSE)
envelope=2*dnorm(x, mean = 2, sd = 5, log = FALSE)
## let's plot them
plot(range(x), range(c(posterior,envelope)), type='n', xlab=expression(paste(theta)), ylab="")
lines(x, posterior, type='l', col='blue',lwd=5)
lines(x, envelope, type='l', col='red',lwd=5)

title("Rejection Sampling Plot")
legend("topright", legend=c(expression(paste("L(", theta, ")",pi,"(", theta, ")")),expression(paste(mu,"g(", theta, ")"))), 
   lty=1, col=c('blue', 'red'),inset = .02)

## British Coal Mining accidents
data(coal)
# years of coal mining accidents
years <- unique(as.integer(coal$date))
# the number of accidents in each year
accident <- integer(length(years))
for (i in 1:length(years)){
accident[i]<-sum(as.integer(coal$date) == years[i])
}

plot(years ,accident, col='blue',lwd=2, xlab="year", ylab="# of disasters")
#rug(coal$date)

