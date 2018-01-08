getPerturbWave<-function(pars,data)
{
a<-pars[1] #intercept
b<-pars[2] #slope
sd_e<-pars[3] #error (residuals)
if(sd_e<=0){return(NaN)}
pred <- a + b * data[,1]
log_likelihood<-sum( dnorm(data[,2],pred,sd_e, log=TRUE) )
prior<- prior_reg(pars)
return(log_likelihood + prior)
}