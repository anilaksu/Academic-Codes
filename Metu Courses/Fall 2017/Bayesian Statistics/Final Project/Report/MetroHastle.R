####################################################
#	                                             #
#            Metropolis Hastling Algorithm for     #
#		      Internal Wave Field              #
#                  by Anil Aksu                    #    
#     						          #
####################################################	
library(mcmc)
library(MCMCpack)
library(MHadaptive)

## functions
source('getPerturbWave.R')
### A LINEAR REGRESSION EXAMPLE ####
## Define a Bayesian linear regression model
li_reg<-function(pars,data)
{
a<-pars[1] #intercept
#b<-pars[2] #slope
c<-pars[2] #second derivative 
sd_e<-pars[3] #error (residuals)
if(sd_e<=0){return(NaN)}
pred <- a +c*data[,1]^2
#+ b * data[,1]
log_likelihood<-sum( dnorm(data[,2],pred,sd_e, log=TRUE) )
prior<- prior_reg(pars)
return(log_likelihood + prior)
}
## Define the Prior distributions
prior_reg<-function(pars)
{
a<-pars[1] #intercept
#b<-pars[2] #slope
c<-pars[2] #second derivative
epsilon<-pars[3] #error
prior_a<-dnorm(a,0,10000,log=TRUE) ## non-informative (flat) priors on all
#prior_b<-dnorm(b,0,10000,log=TRUE) ## parameters.
prior_c<-dnorm(c,0,10000,log=TRUE) ## parameters.
prior_epsilon<-dgamma(epsilon,1,1/10000,log=TRUE)
return(prior_a + prior_c + prior_epsilon)
##Note: use more MCMC chains (i.e NC=10000) for more accurate results.
}
# simulate data
# the number of simulation data
N <- 100
x<-runif(N,0,20)
x_0 <- 10;
A <- 2;
A_0 <- 2;
sig_inc <- 5;
sig_top <- 2;
H <- 20;
# the perturbation envelope
y<-getPerturbWave(A,A_0,x_0,x,sig_inc,sig_top,N)
#y<-2.*x+2+rnorm(30,0,5)
d<-cbind(x,y)
mcmc_r<-Metro_Hastings(li_func=li_reg,pars=c(0,1,1),
par_names=c('a','c','epsilon'),data=d)
## For best results, run again with the previously
## adapted variance-covariance matrix.
mcmc_r<-Metro_Hastings(li_func=li_reg,pars=c(0,1,1),
prop_sigma=mcmc_r$prop_sigma,par_names=c('a','c','epsilon'),data=d)
mcmc_r<-mcmc_thin(mcmc_r)
plotMHdata(mcmc_r)
(mcmc_r)
