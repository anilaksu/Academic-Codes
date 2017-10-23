####################################################
#	                                             #
#      A Sample Triplot by Anil Aksu               #
#   It is developed to show some basics of R       #
#     						         #
####################################################	

## the range of sampling
x=seq(-4,4,length=101)
## this function gets numbers from console
prior=dnorm(x, mean = 0, sd = 1.5, log = FALSE)
likelihood=dnorm(x, mean = 1.5, sd = 0.7, log = FALSE)
posterior=dnorm(x, mean = 1, sd = 0.5, log = FALSE)


## let's plot them
plot(range(x), range(c(likelihood,prior,posterior)), type='n', xlab="x", ylab="f(x)")
lines(x, prior, type='l', col='blue')
lines(x, likelihood, type='l', col='red')
lines(x, posterior, type='l', col='green')

title("Prior, Likelihood and Posterior Distribution")
legend(
  "topright", 
  lty=c(1,1,1), 
  col=c("blue", "red", "green"), 
  legend = c("prior", "likelihood","posterior")
)