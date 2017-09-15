getSigma<- function(Numbers){
	
####################################################
#	                                             #
#  This function calculates the mean of given	   #
#  random numbers                                  #
#	                                             #
#	Input:  Random Numbers	    			   #
#	                                             #
#	Output: standard deviation of random numbers #
#	                                             #
####################################################

## the length of Number Array
numSize <- length(Numbers[,1]);
## Mean of these numbers
Mean <- getMean(Numbers);

## standard deviation
Sigma=0.
for (i in 1:numSize){
	# this converts data into integer as.integer()
	Sigma <- Sigma+Numbers[i,1]^2.
}
Sigma <- sqrt((Sigma/(numSize-1)));
 return(Sigma)	
}