getMean<- function(Numbers){
	
####################################################
#	                                             #
#  This function calculates the mean of given	   #
#  random numbers                                  #
#	                                             #
#	Input:  Random Numbers	    			   #
#	                                             #
#	Output: Mean of random numbers	         #
#	                                             #
####################################################

## the length of Number Array
numSize <- length(Numbers[,1]);
## Mean of these numbers
Mean <- sum(Numbers)/numSize;
 return(Mean)	
}