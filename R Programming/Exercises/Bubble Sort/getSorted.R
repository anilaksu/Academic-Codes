getSorted<- function(Numbers){
	
####################################################
#	                                             #
#    This function is developed to find membership #
#    values for Correction Factor                  #
#	                                             #
####################################################	 
	
####################################################
#	                                             #
#	Input:   Numbers in arbitrary order		   #
#	                                             #
####################################################		
	
####################################################
#	                                             #
#	Output:   Numbers in decreasing order        #
#	                                             #
####################################################

## the length of Number Array
numSize <- length(Numbers)
## let's define sorted array of Numbers 
SortedNumbers <- integer(numSize)

SortedNumbers <- Numbers;

## let's start sorting
for(i in 1:(numSize-1)){
	for (j in i:numSize){
 		if(SortedNumbers[i] < SortedNumbers[j]) {
			# dummy variable
			dummy1 <- SortedNumbers[i]
 	 		SortedNumbers[i] <- SortedNumbers[j];
			SortedNumbers[j] <- dummy1;
 		}
 	 }
}
 return(SortedNumbers)	
}
