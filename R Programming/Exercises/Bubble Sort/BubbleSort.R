####################################################
#	                                               #
#	 Bubble Sort Algorithm by Anil Aksu            #
#   It is developed to show some basics of R       #
#     						         			   #
####################################################	

## functions
source('getSorted.R')
 
## numbers array 
Numbers <- integer(5)

## this function gets numbers from console
for (i in 1:5){
	# this converts data into integer as.integer()
	Numbers[i] <- as.integer(readline("Please enter a number"))
}


## let's output them 
print("The numbers")
print(Numbers)
## let's sort them out

SortedNumbers <- getSorted(Numbers);

## let's output them 
 print("The  sorted numbers")
 print(SortedNumbers)