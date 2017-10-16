####################################################
#	                                             #
#   Mean and Standard Deviation Calculation        #
#                 by Anil Aksu                     #
#   It is developed to show some basics of R       #
#     						         #
####################################################	

## library required to read excel files
require(gdata)
## functions
source('getMean.R')
source('getSigma.R')
## the random data read from excel file
RandomNumbers <- read.xls("RandomNumbers.xlsx", perl = "C:\\Perl\\bin\\perl.exe")
## let's calculate the mean 
Mean <- getMean(RandomNumbers)
## the standard deviation 
Sigma <- getSigma(RandomNumbers)
## this function gets numbers from console

## let's output them 
print("Random Numbers")
print(RandomNumbers[,1])
print("The Mean")
print(Mean)
print("The Standard Deviation")
print(Sigma)

