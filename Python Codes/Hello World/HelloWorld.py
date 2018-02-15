# Hello World Python

# libraries
import math

# functions
def sinus(x,y):
	# global f
 	return math.sin(x+y)
	# return f

def happyBirthdayRabia():
    print("Happy Birthday Rabia!")


def main():
    happyBirthdayRabia()
	
# the main script	
#main()
# x and y arrays
x=[]
y=[]
arr=[]
# input for getFunction
for i in range(1,11):
	#print i
	x.append(i)
	y.append(2.*i)
	print x[i-1],y[i-1],
	# the evaluation of getFunction
	print sinus(x[i-1],y[i-1])
	arr.append(sinus(x[i-1],y[i-1]))
	print arr[i-1]
	
print "merhaba anil"