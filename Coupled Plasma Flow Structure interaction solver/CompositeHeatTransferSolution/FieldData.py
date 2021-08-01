# Libraries
import sys
sys.path.append("NumPy_path")
import numpy as np
import matplotlib.pyplot as plt

# This program adds two numbers

num1 = 1.5
num2 = 6.3

# Add two numbers
sum = num1 + num2

# Display the sum
print('The sum of {0} and {1} is {2}'.format(num1, num2, sum))

data = np.loadtxt( 'Temperature2D.dat' )

Data_size = np.shape(data)

# Display the sum
print('The data has {0} rows and {1} columns.'.format(Data_size[0], Data_size[1]))


x = data[:,0].reshape( (10, 10))
y = data[:,1].reshape( (10, 10))
Temperature = data[:,2].reshape( (10, 10))
#y = np.eeshape(data[:][1], (8, 4))
#Temperature = np.arange(Data_size[0]).reshape(data[:][2], (8, 4))
#Temperature = data[:][2]    # Temperature array 


plt.figure(1)
plt.contourf(x, y, Temperature, 40, cmap='jet')
plt.colorbar(label='Temperature (K)');
ax = plt.gca() #you first need to get the axis handle
ax.set_aspect(2.5) #sets the height to width ratio to 1.5. 
plt.xlabel('x (m)')
plt.ylabel('y (m)')
#plt.clabel('Temperature')
plt.savefig('TemperatureDistribution.png')