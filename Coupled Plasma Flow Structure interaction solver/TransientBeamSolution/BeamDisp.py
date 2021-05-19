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

data_an = np.loadtxt( 'DisplacementAnalytical.dat' )
data_num = np.loadtxt( 'DisplacementNumerical.dat' )
data_err = np.loadtxt( 'SteadyError.dat' )

Data_size = np.shape(data_an)

# Display the sum
print('The data has {0} rows and {1} columns.'.format(Data_size[0], Data_size[1]))

N = data_err[:,0]
x = data_an[:,0]
w_an = data_an[:,1]*1000
w_num = data_num[:,1]*1000
w_err = data_err[:,1]*1000

font1 = {'family':'serif','color':'blue','size':15}
font2 = {'family':'serif','color':'darkred','size':15}

plt.figure(figsize=(8,2))
plt.plot(x, w_an)
plt.plot(x, w_an, "*", color='red')
ax = plt.gca() #you first need to get the axis handle
#ax.set_aspect(1) #sets the height to width ratio to 1.5. 
plt.xlim([-1, 1])
#plt.ylim([-0.15, 0])
plt.grid()
plt.title("Displacement at Steady-State", fontdict = font2)
plt.xlabel('x (m)', fontdict = font1)
plt.ylabel('w (mm)', fontdict = font1)
plt.legend(['Analytical Solution','Numerical Solution'])
#plt.clabel('Temperature')
plt.savefig('AnalyticalDisplacement.png')

plt.figure(figsize=(8,5))
plt.semilogy(N, w_err, "*")
ax = plt.gca() #you first need to get the axis handle
#ax.set_aspect(1) #sets the height to width ratio to 1.5. 
plt.grid(True, which="both")
plt.title("Displacement Error at Steady-State", fontdict = font2)
plt.xlabel('$N_{x}$', fontdict = font1)
plt.ylabel('$w_{error}$ (mm)', fontdict = font1)
plt.legend(['$L_{inf}$ Error'])
#plt.clabel('Temperature')
plt.savefig('DisplacementError.png')
