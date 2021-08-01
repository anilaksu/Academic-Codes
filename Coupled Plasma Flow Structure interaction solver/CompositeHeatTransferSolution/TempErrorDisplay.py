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

data_err = np.loadtxt( 'TempRelErr.dat' )

Data_size = np.shape(data_err)

# Display the sum
print(Data_size)
print('The data has {0} rows and {1} columns.'.format(Data_size[0], Data_size[1]))

N = data_err[:,0]
Err_5 = data_err[:,1]
Err_10 = data_err[:,2]
Err_20 = data_err[:,3]
Err_30 = data_err[:,4]
Err_40 = data_err[:,5]


font1 = {'family':'serif','color':'blue','size':15}
font2 = {'family':'serif','color':'darkred','size':15}

plt.figure(figsize=(8,5))
plt.semilogy(N, Err_5 , "*", color='red')
plt.semilogy(N, Err_10 , "*", color='blue')
plt.semilogy(N, Err_20 , "*", color='green')
#plt.semilogy(N, Err_30 , "*", color='purple')
#plt.semilogy(N, Err_40 , "*", color='pink')
ax = plt.gca() #you first need to get the axis handle
#ax.set_aspect(1) #sets the height to width ratio to 1.5. 
plt.grid(True, which="both")
plt.title("Temperature Error at Steady-State", fontdict = font2)
plt.xlabel('$N_{x}$', fontdict = font1)
plt.ylabel('$L_{inf}(T_{err})$', fontdict = font1)
plt.legend(['$N_{y} = 5$','$N_{y} = 10$','$N_{y} = 20$'])
plt.savefig('TemperatureError.png')
