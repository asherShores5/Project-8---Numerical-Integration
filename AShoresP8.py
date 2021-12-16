#Asher Shores
#CST-305
#Dr. Ricardo Citro
#Project 8 - Numerical Integration
#19 December 2021
#This is my own work

#importing packages
from scipy import pi                #import pi from scipy
from scipy.integrate import quad    #import quad from scipy
from scipy.integrate import simps   #import quad from scipy
import numpy as np                  #import lib numpy as np for math
import matplotlib.pyplot as plt     #import lib matplotlib to graph

#defining the functions to integrate
a1 = lambda x: np.sin(x) + 1                 #eqn from 1a
b1 = lambda x: (3 * x) + (2 * (x ** 2))      #eqn from 1b
c1 = lambda x: np.log(x)                     #eqn from 1c1
c2 = lambda x: (x ** 2) - (x ** 3)           #eqn from 1c2

dd = lambda x: 0.0144*x + 10.66          #eqn from data download

n = 1000  #number of points used to graph

#function to plot riemanns using different vars / eqns
def riemannPlots(f, a, b, N):     #parameters are function, upper/lower limits, intervals
    
    x = np.linspace(a, b, N + 1)  #graphing space
    y = f(x)          #graphing space for specified function
    
    X = np.linspace(a, b, n * N + 1)  #space defined by n * N (recommended size)
    Y = f(X)          #graphing space for function of that size

    plt.figure(figsize=(15,5))  #dpi argument

    plt.subplot(1, 3, 1)        #plot 1 for LRAM
    plt.plot(X, Y, 'b')
    x_left = x[:-1]             #left endpoints
    y_left = y[:-1]

    plt.plot(x_left,y_left, 'b.', markersize = 10)    #ploting from endpoints
    plt.bar(x_left, y_left, width = (b - a) / N, alpha = 0.2, align = 'edge', edgecolor = 'b')
    #line formatting
    plt.title('Left Riemann Sum, N = {}'.format(N))
    #Title for plot

    #same thing repeated but for the second plot using MRAM
    #-------#
    plt.subplot(1, 3, 2)
    plt.plot(X, Y, 'b')
    x_mid = (x[:-1] + x[1:]) / 2
    y_mid = f(x_mid)

    plt.plot(x_mid, y_mid, 'b.', markersize = 10)
    plt.bar(x_mid,y_mid, width = (b - a) / N, alpha = 0.2, edgecolor = 'b')
    plt.title('Midpoint Riemann Sum, N = {}'.format(N))
    #-------#

    #same thing but for the third plot using RRAM
    #-------#
    plt.subplot(1, 3, 3)    #last plot  of 3
    plt.plot(X, Y, 'b')     #graph the slightly larger spaces
    x_right = x[1:]         #right endpoints 
    y_right = y[1:]         #right endpoints
    

    plt.plot(x_right, y_right, 'b.', markersize = 10)
    plt.bar(x_right, y_right, width = -(b - a) / N, alpha = 0.2, align = 'edge', edgecolor = 'b')
    plt.title('Right Riemann Sum, N = {}'.format(N))
    #-------#

    plt.show() #finally plot all three graphs together
    
#------------------------------------------------------------------------#

#Processing the visual data from above function to give summations
def riemannSum(f, a, b, N, m):   #Takes the function, upper/lower limits, intervals
    dx = (b - a)/N                                #number of steps to solve
    x_left = np.linspace(a, b - dx, N)            #step size LRAM
    x_midpoint = np.linspace(dx / 2, b - dx / 2, N)    #step size MRAM
    x_right = np.linspace(dx, b, N)                    #step size RRAM

    left_riemann_sum = np.sum(f(x_left) * dx)          #Calc. LRAM
    print("Left Riemann Sum:\t", left_riemann_sum*m)          #print LRAM

    midpoint_riemann_sum = np.sum(f(x_midpoint) * dx)       #Calc MRAM
    print("Midpoint Riemann Sum:\t", midpoint_riemann_sum*m)  #Print MRAM

    right_riemann_sum = np.sum(f(x_right) * dx)             #Calc RRAM
    print("Right Riemann Sum:\t", right_riemann_sum*m)        #Print RRAM
    
    
    
#------------------------------------------------------------------------#

#Calculate using summation function above,
#check against numpy's built in integral calculator
print("\nThe Riemman sum for 1a is:"); riemannSum(a1, -pi, pi, 4, 1)
print("The correct value for part 1a is: ", quad(a1, -pi, pi)[0])

print("\nThe Riemman sum for Part 1b is:"); riemannSum(b1, 0, 1, 1000,1)
print("The correct value for part 1b is: ", quad(b1, 0, 1)[0])

print("\nThe Riemman sum for Part 1c1 is:"); riemannSum(c1, 1, np.e, 1000,1)
print("The correct value for part 1c1 is: ", quad(c1, 1, np.e)[0])

print("\nThe Riemman sum for Part 1c2 is:"); riemannSum(c2, -1, 0, 1000,1)
print("The correct value for part 1c2 is: ", quad(c2, -1, 0)[0])


#-------------------------------------------------------------------------------------#

#Plot from very first function above using defined f(x) and ranges and subintervals
riemannPlots(a1, -pi, pi, 4)
riemannPlots(b1, 0, 1, 30)
riemannPlots(c1, 1, np.e, 30)
riemannPlots(c2, -1, 0, 30)


#-------------------------------------------------------------------------------------#
"PART 2"

#data from file download
xs = np.linspace(0, 30, 30)
ys = [10.2,9.4,9.7,10.5,11.4,11.2,12.3,12.4,11.1,10.5,9.7,12.4,9.4,10.3,10.6,10.1,11.5,11.3,10.6,11.7,12.2,11.8,10.2,11.1,10.2,11.5,10.6,10.9,11.4,10.3]

print("\n\nThe Riemann sum for data download is:")
riemannSum(dd, 1, 30, 30, 60)
print("\nThe above riemann sums are that of the linearization of the download data")
print("The data was plotted in excel and the \napproximate equation is y = 0.0144x + 10.66")

#Graphing that file data
plt.plot(xs, ys, 'b.', markersize = 1)               #plot the x and y space
plt.fill_between(xs, ys, alpha = 0.2, edgecolor = 'b')  #fills in data
plt.xlabel("Time (minutes)")                            #label x
plt.ylabel("Data Download Speed (Megabits per Second)") #label y
plt.title("Download Speed vs Time")
plt.show()

yys = ys
for i in range(len(yys)): yys[i] *= 60.0 #scale data to minnutes
print("\nThe integral of the data download speeds is:\t", simps(yys, xs), "megabits")
plt.show()     #show graph

riemannPlots(dd, 1, 30, 30)







