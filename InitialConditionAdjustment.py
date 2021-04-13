import numpy as np
import matplotlib.pyplot as plt
m = 11000        #microscopic area
s = 10000        #spread/projected area
e = 0         #internal pool of membrane
L = 7           #time lag before internal pool is recycled to membrane
D = 0.33
t = 120

#k1 = ((m/s))*0.025
k1 = 0.2       #endocytosis rate
k2 = 0.12       #exocytosis rate, remains constant
#k3 = 0.2
c = 0.8         #normalizing constant

#arrays of physical properties  above
M = [0]*t
S = [0]*t
E = [0]*t
M[0] = m
E[0] = e
S[0] = s

main = []            #array to store M/S ratio with time

####################################################################################################################

#lag phase, no recycling
def lagphase(lagtime):
    global k1
    for i in range(lagtime):
        M[i+1] = M[i] - k1*M[i]                        #membrane pool
        E[i+1] = E[i] + k1*M[i]                        #internal pool
        S[i+1] = S[i]                              
        k1 = ((M[i+1]/S[i+1])-1.1)*c
        
#recycling starts
def recycle(lagtime,rcltime):
    global k1
    for i in range(lagtime, rcltime-1):
        E[i+1] = E[i] + k1*M[i] - k2*E[i-L]           
        M[i+1] = M[i] - k1*M[i] + k2*E[i-L]
        S[i+1] = S[i]                                 
        k1 = (M[i+1]/S[i+1])*c
        #plt.scatter(i,k1)
        #plt.pause(0.1)

####################################################################################################################
        
#call functions
lagphase(L)
recycle(L,t)

####################################################################################################################
sub = []

#storing ratio data
for i in range(t):
    main.append((M[i]-S[i])/S[i])
    sub.append(M[i]/S[i])

print(k1)
    
#plt.plot(range(50),M[70:120], label = 'Microscopic area')
#plt.plot(range(50),S[70:120], label = 'Spread area')
plt.plot(range(t),M, label = 'Microscopic area')
plt.plot(range(t),S, label = 'Spread area')
#plt.plot(range(t),E, label = 'endocytosed area (ready to be recycled)')
plt.title("microscopic area and spread area change with time")
plt.xlabel("time(minutes)")
plt.ylabel("Surface areas")
plt.legend()
plt.show()
