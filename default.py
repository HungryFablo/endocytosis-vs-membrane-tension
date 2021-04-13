import numpy as np
import matplotlib.pyplot as plt
n = 60000        #microscopic area
s = 6000        #spread/projected area
e = 10000         #internal pool of membrane
m = n-e
L = 7           #time lag before internal pool is recycled to membrane
D = 0.33
t = 120

k1 = ((m/s)-5)*0.025
#k1 = 0.38       #endocytosis rate
k2 = 0.18       #exocytosis rate, remains constant
k3 = 0.2
c = 0.025        #normalizing constant

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
        k1 = ((M[i+1]/S[i+1])-5)*c
        
#recycling starts
def recycle(lagtime,rcltime):
    global k1
    for i in range(lagtime, rcltime):
        E[i+1] = E[i] + k1*M[i] - k2*E[i-L]           
        M[i+1] = M[i] - k1*M[i] + k2*E[i-L]
        S[i+1] = S[i]                                 
        k1 = ((M[i+1]/S[i+1])-5)*c
        #plt.scatter(i,k1)
        #plt.pause(0.1)
#deadhesion starts
def deadhesion(time1, time2):
    global k1
    for i in range(time1, time2):
        E[i+1] = E[i] + k1*M[i] - k2*E[i-L]           
        M[i+1] = M[i] - k1*M[i] + k2*E[i-L]         
        S[i+1] = S[i] - k3*S[i]
        k1 = ((M[i+1]/S[i+1])-5)*c
        #plt.scatter(i,k1)
        #plt.pause(0.1)
def steadyphase(time1,time):
    global k1
    for i in range(time1, time-1):
        E[i+1] = E[i] + k1*M[i] - k2*E[i-L]           
        M[i+1] = M[i] - k1*M[i] + k2*E[i-L]
        S[i+1] = S[i]
        k1 = ((M[i+1]/S[i+1])-5)*c
        #plt.scatter(i,k1)
        #plt.pause(0.1)
####################################################################################################################
        
#call functions
lagphase(L)
recycle(L,75)
deadhesion(75,85)
steadyphase(85,t)

####################################################################################################################
sub = []

#storing ratio data
for i in range(t):
    main.append((M[i]-S[i])/S[i])
    sub.append(M[i]/S[i])

    
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
plt.plot(range(50),main[70:120])
plt.title("del(A)/A ratio vs time, Lag= "+str(L)+"mins, k1= 0.38, k2= "+str(k2)+", k3= "+str(k3))
#plt.ylim(4,10)
plt.xlabel("time(minutes)")
plt.ylabel("del(A)/A")
#plt.plot(range(50),sub[70:120])
plt.show()

#for i in range(t):
#    print(main[i+1]-main[i])
