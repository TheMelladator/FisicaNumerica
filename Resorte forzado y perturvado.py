import numpy as np
import matplotlib.pylab as mplot


x0=0.03
m=0.5
mu=0.25
g=9.81
k=(mu*g*m)/x0
n=500

deltaT=15/n
t=np.zeros(n)
i=0

while(i<n):
    t[i]=i*deltaT
    i=i+1

def PosiSinPertu(t):
     return -1*m*g*np.cos(np.sqrt(k/m)*t)+mu*g

def PosiConPertu(t):
     return -1*(mu*g+3/(26*m))*np.cos(np.sqrt(k/m)*t)+mu*g+3/(26*m)*np.cos(5*t)


mplot.plot(t,PosiSinPertu(t),label="Sin perturbacion")
mplot.plot(t,PosiConPertu(t),label="Con perturbación")
mplot.legend()
mplot.show()