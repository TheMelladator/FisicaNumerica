import pandas as pd
import numpy as np
import matplotlib.pyplot as mp

def g(x,a):
    deno=a[2]+(x-a[1])**2
    return a[0]/deno

def f1(f,g,E,a,s):
    return (f-g)/((a[2] + (E-a[1])**2)*s**2)

def f2(f,g,E,a,s):
    return ((f-g)*(E-a[1]))/(((a[2] + (E-a[1])**2)**2)*s**2)

def f3(f,g,E,a,s):
    return (f-g)/(((a[2] + (E-a[1])**2)**2)*s**2)

def F1(f,g,E,a,s):
    A=0
    for i in range(len(E)):
        A += f1(f[i],g(E[i],a),E[i],a,s[i])
    return A

def F2(f,g,E,a,s):
    A=0
    for i in range(len(E)):
        A += f2(f[i],g(E[i],a),E[i],a,s[i])
    return A

def F3(f,g,E,a,s):
    A=0
    for i in range(len(E)):
        A += f3(f[i],g(E[i],a),E[i],a,s[i])
    return A

data = pd.read_csv('ResonanciaBW.txt',header=0,delim_whitespace=True)
E=data.iloc[:,0]
f=data.iloc[:,1]
sigma=data.iloc[:,2]

E = np.float64(E)
f = np.float64(f)
s = np.float64(sigma)

a=[1000,76,900]
b=np.zeros(3)
for i in range(3):
    b[i]=a[i]

MaxN=100
err=1e-12
FP=np.zeros((3,3))
dx=0.1

for k in range(MaxN):
    F = [F1(f,g,E,a,s),F2(f,g,E,a,s),F3(f,g,E,a,s)]
    if F[0]<=err and F[1]<=err and F[2]<=err:
        n=k
        break
    for i in range(3):
        for j in range(3):
            if i == 0:
                b[j]=b[j]+dx
                FP[i][j]=(F1(f,g,E,b,s)-F1(f,g,E,a,s))/dx
                b[j]=b[j]-dx
            elif i==1:
                b[j]=b[j]+dx
                FP[i][j]=(F2(f,g,E,b,s)-F2(f,g,E,a,s))/dx
                b[j]=b[j]-dx
            else:
                b[j]=b[j]+dx
                FP[i][j]=(F3(f,g,E,b,s)-F3(f,g,E,a,s))/dx
                b[j]=b[j]-dx
    DA=np.dot(np.linalg.inv(FP),F)
    DA*=-1
    for i in range(3):
        a[i]=a[i]+DA[i]
        b[i]=a[i]
print("fr= ",a[0], " Er= ",a[1]," Gamma= ",np.sqrt(4*a[2]))

x=np.linspace(0,200)
S=a[0]/(a[2]+(x-a[1])**2)
mp.plot(E,f,'ro',label="Datos experimentales")
mp.plot(x,S,'b',label="Ajuste por minimos cuadrados")
mp.xlabel("Energia (MeV)")
mp.ylabel("f(E) (MeV)")
mp.legend(loc="upper right")
mp.show()
