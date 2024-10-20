from scipy.integrate import *
from pylab import *
import numpy as np
from matplotlib import pyplot as plt
tfin=4
t=np.linspace(0,tfin,100)


N=t.size
h=tfin/100
theta=(pi/2)*cos(sqrt(9.81)*t)
thetap=-(pi/2)*sqrt(9.81)*sin(sqrt(9.81)*t)
#Generation de la courbe
plt.plot(t,theta,label='θ(t)')
plt.title('Angle θ en fonction du temps')
plt.xlabel('Temps (s)')
plt.ylabel('Amplitude')
plt.legend()
plt.grid()
plt.show()
Y0=np.array([pi/2,0])

def pendule(y,t):
    yprime=np.zeros(2)
    yprime[0]=y[1]
    yprime[1]=-9.81*y[0]
    return yprime
pendule(Y0,1)



#Euler explicit#
def Euler(f,t,h,Y0):
    N=t.size
    Ye=np.zeros((N,2))
    Ye[0,:]=Y0
    for k in range (N-1):
        Ye[k+1,:]= Ye[k,:]+h*f(Ye[k,:],t[k])
    return t,Ye

t1,Ye=Euler(pendule,t,h,Y0)
#Generation de la courbe
plt.plot(t1,Ye[:,0],label='Euler explicit')
plt.plot(t,theta,label='θ(t)')
plt.title('Angle θ en fonction du temps calculé\n grâce à la methode de Euler explicit')
plt.xlabel('Temps (s)')
plt.ylabel('Amplitude')
plt.grid()
plt.legend()
plt.show()

#Runge Kutta
def RungeKutta(f,t,h,Y0):
    N=t.size
    Yrk = zeros((N,2))
    Yrk[0,:] = Y0
    for k in range(N-1):
        Y = Yrk[k,:]
        k1 = f(Y,t[k])
        k2 = f(Y+h*k1/2,t[k]+h/2)
        k3 = f(Y+h*k2/2,t[k]+h/2)
        k4 = f(Y+h*k3,t[k]+h)
        Yrk[k+1,:] = Y +(h/6)*(k1+2*k2+2*k3+k4)
    return t,Yrk
t2,Yrk=RungeKutta(pendule,t,h,Y0)
#Generation de la courbe
plt.plot(t1,Ye[:,0],label='Euler explicit')
plt.plot(t,theta,label='θ(t)')
plt.plot(t2,Yrk[:,0],label='Runge Kutta')
plt.title('Angle θ en fonction du temps calculé\n grâce à la methode de Runge Kutta')
plt.xlabel('Temps (s)')
plt.ylabel('Amplitude')
plt.grid()
plt.legend()
plt.show()

#Solveur Odeint
Yode = odeint(pendule, Y0,t)
#Generation de la courbe
plt.plot(t1,Ye[:,0],label='Euler explicit')
plt.plot(t,theta,label='θ(t)')
plt.plot(t2,Yrk[:,0],label='Runge Kutta')
plt.plot(t,Yode[:,0],label='odeint')
plt.title('Angle θ en fonction du temps calculé\n grâce au solveur Odeint')
plt.xlabel('Temps (s)')
plt.ylabel('Amplitude')
plt.legend()
plt.grid()
plt.show()
#bonus
#1er partie
#Generation de la courbe
plt.title('Portait de phase des courbes')
plt.plot(Ye[:,0],Ye[:,1],label='Euler explicit')
plt.plot(Yrk[:,0],Yrk[:,1],label='Runge Kutta')
plt.plot(Yode[:,0],Yode[:,1],label='odeint')
plt.plot(theta,thetap,label='Méthode Analytique')
plt.xlabel("Y")
plt.ylabel("Y'")
plt.legend()
plt.grid()
plt.show()
#bonus
#2eme partie
def RungeKutta2(f,t,h,Y0):
    N=t.size
    Yrk = zeros((N,2))
    Yrk[0,:] = Y0
    for k in range(N-1):
        Y = Yrk[k,:]
        k1 = f(Y,t[k])
        k2 = f(Y+h*k1/2,t[k]+h/2)
        Yrk[k+1,:] = Y +(h*k2)
    return t,Yrk
tk2,Yrk2 = RungeKutta2(pendule,t,h,Y0)
#Generation de la courbe
plt.plot(tk2,Yrk2[:,0],label="Runge Kutta d'ordre 2")
plt.plot(t1,Ye[:,0],label='Euler explicit')
plt.plot(t,theta,label='θ(t)')
plt.plot(t2,Yrk[:,0],label='Runge Kutta')
plt.plot(t,Yode[:,0],label='odeint')
plt.title("Angle θ en fonction du temps calculé\n grâce à la methode de Rungé Kunta d'ordre 2")
plt.xlabel('Temps (s)')
plt.ylabel('Amplitude')
plt.legend()
plt.grid()
plt.show()
#Partie suspension 
def suspension (Y,t):
    yprime=zeros(4)
    yprime[0]=Y[2] 
    yprime[1]=Y[3]
    yprime[2]=(-1200*Y[2] + 1200*Y[3] -(55000*Y[0])+ 5000*Y[1])/15
    yprime[3]=(1200*Y[2]-1200*Y[3]+ 5000*Y[0] -5000*Y[1]-1000)/200
    return yprime

y0=array([0,0,0,0])
t=linspace(0,3,100)
yprime=odeint(suspension,y0,t)
#Generation de la courbe
plt.plot(t,yprime[:,0],label='x1(t)/(roue)')
plt.plot(t,yprime[:,1],label='x2(t)/(caisse)')
plt.title("Representation de l'affaissement de la roue et de la caisse\n au cours du temps")
plt.xlabel('Temps (s)')
plt.ylabel('Amplitude')
plt.grid()
plt.legend()
plt.show()















           
