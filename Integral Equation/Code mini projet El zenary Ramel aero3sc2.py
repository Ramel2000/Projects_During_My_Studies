from scipy.integrate import *
from pylab import *
import numpy as np
from matplotlib import pyplot as plt
import math as m

## Équation intégrale
'''
## definition de la fonction Mat
def Mat(a,b,K,f):
    n = 100
    h = (b-a)/(n-1)
    t = np.linspace(a,b,n)
    x = np.linspace(a,b,n)
    F = np.zeros((n,1))
    for i in range (n):
        F[i,:] = f(t[i])
    A = np.zeros((n,n))
    for i in range (n):
        A[i,1] = K(x[i],t[0])+2*K(x[i],t[1])
        A[i,n-1] = K(x[i],t[n-1])
        for j in range (2,n-1):
            A[i,j] = 2*K(x[i],t[j])
    I=np.eye(n)
    M=(I-(h/2)*A)
    print(A)
    return M,F,t


##définition des parametre test
def f1(t):
    fct = cos((pi*t)/2)-(8/pi)
    return(fct)
def u1(t):
    fct = cos((pi*t)/2)
    return(fct)
def K1(x,t):
    K = 2
    return(K)

    
M,F,t = Mat(-1,1,K1,f1)
Mi = np.linalg.inv(M)
#définition de l'approximation
V=Mi@F
V = solve(M,F)
#définition de la solution exacte
U = np.zeros((F.size,1))
for i in range (F.size):
    U[i,:] = u1(t[i])
#représentation graphique    
plt.plot(t,V,color="r",label="Solution approchée")
plt.plot(t,U,color="b",label="Solution exacte")
plt.title(label ="Comparaison de la solution approchée à la solution exacte (n=100)")
plt.legend()
plt.grid()
plt.show()
#représentation de l'erreur
print(("l'erreur de l'approximation est de",np.linalg.norm(U-V)))


##définition des parametre pour la fonction de love
def f2(t):
    fct = 1
    return(fct)
def K2(x,t):
    K = (1/pi)*1/(1+(x-t)**2)
    return(K)

M,F,t = Mat(-1,1,K2,f2)
Mi = np.linalg.inv(M)
#définition de l'approximation
U=Mi@F
#représentation graphique
plt.plot(t,U,color="r",label="Courbe de l'équation de love")
plt.title("Équation de Love en électrostatique")
plt.xlabel("Temps (s)")
plt.ylabel("Amplitude")
plt.legend()
plt.grid()
plt.show()

#Circuit RLC
def rlcprim (Y,t):
    yprime=np.zeros(2)
    yprime[0]=Y[1]/0.000001
    yprime[1]=(10-Y[0]-3*Y[1])/0.5
    return yprime


y0=np.array([0,0])
t=np.linspace(0,2,200)
yprime=odeint(rlcprim,y0,t)
plt.plot(t,yprime[:,0],label="Tension")
plt.title("Représentation de l'évolution de la tension au cours du temps")
plt.grid()
plt.xlabel("Temps (s)")
plt.ylabel("Amplitude")
plt.legend()
plt.show()
plt.plot(t,yprime[:,1],label="Courant")
plt.title("Représentation de l'évolution du courant au cours du temps")
plt.grid()
plt.xlabel("Temps (s)")
plt.ylabel("Amplitude")
plt.legend()
plt.show()





def moteurCC (Y,t):
    u=0
    if (t>=10 and t<=50):
        u=5
    else:
        u=0
    yprime=np.zeros(2)
    yprime[0]=(u-5*Y[0]-0.2*Y[1])/(50*(10**(-3)))
    yprime[1]=(0.1*Y[0]-0.01*Y[1])/0.05
    return yprime


y0=np.array([0,0])
t=np.linspace(0,80,100)
yprime=odeint(moteurCC,y0,t)
plt.plot(t,0.1*yprime[:,0],label = "couple moteur")
plt.title("Représentation du couple moteur au cours du temps")
plt.xlabel("Temps (s)")
plt.ylabel("Couple moteur (N/m)")
plt.grid()
plt.legend()
plt.show()
plt.plot(t,yprime[:,1],label = "vitesse angulaire")
plt.title("Représentation de la vitesse angulaire au cours du temps")
plt.xlabel("Temps (s)")
plt.ylabel("Vitesse angulaire (rad/s)")
plt.legend()
plt.grid()
plt.show()


def fusee(Y, t) :
    D = 4
    a = 8*(10**3)
    g = 9.81
    k = 0.1
    u = 2*(10**3)
    yprime = np.zeros(3)
    if (Y[1]<=80):
        Y[1]=80
        D=0
        
    yprime[0] = ((D*u)/Y[1])-g-k*m.exp(-(Y[2]/a))*(Y[0]**2)/Y[1]
    yprime[1] = -D
    yprime[2] = Y[0]
    return yprime




y0=np.array([0,400,0])
t=np.linspace(0,160,100)
yprime=odeint(fusee,y0,t)
plt.title("Représentation de la vitesse au cours du temps")
plt.plot(t,yprime[:,0],label="courbe de la vitesse")
plt.xlabel("Temps (s)")
plt.ylabel("Vitesse (m/s)")
plt.grid()
plt.legend()
plt.show()
t=np.linspace(0,80,100)
plt.title("Représentation de la trajectoire au cours du temps")
plt.plot(t,yprime[:,2],label="courbe de la trajectoire")
plt.xlabel("Temps (s)")
plt.ylabel("Altitude (m)")
plt.grid()
plt.legend()
plt.show()
'''

t=np.linspace(0,10,480)
h=10/480

def Proie_predateur(Y, t):
    yprime = np.zeros(2)
    yprime[0] = 3*Y[0]-1*Y[0]*Y[1]
    yprime[1] = -2*Y[1]+1*Y[0]*Y[1]
    return yprime


def Euler(f,t,h,Y0):
    N=t.size
    Ye=np.zeros((N,2))
    Ye[0,:]=Y0
    for k in range (0,N-1):
        Ye[k+1,:]= Ye[k,:]+h*f(Ye[k,:],t[k])
    return Ye
#1
y0=np.array([5,0])
Y=Euler(Proie_predateur,t,h,y0)
plt.plot(t,Y[:,0],label="Proie")
plt.xlabel("Temps (Année)")
plt.ylabel("Population")
plt.title("Evolution des proies en l'absence de prédateur")
plt.grid()
plt.legend()
plt.show()
#2
y0=np.array([0,3])
Y=Euler(Proie_predateur,t,h,y0)
plt.plot(t,Y[:,1],label="Prédateur")
plt.xlabel("Temps (Année)")
plt.ylabel("Population")
plt.title("Evolution des prédateurs en l'absence de proie")
plt.grid()
plt.legend()
plt.show()
#3
y0=np.array([5,3])
Y=Euler(Proie_predateur,t,h,y0)
plt.plot(t,Y[:,0],label="Proie")
plt.plot(t,Y[:,1],label="Prédateur")
plt.xlabel("Temps (Année)")
plt.ylabel("Population")
plt.title("Evolution des prédateurs et proies grâce à la méthode d'Euler")
plt.grid()
plt.legend()
plt.show()

#4
yprime=odeint(Proie_predateur,y0,t)


plt.plot(t,yprime[:,0],label="Proie")
plt.plot(t,yprime[:,1],label="Prédateur")
plt.xlabel("Temps (Année)")
plt.ylabel("Population")
plt.title("Evolution des prédateurs et proies grâce au solveur Odeint")
plt.grid()
plt.legend()
plt.show()
#5
plt.title("Représentation du portait de phase")
plt.plot(Y[:,0],Y[:,1],label="Portait par la méthode d'Euler")
plt.plot(yprime[:,0],yprime[:,1],label="Portait par le solveur Odeint")
plt.grid()
plt.legend()
plt.show()

#6
def Proie_predateur(Y, t):
    yprime = np.zeros(2)
    yprime[0] = 3*Y[0]-1*Y[0]*Y[1]
    yprime[1] = -2*Y[1]+1*Y[0]*Y[1]
    return yprime
y0=np.array([30,3])

yprime=odeint(Proie_predateur,y0,t)
print(yprime)
plt.plot(t,yprime[:,0],label="Proie")
plt.plot(t,yprime[:,1],label="Prédateur")
plt.xlabel("Temps (Année)")
plt.ylabel("Population")
plt.title("Evolution des prédateurs et proies grâce au solveur Odeint")
plt.grid()
plt.legend()
plt.show()

