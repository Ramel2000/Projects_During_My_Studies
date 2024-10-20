import numpy as np
import matplotlib.pyplot as plt
import time
from math import *
from pprint import pprint
from typing import Union

#Décomposition QR avec houseolder
def householder_vectorized(a):
    v = a / (a[0] + np.copysign(np.linalg.norm(a), a[0]))
    v[0] = 1
    tau = 2 / (v.T @ v)
    return v,tau
def qr_decomposition(A: np.ndarray) -> Union[np.ndarray, np.ndarray]:
    m,n = A.shape
    R = A.copy()
    Q = np.identity(m)
    for j in range(0, n):
        v, tau = householder_vectorized(R[j:, j, np.newaxis])
        H = np.identity(m)
        H[j:, j:] -= tau * (v @ v.T)
        R = H @ R
        Q = H @ Q
    Q = Q[:n].T
    R = np.triu(R[:n])
    Qt = np.transpose(Q)
    return Q,R,Qt

## Alternative-au-cocatenate

def concate (a,b):
    li,co = np.shape( a )
    a=np.ravel(a)
    b=np.ravel(b)
    V=[ ]
    X=['a']
    for x in a :
        X.append(x)
    for x in b :
        V.append(x)
    T=[]
    c=0 
    for  x in range (1,li*co+1) :
        if x%co!=0:
            T.append(X[x])
        elif x%co==0:
            T.append(X[x])
            T.append(V[c])
            c=c+1
    a = np.array ([T])       
    a.resize(co,co+1)
    Ag=a
    Ag = np.array(Ag,dtype=float)
    return Ag

## Generation de matrices aléatoires
def geneMa(n):
    
    p = 100
    if n >= p :
        p = n
    A = np.random.randint ( 1 , p , size = ( n , n ) )
    A=np.array(A,dtype=float)
    det = np.linalg.det(A)
    while det == 0 :
        A = np.random.randint ( 1 , p , size = ( n , n ) )
        A = np.array(A,dtype=float)
        det = np.linalg.det(A)
    At = np.transpose(A)
    A=At@A
    b = np.random.randint ( 1 , p , size = ( 1 , n ) )
    b = np.array(b,dtype=float)
    return A,b

## Programme de la resolution de système

def resolutionSysTriSup (mat):
    P = np.array(mat)
    li , co = np.shape( P )
    global x
    x = np.zeros( li )
    for i in range ( li-1 , -1 , -1 ) :    
        somme = 0
        for  k  in  range  ( li-1,i, -1 ) :	
            somme = somme + P [ i , k ] * x[ k ]  
        x [ i ] =( P [ i , li ] - somme ) / P [ i , i ]
    return x

def resolutionDiag ( mat ) :
    P = np.array ( mat )
    li , co = np.shape( P )
    y = np.zeros ( ( li ) )
    for i in range ( li ) :
        y [ i ] = ( P [ i , li ] ) / P [ i , i ]
    return y

def resolutionSysTriInf ( mat ) :
    P = np.array ( mat )
    li , co = np.shape( P )
    y = np.zeros ( ( li ) )
    for i in range ( li ) :
        somme = 0
        for k in range ( 0 , i + 1 ) :
            somme = somme + y [ k ] * P [ i , k ]
        y [ i ] = ( P [ i , li ] - somme ) / P [ i , i ]
    return y

##  different-type-de-décomposition-de-matrice

def Cholesky(A) :
    li , co = np.shape ( A )
    L = np.zeros (( li,co ))
    for i in range ( 0 , li ) : 
        somme = 0
        for x in range (i):   
            somme = somme + (L[i,x])**2     
        L [i,i] = (A[i,i] - somme)**0.5
        for c in range ( i+1 , li ) :   
            somme1=0   
            for x in range (i):
                somme1 = somme1 + L[c,x]*L[i,x]    
            L [c,i] = (A[c,i] - somme1)/L[i,i]         
    A = np.array(A,dtype=float)
    L = np.array(L,dtype=float)
    Lt = np.transpose(L)
    Lt = np.array(Lt,dtype=float)
    return  L,Lt

def Cholesky_Python (A):
    L = np.linalg.cholesky(A)
    Lt = np.transpose(L)
    return L,Lt

#----different-type-de-résolution-d'équation-de-matrice---------

def ResolCholesky (A,B) :
    L,Lt = Cholesky(A)
    Ly = concate (L,B)
    y = resolutionSysTriInf ( Ly )
    aLy = concate (Lt,y)
    x = resolutionSysTriSup ( aLy )
    return x

def ResolCholesky_Python (A,B) :
    L,Lt = Cholesky_Python(A)
    Ly = concate (L,B)
    y = resolutionSysTriInf ( Ly )
    aLy = concate (Lt,y)
    x = resolutionSysTriSup ( aLy )
    return x

#------------------Gauss---------------------------
def ReductionGauss ( mat ) :
    li,co = np.shape( mat )
    for k in range ( 0 , li-1 ) :  
        for i in range ( k+1 , li ) :
            g = mat [ i , k ] / mat[ k , k]
            mat[ i , : ] = mat [ i , : ] - g*mat[ k , : ]    
    return mat 

def Gauss ( a , b ) :
    li,co = np.shape( a )
    Ag = concate ( a , b )   
    ReductionGauss (  Ag )
    resolutionSysTriSup ( Ag )
    return x

def Mg(n):
    Rd1,Bd=geneMa(n)
    x = Gauss (Rd1,Bd)
    return x 

#---Génération-de-matrice-aléatoire-pour-la-résolution--------

def AleaResolCholesky (n) :
    Rd1,Bd=geneMa(n)
    x = ResolCholesky (Rd1,Bd)
    return x

def AleaResolCholesky_Python (n) :
    Rd1,Bd=geneMa(n)
    x = ResolCholesky_Python (Rd1,Bd)
    return x

def AllResolv (n):
    Rd1,Bd=geneMa(n)
    x = ResolCholesky(Rd1,Bd)
    xP = ResolCholesky_Python(Rd1,Bd)
    return x,xP

## Exercice 1

def DecompositionGS (A):
    n=len(A)
    R = np.zeros ([n,n])
    Q = np.zeros ([n,n])
    W = np.zeros ([n,n])
    S = 0
    R[0,0] = np.linalg.norm(A[:,0]) 
    Q[:,0] = A[:,0]*(1/R[0,0])
    for j in range (1,n):
        S=0
        for i in range (0,j):
            R[i,j] =  np.vdot(Q[:,i],A[:,j])
            S = S + R[i,j]*Q[:,i]
        W[:,j] = A[:,j]- S     
        R[j,j] = np.linalg.norm(W[:,j]) 
        Q[:,j] = W[:,j]*(1/R[j,j])
    Qt = np.transpose(Q)
    return Q,R,Qt
 


## Exercice 2

def ResolGS(A,b) :
    Q,R,Qt = DecompositionGS (A)
    b = np.transpose(b)
    y = np.dot(Qt,b)
    X = concate(R,y)
    x = resolutionSysTriSup (X)
    return x

def RGS(n):
    A,b = geneMa(n)
    x = ResolGS(A,b)
    E = np.linalg.norm(A@x-np.ravel(b))
    print(E)


## Resolution avec houseolder
    
def ResolHouse(A,b) :
    Q,R,Qt = qr_decomposition(A)
    b = np.transpose(b)
    y = np.dot(Qt,b)
    X = concate(R,y)
    x = resolutionSysTriSup (X)
    return x

def RHouse(n) :
    A,b = geneMa(n)
    x = ResolHouse(A,b)
    E = np.linalg.norm(A@x-np.ravel(b))
    print(E)

## REsolution Qr python

def QrPython(A,b) :
    Q, R = np.linalg.qr(A)
    Qt = np.transpose(Q)
    b = np.transpose(b)
    y = np.dot(Qt,b)
    X = concate(R,y)
    x = resolutionSysTriSup (X)
    return x

def AllResolvListe (n):
    
    Rd1,Bd=geneMa(n)
    
    T0=time.time()
    x = ResolCholesky(Rd1,Bd)
    T1=time.time()
    T = round ( T1-T0 , 5 )
    E = np.linalg.norm(Rd1@x-np.ravel(Bd))

    T0=time.time()
    x = ResolCholesky_Python(Rd1,Bd)
    T1=time.time()
    Ta = round ( T1-T0 , 5 )
    Ea = np.linalg.norm(Rd1@x-np.ravel(Bd))

    T0=time.time()
    x = ResolGS(Rd1,Bd)
    T1=time.time()
    Tb = round ( T1-T0 , 5 )
    Eb = np.linalg.norm(Rd1@x-np.ravel(Bd))

    T0=time.time()
    xg = Gauss(Rd1,Bd)
    T1=time.time()
    Tg = round ( T1-T0 , 5 )
    Eg = np.linalg.norm(Rd1@xg-np.ravel(Bd))

    T0=time.time()
    xh = ResolHouse(Rd1,Bd)
    T1=time.time()
    Th = round ( T1-T0 , 5 )
    Eh = np.linalg.norm(Rd1@xg-np.ravel(Bd))

    T0=time.time()
    xp = QrPython(Rd1,Bd)
    T1=time.time()
    Tp = round ( T1-T0 , 5 )
    Ep = np.linalg.norm(Rd1@xg-np.ravel(Bd))
   
    return T,Ta,Tb,Tg,Th,Tp,E,Ea,Eb,Eg,Eh,Ep




L=[]
L1=[]
L2=[]
L3=[]
L4=[]
L5=[]
L6=[]
L7=[]
L8=[]
L9=[]
L10=[]
L11=[]
L12=[]


for i in range (4,501,25) :
    L.append(float(i))
    T,Ta,Tb,Tg,Th,Tp,E,Ea,Eb,Eg,Eh,Ep = AllResolvListe(i)
    L1.append(T)
    L2.append(E)
    L3.append(Ta)
    L4.append(Ea)
    L5.append(Tb)
    L6.append(Eb)
    L7.append(Tg)
    L8.append(Eg)
    L9.append(Th)
    L10.append(Eh)
    L11.append(Tp)
    L12.append(Ep)

x = L
y1 = L1
y2 = L3
y3 = L5
y4 = L7
y5 = L9
y6 = L11

plt.title("Temps de traitement en fonction de la taille de la matrice")
plt.xlabel("Taille de la matrice")
plt.ylabel("Temps")
plt.plot(x,y1,label="Cholesky")
plt.plot(x,y2,label="Cholesky Python")
plt.plot(x,y3,label="Decomposition QR")
plt.plot(x,y4,label="Gauss")
plt.plot(x,y5,label="QR Houseolder")
plt.plot(x,y6,label="QR Python")
plt.legend()
plt.show()

x = L
y1 = L2
y2 = L4
y3 = L6
y4 = L8
y5 = L10
y6 = L12

plt.title("Erreur en fonction de la taille de la matrice")
plt.xlabel("Taille de la matrice")
plt.ylabel("Erreurs")
plt.plot(x,y1,label="Cholesky")
plt.plot(x,y2,label="Cholesky Python")
plt.plot(x,y3,label="Decomposition QR")
plt.plot(x,y4,label="Gauss")
plt.plot(x,y5,label="QR Houseolder")
plt.plot(x,y6,label="QR Python")

plt.legend()
plt.show()



