# Importation des bibliothèque
import numpy as np
import time 
from math import * 
from matplotlib import pyplot as plt

def geneMa(n):
    p = 100
    if n >= p :
        p = n
    A = np.random.randint ( 1 , 10 , size = ( n , n ) )
    Da = np.random.randint ( 30, 50 , size = ( n , n ) )
    Da = np.diag(np.diag(Da))
    A = A + Da
    A=np.array(A,dtype=float)
    det = np.linalg.det(A)
    while det == 0 :
        A = np.random.randint ( 1 , 10 , size = ( n , n ) )
        A = np.array(A,dtype=float)
        det = np.linalg.det(A)
    At = np.transpose(A)
    A=At@A
    b = np.random.randint ( 1 , 10 , size = ( 1 , n ) )
    b = np.array(b,dtype=float)
    x = np.random.randint ( 1 , 10 , size = ( 1 , n ) )
    x = np.array(x,dtype=float)
    b = np.transpose(b)
    x = np.transpose(x)
    return A,b,x
#----------Méthode direct pour comparaison-------------------
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

#----------------------------deomposition matrice A en M N --------------   

def Rendu_M (A) :
    li , co = np.shape ( A )
    M = np.zeros (( li,co ))
    for i in range ( 0 , li ) :    
        M [i,i] = A[i,i]
        for c in range ( i+1 , li ) :     
            M[i,c] = A[i,c]        
    A = np.array(A,dtype=float)
    M = np.array(M,dtype=float)
    N = M-A
    
    return  M,N

def DecompJacobi (A) :
    li , co = np.shape ( A )
    F = np.zeros (( li,co ))
    for i in range ( 0 , li ) :    
        F [i,i] = A[i,i]
        for c in range ( i+1 , li ) :     
            F[i,c] = A[i,c]         
    E = F-A
    M = np.diag(np.diag(F))
    F = M-F
    N = E + F

    return  M,N


def Decomp_Gauss_Seidel (A) :
    li , co = np.shape ( A )
    F = np.zeros (( li,co ))
    for i in range ( 0 , li ) :    
        F [i,i] = A[i,i]
        for c in range ( i+1 , li ) :     
            F[i,c] = A[i,c]         
    E = F-A
    M = np.diag(np.diag(F))
    F = M-F
    M = M-E
    N =  F
    N = np.array(N,dtype=float)
    M = np.array(M,dtype=float)

    return  M,N

# Question 1

def MIGenerale(M,N,b,x0,epsilon,Nitermax):
    i=0
    e = 1000
    xpas = np.zeros (( len(M),1 ))
    while (i <= Nitermax) and (e>epsilon)  :
        xpas=x0
        x0 = np.linalg.solve(M,N@x0+b)
        e = np.linalg.norm(xpas-x0)
        i=i+1
    if (i >= Nitermax):
        print("ça a divergé")
    return xpas,i,e

def MIrelaxation(A,ome,b,x0,epsilon,Nitermax):
    li , co = np.shape ( A )
    F = np.zeros (( li,co ))
    for i in range ( 0 , li ) :    
        F [i,i] = A[i,i]
        for c in range ( i+1 , li ) :     
            F[i,c] = A[i,c]         
    E = F-A
    M = np.diag(np.diag(F))
    F = M-F
    N =  (((1/ome)-1)*M)+F
    M = ((1/ome)*M-E)
    
    N = np.array(N,dtype=float)
    M = np.array(M,dtype=float)
    i=0
    e = 1000
    xpas = np.zeros (( len(M),1 ))

    while (i <= Nitermax) and (e>epsilon)  :
        xpas=x0
        x0 = np.linalg.solve(M,N@x0+b)
        e = np.linalg.norm(xpas-x0)
        i=i+1
    if (i >= Nitermax):
        print("ça a divergé")
    return xpas,i,e


####Partie 2 


#Question 1

def GeneMat_A_Q1(n):
    A = np.zeros (( n,n ))
    m = np.zeros (( n,n ))
    b = np.zeros (( 1,n ))
    
    for i in range ( 0 , n ) :
        b[0,i] = cos(i/8)
        for j in range ( 0 , n ) :     
            A[i,j] = 1/(12+(3*i-5*j)**2)
            m[i,j] = 3
    A = A-np.diag(np.diag(A))
    m = np.diag(np.diag(m))
    x = np.random.randint ( 1 , 10 , size = ( 1 , n ) )
    x = np.array(x,dtype=float)
    A = A + m
    A = np.array(A,dtype=float)
    b = np.transpose(b)
    x = np.transpose(x)
    return A,b,x




#Question 2

def GeneMat_A_Q2(n):
    A = np.zeros (( n,n ))
    b = np.zeros (( 1,n ))
    
    for i in range ( 0 , n ) :
        b[0,i] = cos(i/8)
        for j in range ( 0 , n ) :     
            A[i,j] = 1/(1+3*abs(i-j))
    x = np.random.randint ( 1 , 10 , size = ( 1 ,n ) )
    x = np.array(x,dtype=float)
    b = np.transpose(b)
    x = np.transpose(x)
    return A,b,x




def Jacobi_compar_Gausseidel_Q1(n,epsilon):
    A,b,x = GeneMat_A_Q1(n)
    M,N = DecompJacobi(A)
    D=[]
    D1=[]
    D2=[]
    sol,nb_ite,ecart = MIGenerale(M,N,b,x,epsilon,Nitermax = 10**100)
    D.append(sol)
    D1.append(nb_ite)
    D2.append(ecart)
    
    M1,N1 = Decomp_Gauss_Seidel(A)
    Q=[]
    Q1=[]
    Q2=[]
    sol1,nb_ite1,ecart1 = MIGenerale(M1,N1,b,x,epsilon,Nitermax = 10**100)
    Q.append(sol1)
    Q1.append(nb_ite1)
    Q2.append(ecart1)

    return D,D1,D2,Q,Q1,Q2


def Jacobi_compar_Gausseidel_Q2(n,epsilon):
    A,b,x = GeneMat_A_Q2(n)
    M,N = DecompJacobi(A)
    D=[]
    D1=[]
    D2=[]
    sol,nb_ite,ecart = MIGenerale(M,N,b,x,epsilon,Nitermax = 10**100)
    D.append(sol)
    D1.append(nb_ite)
    D2.append(ecart)
    
    M1,N1 = Decomp_Gauss_Seidel(A)
    Q=[]
    Q1=[]
    Q2=[]
    sol1,nb_ite1,ecart1 = MIGenerale(M1,N1,b,x,epsilon,Nitermax = 10**100)
    Q.append(sol1)
    Q1.append(nb_ite1)
    Q2.append(ecart1)

    return D1,D2,Q1,Q2



#------------------GRAPHIQUE QUESTION  ------------------------------


def graphique_Q1(ome):
    l1=[]
    l2=[]
    l3=[]
    l4=[]
    l5=[]
    l6=[]
    l7=[]
    l8=[]
    l9=[]
    le=[]
    logle=[]
    t1=[]
    t2=[]
    t3=[]
    tr=[]
    g1=[]
    g2=[]
    g3=[]
    gr=[]
    for i in range (4,15):
        
        a,b,x = GeneMat_A_Q1(100)
        
        T0=time.time()
        M1,N1 = DecompJacobi(a)
        S1,i1,e1 = MIGenerale(M1,N1,b,x,(1*(10**(-i))),Nitermax =100000)
        T1=time.time()
        Ta = round ( T1-T0 , 5 )
        E1 = np.linalg.norm(a@S1-b)
        
        T0=time.time()
        M2,N2 = Decomp_Gauss_Seidel(a)
        S2,i2,e2 = MIGenerale(M2,N2,b,x,(1*(10**(-i))),Nitermax =100000)
        T1=time.time()
        Tb = round ( T1-T0 , 5 )
        E2 = np.linalg.norm(a@S2-b)

        T0=time.time()
        S3,i3,e3 = MIrelaxation(a,ome,b,x,(1*(10**(-i))),Nitermax = 100000)
        T1=time.time()
        Tg = round ( T1-T0 , 5 )
        E3 = np.linalg.norm(a@S3-b)

        T0=time.time()
        xg = Gauss(a,b)
        T1=time.time()
        TR = round ( T1-T0 , 5 )
        E4 = np.linalg.norm(a@xg-np.ravel(b))
        
        l1.append(S1)
        l2.append(i1)
        l3.append(e1)
        l4.append(S2)
        l5.append(i2)
        l6.append(e2)
        l7.append(S3)
        l8.append(i3)
        l9.append(e3)
        le.append((1*(10**(-i))))
        logle.append(log10((1*(10**(-i)))))
        t1.append(Ta)
        t2.append(Tb)
        t3.append(Tg)
        tr.append(TR)
        g1.append(E1)
        g2.append(E2)
        g3.append(E3)
        gr.append(E4)
        print(i)


    x = le
    y1 = l3
    y2 = l6
    y3 = l9
    plt.title("Erreur en fonction de epsilon pour les 3 méthodes du TP")
    plt.xlabel("Logarithme de epsilon")
    plt.ylabel("Logarithme de l'erreur")
    plt.loglog(x,y1,label="Jacobi")
    plt.loglog(x,y2,label="Gauss-seidel")
    plt.loglog(x,y3,label="Relaxation")
    plt.legend()
    plt.show()

    x = le
    y1 = l2
    y2 = l5
    y3 = l8
    plt.title("Itération en fonction de epsilon pour les 3 méthodes du TP")
    plt.xlabel("Logarithme de epsilon")
    plt.ylabel("Itérations")
    plt.semilogx(x,y1,label="Jacobi")
    plt.semilogx(x,y2,label="Gauss-seidel")
    plt.semilogx(x,y3,label="Relaxation")
    plt.legend()
    plt.show()

    x = le
    y1 = t1
    y2 = t2
    y3 = t3
    y4 = tr
    plt.title("Temps en fonction de epsilon pour 4 méthodes")
    plt.xlabel("Logarithme de epsilon")
    plt.ylabel("Temps (s)")
    plt.semilogx(x,y1,label="Jacobi")
    plt.semilogx(x,y2,label="Gauss-seidel")
    plt.semilogx(x,y3,label="Relaxation")
    plt.plot(x,y4,label="Gauss_direct")
    plt.legend()
    plt.show()

    x = le
    y1 = g1
    y2 = g2
    y3 = g3
    y4 = gr
    plt.title("Norme de (A*x-b) en fonction de epsilon pour les 3 méthodes")
    plt.xlabel("Epsilon")
    plt.ylabel("Norme")
    plt.loglog(x,y1,label="Jacobi")
    plt.loglog(x,y2,label="Gauss-seidel")
    plt.loglog(x,y3,label="Relaxation")
    plt.legend()
    plt.show()

def graphique_Q2(ome):
    l1=[]
    l2=[]
    l3=[]
    l4=[]
    l5=[]
    l6=[]
    l7=[]
    l8=[]
    l9=[]
    le=[]
    logle=[]
    t1=[]
    t2=[]
    t3=[]
    tr=[]
    for i in range (4,15):
        
        a,b,x = GeneMat_A_Q2(100)
        
        T0=time.time()
        M1,N1 = DecompJacobi(a)
        S1,i1,e1 = MIGenerale(M1,N1,b,x,(1*(10**(-i))),Nitermax =100000)
        T1=time.time()
        Ta = round ( T1-T0 , 5 )

        T0=time.time()
        M2,N2 = Decomp_Gauss_Seidel(a)
        S2,i2,e2 = MIGenerale(M2,N2,b,x,(1*(10**(-i))),Nitermax =100000)
        T1=time.time()
        Tb = round ( T1-T0 , 5 )

        T0=time.time()
        S3,i3,e3 = MIrelaxation(a,ome,b,x,(1*(10**(-i))),Nitermax = 100000)
        T1=time.time()
        Tg = round ( T1-T0 , 5 )
        
        T0=time.time()
        xg = Gauss(a,b)
        T1=time.time()
        TR = round ( T1-T0 , 5 )
        
        l1.append(S1)
        l2.append(i1)
        l3.append(e1)
        l4.append(S2)
        l5.append(i2)
        l6.append(e2)
        l7.append(S3)
        l8.append(i3)
        l9.append(e3)
        le.append((1*(10**(-i))))
        t1.append(Ta)
        t2.append(Tb)
        t3.append(Tg)
        tr.append(TR)
        print(i)


    x = le
    y1 = l3
    y2 = l6
    y3 = l9
    plt.title("Erreur en fonction de epsilon pour les 3 méthodes du TP")
    plt.xlabel("Logarithme de epsilon")
    plt.ylabel("Logarithme de l'erreur")
    plt.loglog(x,y1,label="Jacobi")
    plt.loglog(x,y2,label="Gauss-seidel")
    plt.loglog(x,y3,label="Relaxation")
    plt.legend()
    plt.show()

    x = le
    y1 = l2
    y2 = l5
    y3 = l8
    plt.title("Iteration en fonction de epsilon pour les 3 méthodes du TP")
    plt.xlabel("Logarithme de epsilon")
    plt.ylabel("Iterations")
    plt.semilogx(x,y1,label="Jacobi")
    plt.semilogx(x,y2,label="Gauss-seidel")
    plt.semilogx(x,y3,label="Relaxation")
    plt.legend()
    plt.show()

    x = le
    y1 = t1
    y2 = t2
    y3 = t3
    y4 = tr
    plt.title("Temps en fonction de epsilon pour 4 méthodes")
    plt.xlabel("Logarithme de epsilon")
    plt.ylabel("Temps (s)")
    plt.semilogx(x,y1,label="Jacobi")
    plt.semilogx(x,y2,label="Gauss-seidel")
    plt.semilogx(x,y3,label="Relaxation")
    plt.semilogx(x,y4,label="Gauss_direct")
    plt.legend()
    plt.show()
    
def Relexation_omega():   
    l1=[]
    l2=[]
    l3=[]
    l4=[]
    X = np.linspace(0.5,1.5,25)
    for i in X:
        ome = (i)
        a,b,x = GeneMat_A_Q1(100)
        s,it,e = MIrelaxation(a,ome,b,x,(0.0001),Nitermax = 100000)
        l1.append(s)
        l2.append(it)
        l3.append(e)
        l4.append(ome)


    x = l4
    y1 = l2
    plt.title("Vitesse de resolution en fonction de omega")
    plt.xlabel("Omega")
    plt.ylabel("Logarithme de itérations")
    plt.plot(x,y1,label="Omega")
    plt.legend()
    plt.show()

    l1=[]
    l2=[]
    l3=[]
    l4=[]
    X = np.linspace(0.25,1.25,50)
    for i in X:
        ome = (i)
        a,b,x = GeneMat_A_Q2(100)
        s,it,e = MIrelaxation(a,ome,b,x,(0.0001),Nitermax = 100000)
        l1.append(s)
        l2.append(it)
        l3.append(e)
        l4.append(ome)


    x = l4
    y1 = l2
    plt.title("Vitesse de resolution en fonction de omega")
    plt.xlabel("Omega")
    plt.ylabel("Itérations")
    plt.plot(x,y1,label="Omega")
    plt.legend()
    plt.show()

    
    
def graphique_MATRICE_Aletatoire(ome):
    l1=[]
    l2=[]
    l3=[]
    l4=[]
    l5=[]
    l6=[]
    l7=[]
    l8=[]
    l9=[]
    le=[]
    logle=[]
    t1=[]
    t2=[]
    t3=[]
    tr=[]
    for i in range (4,15):
        
        a,b,x = geneMa(50)
        
        T0=time.time()
        M1,N1 = DecompJacobi(a)
        S1,i1,e1 = MIGenerale(M1,N1,b,x,(1*(10**(-i))),Nitermax =100000)
        T1=time.time()
        Ta = round ( T1-T0 , 5 )

        T0=time.time()
        M2,N2 = Decomp_Gauss_Seidel(a)
        S2,i2,e2 = MIGenerale(M2,N2,b,x,(1*(10**(-i))),Nitermax =100000)
        T1=time.time()
        Tb = round ( T1-T0 , 5 )

        T0=time.time()
        S3,i3,e3 = MIrelaxation(a,ome,b,x,(1*(10**(-i))),Nitermax = 100000)
        T1=time.time()
        Tg = round ( T1-T0 , 5 )

        T0=time.time()
        xg = Gauss(a,b)
        T1=time.time()
        TR = round ( T1-T0 , 5 )
        
        l1.append(S1)
        l2.append(i1)
        l3.append(e1)
        l4.append(S2)
        l5.append(i2)
        l6.append(e2)
        l7.append(S3)
        l8.append(i3)
        l9.append(e3)
        le.append((1*(10**(-i))))
        t1.append(Ta)
        t2.append(Tb)
        t3.append(Tg)
        tr.append(TR)
        print(i)


    x = le
    y1 = l3
    y2 = l6
    y3 = l9
    plt.title("Erreur en fonction de epsilon pour les 3 methodes")
    plt.xlabel("Logarithme de epsilon")
    plt.ylabel("Logarithme de l'erreur")
    plt.loglog(x,y1,label="Jacobi")
    plt.loglog(x,y2,label="Gauss")
    plt.loglog(x,y3,label="Relaxation")
    plt.legend()
    plt.show()

    x = le
    y1 = l2
    y2 = l5
    y3 = l8
    plt.title("Iteration en fonction de epsilon pour les 3 methodes")
    plt.xlabel("Logarithme de epsilon")
    plt.ylabel("Itération")
    plt.semilogx(x,y1,label="Jacobi")
    plt.semilogx(x,y2,label="Gauss_seidel")
    plt.semilogx(x,y3,label="Relaxation")
    plt.legend()
    plt.show()

    x = le
    y1 = t1
    y2 = t2
    y3 = t3
    y4 = tr
    plt.title("Temps en fonction de epsilon pour les 3 methodes")
    plt.xlabel("Logarithme de epsilon")
    plt.ylabel("Temps(s)")
    plt.semilogx(x,y1,label="Jacobi")
    plt.semilogx(x,y2,label="Gauss_seidel")
    plt.semilogx(x,y3,label="Relexation")
    plt.semilogx(x,y4,label="Gauss_direct")
    plt.legend()
    plt.show()   


#graphique_Q1(1.2)
#graphique_Q2(1.2)
##Relexation_omega()
#graphique_MATRICE_Aletatoire(1.2)
#graphique_Q2(0.8)























    
