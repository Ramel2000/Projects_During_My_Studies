import numpy as np
import matplotlib.pyplot as plt
import time
global TC
global Ttotal
global y 
import copy 

Aaug= np.array([[1,1,1,1],[2,4,-7,2],[-1,-1,0,-3],[1,-1,4,1]])
B = np.array([[1,1,2,-8]])

##li,co=np.shape(Aaug)
##print(Aaug)
##print(li,co)
##Rd = np.random.randint(low=1, high=100, size=(10,10))
##Bd = np.random.randint(low=1, high=100, size=(1,10))
##Ag=np.concatenate((Rd, Bd.T), axis=1)
##
##ReductionGauss (Ag)
#numpy.linalg.norm

##Rd = np.random.randint(low=1, high=100, size=(10,10))
##Bd = np.random.randint(low=1, high=100, size=(1,10))
##print(Rd)
##print(Bd)
##Ad = np.concatenate((Rd, Bd.T), axis=1)
##
##print (Ad)

    
def ReductionGauss ( mat ) :
	
    li,co = np.shape( mat )

    for k in range ( 0 , li-1 ) :
        
        for i in range ( k+1 , li ) :

            g = mat [ i , k ] / mat[ k , k]

            mat[ i , : ] = mat [ i , : ] - g*mat[ k , : ]
            
    return mat 





def resolutionSysTriSup ( mat ) :
    
    P = np.array(mat)
	
    li , co = np.shape( P )
    
    global x
    
    x = np.zeros( li )
    
    for i in range ( li-1 , -1 , -1 ) :
	    
        somme = 0
        
        for  k  in  range  ( li-1,i, -1 ) :
		
            somme = somme + P [ i , k ] * x[ k ]
            
        x [ i ] =( P [ i , li ] - somme ) / P [ i , i ]
    
    print ( 'Les solutions sont' , x )
    return x




def resolutionSysTriInf ( mat ) :
    P = np.array ( mat )
    li , co = np.shape( P )
    global y
    y = np.zeros((li))
    for i in range ( li ) :
        somme = 0
        for k in range ( 0 , i+1 ) :
            somme = somme + y[k]*P[i,k]
        y[i]=(P[i,li]-somme)/P[i,i]
    print(y)
    return y

##def resolutionSysTriInf(mat):
##
##    n=mat.shape[0]
##    y=np.zeros(n)
##    for i in range(n):
##        somme=0
##        for k in range(i):
##            somme=somme+y[k]*mat[i,k]
##        y[i]=(mat[i,-1]-somme)/mat[i,i]
##    return y

            



##def resolutionSysTriSup(mat):
##    global x
##    li,co=mat.shape
##    if co !=li+1:
##        print('pas une matrice augmentée')
##        return
##    x=np.zeros(li)
##    for i in range(li-1,-1,-1):
##        somme=0
##        for k in range(i,li):
##            somme=somme+x[k]*mat[i,k]
##        x[i]=(mat[i,-1]-somme)/mat[i,i]
##    print ( 'Les solutions sont' , x )
##    return x        




def Gauss ( a , b ) :
    global x
    li,co = np.shape( a )
##    print(a)
##    a=np.ravel(a)
##    print(a)
##    print( b)
##    b=np.ravel(b)
##    V=[]
##    X=['a']
##    for x in a :
##        X.append(x)
##    for x in b :
##        V.append(x)
##    print(X)
##    print(V)
##    T=[]
##    print(li)
##
##    c=0 
##    for  x in range (1,li*co+1) :
##
##        if x%co!=0:
##            T.append(X[x])
##
##        elif x%co==0:
##            T.append(X[x])
##            T.append(V[c])
##            c=c+1
##
##    a = np.array ([T])       
##    print(T)
##    print(a)
##    a.resize(co,co+1)
##    print(a)
##    Ag=a
    Ag = np.concatenate ( ( a , b.T ) , axis = 1 )  
    ReductionGauss (  Ag )
    resolutionSysTriSup ( Ag )
    print('Les resultat pour x avec Gauss sont:',x)
    Erreur = np.linalg.norm(a@x-np.ravel(b))
    print(Erreur)
    return x


def Q4 ( n ) :
    global Rd1
    global Bd1
	
    p = 100
    
    if n >= p :
	    
        p = n
        
    print ( 'La taille de la matrice est' , n )
    
    Rd = np.random.randint ( 1 , p , size = ( n , n ) )
    Rd1=np.array(Rd,dtype=float)
    
    Bd = np.random.randint ( 1 , p , size = ( 1 , n ) )
    Bd1=np.array(Bd,dtype=float)
    
    print ( '""""""""""""""""' )
    
    Gauss( Rd , Bd )










def TQ( n ) :
	global TC
	TC = []
	for i in range ( 10 ) :
		TM=[]
		T0=time.time()
		Q4( n )
		T1 = time.time()
		T = round ( T1-T0 , 5 )
		print( T )
		TM.append( n )
		TC.append( T )
		TT=[TM,TC]
		print('La matrice de taille' , TT[[0][0]] , 'prend' , TT[[1][0]] , 'seconde à etre calculé' )    
def FMoyenneTemps () :
	global TC
	global Ttotal
	Ttot = 0
	for i in TC :	
		Ttot = Ttot + i	
	Ttotal = Ttot / len ( TC )   
def TempsMoyen(n):
	global TC
	global Ttotal
	TQ(n)
	FMoyenneTemps()
	print(Ttotal)
	return Ttotal


##MTT10=TempsMoyen(10)
##MTT100=TempsMoyen(100)
##MTT500=TempsMoyen(500)
##MTT1000=TempsMoyen(1000)
##MTT1500=TempsMoyen(1500)
##L=[]
##L.append(MTT10)
##L.append(MTT100)
##L.append(MTT500)
##L.append(MTT1000)
##L.append(MTT1500)
##
##F=[10,100,500,1000,1500]
##
##x = L
##y = F
##
##plt.title("Temps de traitement en fonction de la taille de la matrice")
##plt.ylabel("Taille de la matrice")
##plt.xlabel("Temps")
##plt.plot(x,y)
##plt.show()




def erreur ( n ) :
    global Rd1
    global Bd1
    global x
    p = 100
    if n >= p :
    	p = n
    print('La taille de la matrice est',n)
    Q4(n)
    Erreur = np.linalg.norm(Rd1@Gauss(Rd1,Bd1)-np.ravel(Bd1))

    return Erreur


##erreur (4)



##MTT10=erreur(10)
##MTT100=erreur(100)
##MTT500=erreur(500)
##MTT1000=erreur(1000)
##MTT1500=erreur(1500)
##L=[]
##L.append(MTT10)
##L.append(MTT100)
##L.append(MTT500)
##L.append(MTT1000)
##L.append(MTT1500)
##print(L)
##F=[10,100,500,1000,1500]
##
##x = L
##y = F
##
##plt.title("Erreur en fonction de la taille de la matrice")
##plt.ylabel("Taille de la matrice")
##plt.xlabel("Erreur commises")
##plt.plot(x,y)
##plt.show()


def DecompositionLU ( mat ) :
    li , co = np.shape ( mat )
    L0=  []
    for k in range ( 0 , li-1 ) :
        for i in range ( k+1 , li ) :
            g = mat [ i , k ] / mat[ k , k]
            mat[ i , : ] = mat [ i , : ] - g*mat [ k , : ]
            L0.append ( g )
    I = np.identity ( li )
    d = 0
    m = 0
    for i in range ( co ) :
        m = m + 1        
        for k in range ( li-1-i ) :  
            I [ k + m , i ] = float ( L0 [ d ] ) 
            d = d + 1
    mat = np.array(mat,dtype=float)
    return mat , I




def concate (a,b):
    
    li,co = np.shape( a )
    a=np.ravel(a)
    b=np.ravel(b)
    V=[]
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

def ResolutionLU ( L , U , B ) :
    global y
    global x
    Ag = np.concatenate ( ( L , B.T ) , axis = 1 )
    li,co = np.shape(L)
    y = resolutionSysTriInf ( Ag )
    Ag1 = concate ( U , y )
    resolutionSysTriSup ( Ag1 )
    print('Les resultat pour x avec LU sont:',x)

def Resolution_LU(n) :
    global x
    global y
    p = 100
    if n >= p :    
        p = n   
    print ( 'La taille de la matrice est' , n )
    Rd = np.random.randint ( 1 , p , size = ( n , n ) )
    Rd1=np.array(Rd,dtype=float)
    Bd = np.random.randint ( 1 , p , size = ( 1 , n ) )
    Bd1=np.array(Bd,dtype=float)

    U1 , L1 = DecompositionLU ( Rd1 )

    ResolutionLU ( L1 , U1 , Bd1 )
    Erreur = np.linalg.norm(L1@U1@x-np.ravel(Bd1))

    print(Erreur)
    print('x:',x)
    return (Erreur)

def TLU( n ) :
    T0=time.time()
    Resolution_LU( n )
    T1 = time.time()
    T = round ( T1-T0 , 5 )
    return T

##MTT10=TLU(10)
##MTT100=TLU(100)
##MTT500=TLU(500)
##MTT1000=TLU(1000)
##MTT1500=TLU(1500)
##L=[]
##L.append(MTT10)
##L.append(MTT100)
##L.append(MTT500)
##L.append(MTT1000)
##L.append(MTT1500)
##print(L)
##F=[10,100,500,1000,1500]
##
##x = L
##y = F

##plt.title("Temps de traitement en fonction de la taille de la matrice")
##plt.ylabel("Taille de la matrice")
##plt.xlabel("Temps de traitement")
##plt.plot(x,y)
##plt.show()



def ReductionGaussPartiel ( mat ) :
    print(mat)
	
    li,co = np.shape( mat )
    
    for k in range ( 0 , li - 1 ) :
        
        for i in range ( k + 1 , li ) :

            for q in range (k,li-1) :
                sup = k
                if abs(mat[q,k]) < abs(mat [q+1,k]):
                    sup = q+1
            if sup!=k :
                Nmat=np.array(mat)
                Nmat1=np.array(mat)

                mat [k,:]=Nmat[sup,:]
                mat[sup,:]=Nmat1[k,:]
                sup=k
                
                

            g = mat [ i , k ] / mat[ k , k]

            mat[ i , : ] = mat [ i , : ] - g*mat[ k , : ]
    print(mat)
            
    return mat 


def GaussChoixPivotPartiel ( A , B ) :
    global x
    Ag = np.concatenate ( ( A , B.T ) , axis = 1 )
    Ag=np.array(Ag,dtype=float)
    ReductionGaussPartiel (  Ag )
    print(Ag)
    resolutionSysTriSup ( Ag )
    
    print('Les resultat pour x avec Gausspartiel sont:',x)


    return x

    
def Resolution_PIVOTPARTIEL(n) :
    global x
    global y
    p = 100
    if n >= p :    
        p = n   
    print ( 'La taille de la matrice est' , n )
    Rd = np.random.randint ( 1 , p , size = ( n , n ) )
    Rd1=np.array(Rd,dtype=float)
    Bd = np.random.randint ( 1 , p , size = ( 1 , n ) )
    Bd1=np.array(Bd,dtype=float)

    GaussChoixPivotPartiel(Rd1 , Bd1 )
    Erreur = np.linalg.norm(Rd1@x-np.ravel(Bd1))

    print(Erreur)
    print('x:',x)
    return (Erreur)

def TPP( n ) :
    T0=time.time()
    Resolution_PIVOTPARTIEL( n )
    T1 = time.time()
    T = round ( T1-T0 , 5 )
    return T

##MTT10=TPP(10)
##MTT100=TPP(100)
##MTT500=TPP(200)
##MTT1000=TPP(300)
##MTT1500=TPP(400)
##L=[]
##L.append(MTT10)
##L.append(MTT100)
##L.append(MTT500)
##L.append(MTT1000)
##L.append(MTT1500)
##print(L)
##F=[10,100,200,300,400]
##
##x = L
##y = F
##
##plt.title("Temps de traitement en fonction de la taille de la matrice")
##plt.ylabel("Taille de la matrice")
##plt.xlabel("Temps de traitement")
##plt.plot(x,y)
##plt.show()





def ReductionGaussTotal ( mat ) :
    global Ini
    Ini=[]
    li,co = np.shape( mat )
    for k in range ( 0 , li - 1 ) :
        for i in range ( k +1 , li ) :
            for q in range (k,li-1) :
                sup1=q
                if abs(mat[k,q]) < abs(mat [k,q+1]):
                    sup1 = q+1
                if sup1 !=q :
                    Ini.append([q,sup1])
                    Nmat1=np.array(mat)
                    Nmatt1=np.array(mat)
                    mat [:,q]=Nmat1[:,sup1]
                    mat[:,sup1]=Nmatt1[:,q]
                    sup1=q    
                sup = k
                if abs(mat[q,k]) < abs(mat [q+1,k]):
                    sup = q+1     
            if sup!=k :
                Nmat=np.array(mat)
                Nmatt=np.array(mat)
                mat [k,:]=Nmat[sup,:]
                mat[sup,:]=Nmatt[k,:]
                sup=k
            g = mat [ i , k ] / mat[ k , k]
            mat[ i , : ] = mat [ i , : ] - g*mat[ k , : ]
    return mat

def GaussChoixPivotTotal ( A , B ) :
    global x
    global Ini
    global Fin
    Ag = np.concatenate ( ( A , B.T ) , axis = 1 )
    Ag=np.array(Ag,dtype=float)
    ReductionGaussTotal (  Ag )
    print(Ag)
    resolutionSysTriSup ( Ag )
    for e in reversed(Ini):
        h=copy.copy(x[e[0]])
        x[e[0]]=x[e[1]]
        x[e[1]]=h        
    print('Les resultat pour x avec Gausstotal sont:',x)
    Erreur = np.linalg.norm(A@x-np.ravel(B))
    print(Erreur)
    return x



def Resolution_PIVOTTOTAL(n) :
    global x
    global y
    p = 100
    if n >= p :    
        p = n   
    print ( 'La taille de la matrice est' , n )
    Rd = np.random.randint ( 1 , p , size = ( n , n ) )
    Rd1=np.array(Rd,dtype=float)
    Bd = np.random.randint ( 1 , p , size = ( 1 , n ) )
    Bd1=np.array(Bd,dtype=float)

    GaussChoixPivotTotal(Rd1 , Bd1 )
    Erreur = np.linalg.norm(Rd1@x-np.ravel(Bd1))

    print(Erreur)
    print('x:',x)
    return (Erreur)

MTT10=Resolution_PIVOTTOTAL(10)
MTT100=Resolution_PIVOTTOTAL(100)
MTT500=Resolution_PIVOTTOTAL(200)
MTT1000=Resolution_PIVOTTOTAL(250)

L=[]
L.append(MTT10)
L.append(MTT100)
L.append(MTT500)
L.append(MTT1000)

print(L)
F=[10,100,200,250]

x = L
y = F

plt.title("Temps de traitement en fonction de la taille de la matrice")
plt.ylabel("Taille de la matrice")
plt.xlabel("Temps de traitement")
plt.plot(x,y)
plt.show()

def TPP( n ) :
    T0=time.time()
    Resolution_PIVOTTOTAL( n )
    T1 = time.time()
    T = round ( T1-T0 , 5 )
    return T

##MTT10=TPP(10)
##MTT100=TPP(100)
##MTT500=TPP(200)
##MTT1000=TPP(250)
##
##L=[]
##L.append(MTT10)
##L.append(MTT100)
##L.append(MTT500)
##L.append(MTT1000)
##
##print(L)
##F=[10,100,200,250]
##
##x = L
##y = F
##
##plt.title("Temps de traitement en fonction de la taille de la matrice")
##plt.ylabel("Taille de la matrice")
##plt.xlabel("Temps de traitement")
##plt.plot(x,y)
##plt.show()
