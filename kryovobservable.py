import numpy as np
from scipy.linalg import solve_triangular
import csv
from joblib import Parallel, delayed
import time
import math
def add(a,b):
    c = {}
    nz_idx = set(a.keys()).union(set(b.keys())) 
    for i in nz_idx:
        c[i] = 0
        if i in a:
            c[i] += a[i]
        if i in b:
            c[i] += b[i]
    return c

def innerprod(a,b):
    c = 0.0
    nz_idx = set(a.keys()).intersection(set(b.keys())) 
    for  i in nz_idx:
        c += a[i]*b[i]
    return c
def norm(a):
    sqa=math.sqrt(innerprod(a,a))
    return sqa
def cleandict(a):
    result={}
    for x, y in a.items(): 
        if(abs(y) > 1e-7):result.update({x:y})
    return result
#####################################################################################################
def getAofE(A,Cm,Hm,N,maxm):## A is a vector as H0. Cm np (m,1) where m =5 matrix and Hm ={H0...H4} where the dictionaries in dictionary Hm 
    sumA=0.0
    m=len(Hm)
#    print(m)
 #   sumAnp=np.zeros((1,maxm+1))
    for i in range(maxm+1):
        #print(i)
#        print(innerprod(A,Hm['H'+str(i)]))
        sumA=sumA+Cm[i][0]*innerprod(A,Hm['H'+str(i)])

    return sumA
        
def getAofEnp(A,Cm,Hm,N,maxm):## A is a vector as H0. Cm np (m,1) where m =5 matrix and Hm ={H0...H4} where the dictionaries in dictionary Hm 
    sumA=0.0
    m=len(Hm)
    sumAnp=np.zeros((1,maxm+1))
    for i in range(maxm+1):
        print(i)
        print(innerprod(A,Hm['H'+str(i)]))
        sumA=sumA+Cm[i][0]*innerprod(A,Hm['H'+str(i)])
        sumAnp[0][i]=sumA

    return sumAnp
        
def obsnp(Cm,listtr,N,m):
   #listtr={'listtrzH':listtrzH,'listtrzzH'....'listtrxH':...} where listtrzH=[nparray1..m] noted 1st to Hm0 no need
    listtrzH=listtr['listtrzH']
    listtrzzH=listtr['listtrzzH']
#    print('\nlisttrzzH={}\n'.format(listtrzzH))
    listtrxH=listtr['listtrxH']
#    print('\nlenoflisttrzH{}\n'.format(listtrzzH[0]))
    listMzCm=[]
    listMxCm=[]
    listCzzstCm=[]
    listsumMz=[]
    listsumMx=[]
    listsumCzzst=[]
    sumMz=0.0
    sumMx=0.0
    sumCzzst=0.0
#    print('\nTr(x*H1)={},C1={},C1.*Tr()={},sum(previous)={}\n'.format(listtrxH[1],Cm[1][0],Cm[1][0]*listtrxH[1],np.sum(Cm[1][0]*listtrxH[1])))
#    print('\nTr(x*H2)={},C2={},C2.*Tr()={},sum(previous)={}\n'.format(listtrxH[2],Cm[2][0],Cm[2][0]*listtrxH[2],np.sum(Cm[2][0]*listtrxH[2])))

    for i in range(m):
        listMzCm.append(np.sum(Cm[i][0]*listtrzH[i])) ## /sum_sites CmTr(Hm*)
        listMxCm.append(np.sum(Cm[i][0]*listtrxH[i]))
        listCzzstCm.append(np.sum(Cm[i][0]*listtrzzH[i]))
        
    for i in range(m):
        sumMz=sumMz+listMzCm[i]
#        print('\ni={},sumMz={}\n'.format(i,sumMz/N))
        listsumMz.append(sumMz/N)
        sumMx=sumMx+listMxCm[i]
        listsumMx.append(sumMx/N)
        sumCzzst=sumCzzst+listCzzstCm[i]
        listsumCzzst.append(sumCzzst/N)
    return listsumMx,listsumMz,listsumCzzst 
    
    


# In[120]:


def Mx(Cm,Hm,N,maxm):
    sumMx=0.0
    allI=''
    for p in range(N):allI=allI+'I'
    
    for i in range(N):
        if i==0:A={'x'+allI[i+1:]:1.0}
        elif i==(N-1):A={allI[:i]+'x':1.0}
        else:A={allI[:i]+'x'+allI[i+1:]:1.0}
        sumMx=sumMx+getAofE(A,Cm,Hm,N,maxm)
#        print(A)
#        print(sumMx)
    sumMx=sumMx/N
    return sumMx
        

    


# In[121]:

def obsoneforall(Cm,Hm,N,Nx,maxm,dr):
    Ny=int(N/Nx)
    sumMx=0.0
    sumMz=0.0
    sumCzzc=0.0
    sumCzz=0.0
    sitesxy=np.zeros((Ny, Nx))
    allI=''
    for p in range(N):allI=allI+'I'
    for j in range(0,N-Nx+1,Nx):
        for i in range(j,Nx+j,1):
            if (i-dr)< j and (i+dr)>(Nx+j-1):exp=0.0;exp1=0.0
            elif (i-dr) < j:
                if(dr==0):Ax={'x'+allI[i+1:]:1.0}
                if(dr==0):Az={'z'+allI[i+1:]:1.0}
                A={allI[:i]+'z'+allI[i+1:]:1.0}
                B={allI[:i+dr]+'z'+allI[i+dr+1:]:1.0}
                C={allI[:i]+'z'+allI[:dr-1]+'z'+allI[i+dr+1:]:1.0}
#             print(A.keys())
#             print(B.keys())
#             print(C.keys())
                exp1=getAofE(C,Cm,Hm,N,maxm)
                exp=exp1-getAofE(A,Cm,Hm,N,maxm)*getAofE(B,Cm,Hm,N,maxm)
            elif (i+dr)>(Nx+j-1):
                if(dr==0):Ax={allI[:i]+'x':1.0}
                if(dr==0):Az={allI[:i]+'z':1.0}
                A={allI[:i-dr]+'z'+allI[i-dr+1:]:1.0}
                B={allI[:i]+'z'+allI[i+1:]:1.0}
                C={allI[:i-dr]+'z'+allI[:dr-1]+'z'+allI[i+1:]:1.0}
                exp1=getAofE(C,Cm,Hm,N,maxm)
                exp=exp1-getAofE(A,Cm,Hm,N,maxm)*getAofE(B,Cm,Hm,N,maxm)
            else:
                if(dr==0):Ax={allI[:i]+'x'+allI[i+1:]:1.0}
                if(dr==0):Az={allI[:i]+'z'+allI[i+1:]:1.0}
                Af={allI[:i]+'z'+allI[i+1:]:1.0}
                Bf={allI[:i+dr]+'z'+allI[i+dr+1:]:1.0}
                Cf={allI[:i]+'z'+allI[:dr-1]+'z'+allI[i+dr+1:]:1.0}
                Bb={allI[:i-dr]+'z'+allI[i-dr+1:]:1.0}
                Cb={allI[:i-dr]+'z'+allI[:dr-1]+'z'+allI[i+1:]:1.0}           
                exp1=getAofE(Cf,Cm,Hm,N,maxm)+getAofE(Cb,Cm,Hm,N,maxm)
                exp=exp1-getAofE(Af,Cm,Hm,N,maxm)*getAofE(Bf,Cm,Hm,N,maxm)-getAofE(Af,Cm,Hm,N,maxm)*getAofE(Bb,Cm,Hm,N,maxm)
            if(dr==0):sumMx=sumMx+getAofE(Ax,Cm,Hm,N,maxm)
            if(dr==0):sitesxy[int(i/Nx),int(i%Nx)]=getAofE(Az,Cm,Hm,N,maxm)
            if(dr==0):sumMz=sumMz+sitesxy[int(i/Nx),int(i%Nx)]
            sumCzzc=sumCzzc+exp
            sumCzz=sumCzz+exp1
    sumMx=sumMx/N
    sumMz=sumMz/N
    sumCzzc=sumCzzc/(2.0*(N-dr))
    sumCzz=sumCzz/(2.0*(N-dr))
    return sumMx,sumMz,sitesxy,sumCzzc,sumCzz
def Czzdr(Cm,Hm,N,Nx,dr,maxm):## to be completed
    sumCzzc=0.0
    sumCzz=0.0
    allI=''
    for p in range(N):allI=allI+'I'
    
    for j in range(0,N-Nx+1,Nx):
        for i in range(j,Nx+j,1): ## 
#            print(j)
#            print(i)
            if (i-dr)< j and (i+dr)>(Nx+j-1):exp=0.0;exp1=0.0
            elif (i-dr) < j:
                A={allI[:i]+'z'+allI[i+1:]:1.0}
                B={allI[:i+dr]+'z'+allI[i+dr+1:]:1.0}
                C={allI[:i]+'z'+allI[:dr-1]+'z'+allI[i+dr+1:]:1.0}
#             print(A.keys())
#             print(B.keys())
#             print(C.keys())
                exp1=getAofE(C,Cm,Hm,N,maxm)
                exp=exp1-getAofE(A,Cm,Hm,N,maxm)*getAofE(B,Cm,Hm,N,maxm)
            elif (i+dr)>(Nx+j-1):
                A={allI[:i-dr]+'z'+allI[i-dr+1:]:1.0}
                B={allI[:i]+'z'+allI[i+1:]:1.0}
                C={allI[:i-dr]+'z'+allI[:dr-1]+'z'+allI[i+1:]:1.0}
                exp1=getAofE(C,Cm,Hm,N,maxm)
                exp=exp1-getAofE(A,Cm,Hm,N,maxm)*getAofE(B,Cm,Hm,N,maxm)
#             print(A.keys())
#             print(B.keys())
#             print(C.keys())

            else:
                Af={allI[:i]+'z'+allI[i+1:]:1.0}
                Bf={allI[:i+dr]+'z'+allI[i+dr+1:]:1.0}
                Cf={allI[:i]+'z'+allI[:dr-1]+'z'+allI[i+dr+1:]:1.0}
                Bb={allI[:i-dr]+'z'+allI[i-dr+1:]:1.0}
                Cb={allI[:i-dr]+'z'+allI[:dr-1]+'z'+allI[i+1:]:1.0}
#             print(Af.keys())
#             print(Bf.keys())
#             print(Cf.keys())
#             print(Bb.keys())
#             print(Cb.keys())             
                exp1=getAofE(Cf,Cm,Hm,N,maxm)+getAofE(Cb,Cm,Hm,N,maxm)
                exp=exp1-getAofE(Af,Cm,Hm,N,maxm)*getAofE(Bf,Cm,Hm,N,maxm)-getAofE(Af,Cm,Hm,N,maxm)*getAofE(Bb,Cm,Hm,N,maxm)
            sumCzzc=sumCzzc+exp
            sumCzz=sumCzz+exp1
#         print(exp)
#        print(sumMx)
    sumCzzc=sumCzzc/(2.0*(N-dr))
    sumCzz=sumCzz/(2.0*(N-dr))
    return sumCzzc,sumCzz
def Mz(Cm,Hm,N,Nx,Ny,maxm):
    sumMx=0.0
    allI=''
    for p in range(N):allI=allI+'I'
    sitesxy=np.zeros((Nx, Ny))
    print('\nwriting on-site\n')
#        print(getAofE(A,Cm,Hm,N))
    for i in range(N):
        if i==0:A={'z'+allI[i+1:]:1.0}
        elif i==(N-1):A={allI[:i]+'z':1.0}
        else:A={allI[:i]+'z'+allI[i+1:]:1.0}
        sitesxy[int(i%Nx),int(i/Nx)]=getAofE(A,Cm,Hm,N,maxm)
        sumMx=sumMx+sitesxy[int(i%Nx),int(i/Nx)]

        
    sumMx=sumMx/N
    
    return sumMx,sitesxy
def Czz(Cm,Hm,N):
    sumCzzc=0.0
    sumCzz=0.0
    allI=''
    for p in range(N):allI=allI+'I'
    
    for i in range(N):
        if i==0:
            A={'z'+allI[i+1:]:1.0}
            B={'Iz'+allI[i+2:]:1.0}
            C={'zz'+allI[i+2:]:1.0}
            exp1=getAofE(C,Cm,Hm,N,maxm)
            exp=exp1-getAofE(A,Cm,Hm,N,maxm)*getAofE(B,Cm,Hm,N,maxm)
        elif i==(N-1):
            A={allI[:i-1]+'zI':1.0}
            B={allI[:i-1]+'Iz':1.0}
            C={allI[:i-1]+'zz':1.0}
            exp1=getAofE(C,Cm,Hm,N,maxm)
            exp=exp1-getAofE(A,Cm,Hm,N,maxm)*getAofE(B,Cm,Hm,N,maxm)
        else:
            Af={allI[:i]+'z'+allI[i+1:]:1.0}
            Bf={allI[:i+1]+'z'+allI[i+2:]:1.0}
            Cf={allI[:i]+'zz'+allI[i+2:]:1.0}
            Bb={allI[:i-1]+'z'+allI[i:]:1.0}
            Cb={allI[:i-1]+'zz'+allI[i+1:]:1.0}
            exp1=getAofE(Cf,Cm,Hm,N,maxm)+getAofE(Cb,Cm,Hm,N,maxm)
            exp=exp1-getAofE(Af,Cm,Hm,N,maxm)*getAofE(Bf,Cm,Hm,N,maxm)-getAofE(Af,Cm,Hm,N,maxm)*getAofE(Bb,Cm,Hm,N,maxm)
  
            ## Ab=Af
            
            
        sumCzzc=sumCzzc+exp
        sumCzz=sumCzz+exp1
 #       print(exp)
#        print(sumMx)
    sumCzzc=sumCzzc/(2.0*(N-1))
    sumCzz=sumCzz/(2.0*(N-1))
    return sumCzzc,sumCzz

#12 13  14  15
#8  9  10  11
#4  5   6   7
#0  1   2   3
    
def Czzdrmid(Cm,Hm,N,Nx,maxm):## to be completed 
  
    sumCzzc=0.0
    sumCzz=0.0
    Ny=int(N/Nx)
    midy=(math.ceil(Ny/2)-1)*Nx
    listCzzcdr=[]
    listCzzdr=[]
    allI=''
    for p in range(N):allI=allI+'I'
    i=midy+2
    for dr in range(1,Nx-2,1): ## 

        A={allI[:i]+'z'+allI[i+1:]:1.0}
        B={allI[:i+dr]+'z'+allI[i+dr+1:]:1.0}
        C={allI[:i]+'z'+allI[:dr-1]+'z'+allI[i+dr+1:]:1.0}

        
        exp1=getAofE(C,Cm,Hm,N,maxm)
        listCzzdr.append(exp1)
        exp=exp1-getAofE(A,Cm,Hm,N,maxm)*getAofE(B,Cm,Hm,N,maxm)
        listCzzcdr.append(exp)
 #       print(exp)
         

    return listCzzcdr,listCzzdr

    