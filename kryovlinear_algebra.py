import numpy as np
from scipy.linalg import solve_triangular
import csv
from joblib import Parallel, delayed
import time
import math
import os
import json
#from memory_profiler import profile
#@profile

def add(a,b):
    c = {}
    nz_idx = set(a.keys()).union(set(b.keys())) 
    for i in nz_idx:
        c[i] = 0
        if i in a:
            c[i] += a[i]
        if i in b:
            c[i] += b[i]
    nz_idx.clear()
    return c

def innerprod(a,b):
    c = 0.0
    nz_idx = set(a.keys()).intersection(set(b.keys())) 
    for  i in nz_idx:
        c += a[i]*b[i]
    nz_idx.clear()
    return c
def norm(a):
    sqa=math.sqrt(innerprod(a,a))
    return sqa
def cleandict(a):
    result={}
    for x, y in a.items(): 
        if(abs(y) > 1e-7):result.update({x:y})
    return result
def get_Q_Hm(Hk,R,N):#set(a.keys()).intersection(set(b.keys()))
    m=len(R)
    print('\nm={}from len(R) in get_Q_Hm sub\n'.format(m))
    allI=''
    for p in range(N):allI=allI+'I'
    H0={allI:1.0}
    Hk['H0']=H0
    U=set()
    Hm={}
    listtrzH=[]
    listtrzzH=[]
    listtrxH=[]

    for i in range(m):
        listtrzH.append([])
        listtrzzH.append([])
        listtrxH.append([])
    
    for key in Hk.keys():Hm[key]={}

    for i in range(m):U=set.union(set(Hk['H'+str(i)].keys()),U)

    Rnp=np.zeros((m, m))

    for i in range(m):
        for Hkey,Hvalue in R['R'+str(i)].items():
            if (len(Hvalue.values())==0):Rnp[i][int(Hkey[1])]=0.0
            else:Rnp[i][int(Hkey[1])]=list(Hvalue.values())[0]
    Rnpt=np.transpose(Rnp)
    for row in U:
        Hknp=np.zeros((m, 1))

        for i in range(m):
            if row in Hk['H'+str(i)].keys():Hknp[i][0]=Hk['H'+str(i)][row] ## Hknp is not correct


        x = solve_triangular(Rnpt, Hknp, lower=True)
        for i in range(m):
            if(abs(x[i][0])>1e-7):
                Hm['H'+str(i)].update({row:x[i][0]})
                if(row.count('z')==1 and row.count('I')==N-1):listtrzH[i].append([x[i][0]])
                if(row.count('x')==1 and row.count('I')==N-1):listtrxH[i].append([x[i][0]])
                if(row.count('z')==2 and row.count('I')==N-2):listtrzzH[i].append([x[i][0]])
    ################################
    listtr={'listtrzH':listtrzH,'listtrzzH':listtrzzH,'listtrxH':listtrxH}                  
    with open('listtrjason.txt', 'w') as outfile:
        json.dump(listtr, outfile)
    listtr.clear()
###########load###############################
#    with open('listtrjason.txt') as json_file:
#        listtr = json.load(json_file)
#    listtrzH=listtr['listtrzH']
#    listtrzzH=listtr['listtrzzH']
#    listtrxH=listtr['listtrxH']
    ######################################
    for i in range(m):
        listtrzH[i]=np.array(listtrzH[i])
        listtrxH[i]=np.array(listtrxH[i])
        listtrzzH[i]=np.array(listtrzzH[i])
    listtr={'listtrzH':listtrzH,'listtrzzH':listtrzzH,'listtrxH':listtrxH}        
    return Hm, listtr             
def get_Cm(R,E):## R is a dictionaries R[Ri:Hj-1] = Rnp[i] [j]
    m=len(R)
#    print(m)
    Enp=np.zeros((m, 1))
#    print(m)
    Rnp=np.zeros((m, m))

    for i in range(m):
        Enp[i][0]=pow(E,i)
        for Hkey,Hvalue in R['R'+str(i)].items():
            if (len(Hvalue.values())==0):Rnp[i][int(Hkey[1])]=0.0
            else:Rnp[i][int(Hkey[1])]=list(Hvalue.values())[0]
#    print(Rnp)
#    print(Enp)
    Rnpt=np.transpose(Rnp)
    x = solve_triangular(Rnpt, Enp,lower=True)
 #   print(Rnp)
 #   print(Enp)
 #   print(Rnpt)
 #   print(x)
    return x       

##################################################################################

def QRgetR(Hk,N):
    m=len(Hk)
    Hm={}
    Hp={}
    Hp2={}
 #   v={}
    R={}
    R5={}
    R1={}
    R2={}
    R3={}
    R4={}
    R0={}
#    print(m)
    allI=''
    for p in range(N):allI=allI+'I'
    R0['H0']={allI:1.0}
    for l in range(m):
        if allI in Hk['H'+str(l+1)]:
            R0['H'+str(l+1)]={allI:Hk['H'+str(l+1)].get(allI)}
        Hk['H'+str(l+1)].pop(allI, None)
        Hp['H'+str(l+1)]= Hk['H'+str(l+1)]

    Hk.clear()## to save memory from MemoryError. But you must read back Hk before getHm function
    for k in range(m):
        if k==1:
            for l in range(m):Hp['H'+str(l+1)].pop(allI[0:N-1]+'x', None)# next matrix
        if k==2:
            for l in range(m):Hp['H'+str(l+1)].pop(allI[0:N-1]+'y', None)# next matrix
        if k==3:
            for l in range(m):Hp['H'+str(l+1)].pop(allI[0:N-1]+'z', None)# next matrix
        if k==4:
            for l in range(m):Hp['H'+str(l+1)].pop(allI[0:N-2]+'xI', None)# next matrix
        if k==5:
            for l in range(m):Hp['H'+str(l+1)].pop(allI[0:N-2]+'yI', None)# next matrix             
        alpha=-norm(Hp['H'+str(k+1)]) # assuming pivot is usually zero and hope no loss of significance
        # maybe we need gc.collect() here to release memory used in norm which is innerprod subroutine
        if k==0:alphavec={allI[0:N-1]+'x':alpha}
        if k==1:alphavec={allI[0:N-1]+'y':alpha}
        if k==2:alphavec={allI[0:N-1]+'z':alpha}
        if k==3:alphavec={allI[0:N-2]+'xI':alpha}
        if k==4:alphavec={allI[0:N-2]+'yI':alpha}
        if k==5:alphavec={allI[0:N-2]+'zI':alpha}
    
#        if (k+1)%4==1:alphavec={allI[0:N-1]+'x':alpha}
#        if (k+1)%4==2:alphavec={allI[0:N-1]+'y':alpha}
#        if (k+1)%4==3:alphavec={allI[0:N-1]+'z':alpha}
#        if (k+1)%4==0:alphavec={allI[0:N-2]+'x'+'I':alpha}## for now, work only for m <= 4
        #u=add(Hp['H'+str(k+1)],alphavec)
#        for Hkey,Hvalue in addvector.items():
# addkey=Hkey
# addvalue=Hvalue

#if addkey in thisdict:
# thisdict.update({addkey:(thisdict.get(addkey)+addvalue)})
#else:
# thisdict.update({addkey:addvalue})
        #############==  u=add(Hp['H'+str(k+1)],alphavec)########
#        u={}
#        for Hkey,Hvalue in alphavec.items():
#            addkey=Hkey
#            addvalue=Hvalue
#        if addkey in Hp['H'+str(k+1)]:
#            u.update({addkey:(Hp['H'+str(l+1)].get(addkey)+addvalue)})
#        else:
#            Hp['H'+str(k+1)].update({addkey:addvalue})
####################################################################
        u=add(Hp['H'+str(k+1)],alphavec)
        nu=norm(u)
        v={x:y/nu for x, y in u.items()}## v = u/norm(u)
        u.clear()
        for l in range(m-k):
            a=(-2.0)*innerprod(v,Hp['H'+str(l+1+k)])
            vp={x:y*a for x, y in v.items()}
            Hp['H'+str(l+1+k)]=add(Hp['H'+str(l+1+k)],vp)# continue for next manupilation
            Hp['H'+str(l+1+k)]=cleandict(Hp['H'+str(l+1+k)])
            #print('\n',k,l,Hp['H'+str(l+1+k)])
            if k==0:R1['H'+str(l+1+k)]={Hkey:Hvalue for Hkey,Hvalue in Hp['H'+str(l+1+k)].items() if Hkey == allI[0:N-1]+'x'} # save the row in R
            if k==1:R2['H'+str(l+1+k)]={Hkey:Hvalue for Hkey,Hvalue in Hp['H'+str(l+1+k)].items() if Hkey == allI[0:N-1]+'y'}
            if k==2:R3['H'+str(l+1+k)]={Hkey:Hvalue for Hkey,Hvalue in Hp['H'+str(l+1+k)].items() if Hkey == allI[0:N-1]+'z'}
            if k==3:R4['H'+str(l+1+k)]={Hkey:Hvalue for Hkey,Hvalue in Hp['H'+str(l+1+k)].items() if Hkey == allI[0:N-2]+'xI'}
            if k==4:R5['H'+str(l+1+k)]={Hkey:Hvalue for Hkey,Hvalue in Hp['H'+str(l+1+k)].items() if Hkey == allI[0:N-2]+'yI'}
#            if k==5:R5['H'+str(l+1+k)]={Hkey:Hvalue for Hkey,Hvalue in Hp['H'+str(l+1+k)].items() if Hkey == allI[0:N-2]+'zI'}
            vp.clear()
 #   print('\n\n',R0)
 #   print('\n\n',R1)
 #   print('\n\n',R2)
 #   print('\n\n',R3)
#    print('\n\n',R4)
    R['R0']=R0
    R['R1']=R1
    R['R2']=R2
    R['R3']=R3
    R['R4']=R4
    R['R5']=R5
  
    return R  
