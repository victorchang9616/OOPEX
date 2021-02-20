
import numpy as np
from scipy.linalg import solve_triangular
import csv
from joblib import Parallel, delayed
import time
import gc
import os
import json


#from os.path import dirname, basename, isfile, join
#import glob
#modules = glob.glob(join(dirname(__file__), "*.py"))
#__all__ = [ basename(f)[:-3] for f in modules if isfile(f) and not f.endswith('__init__.py')]
import sys 
  
#print("This is the name of the program:", sys.argv[0])
#print("Argument List:", str(sys.argv))

#from sys.argv[2] import * 
#for i in range(len(sys.argv-1)):from sys.argv[i+1] import *
from kryovdiction import *
from kryovmultiply import *
from kryovlinear_algebra import *
from kryovobservable import *
#pip install U memory_profiler
J=0.5;hx=-1.6;hz=0.0;Nx=10;Ny=5
#J=1;hx=-1.05;hz=0.5;Nx=10;Ny=1
print('\nparameters for Hamiltonian\n')
print('J' + repr(J) + 'hx' + repr(hx)+'hz' + repr(hz)+'Nx' + repr(Nx)+'Ny' + repr(Ny)+'\n')
N=Nx*Ny
writem=False
writeHk=True
writeHm=True
writronsiteMz=False
writeMz=False
writeMx=False
writeCzz17=False
writeoneforall=False
writenpobs=True
writeCzzmid=True
parameters={
   'J':J,
   'hz':hz,
    'hx':hx
    }
Hk={}
allI=''
for p in range(N):allI=allI+'I'
H0={allI:1.0}
start_time = time.time()
Hk['H1']=autoHvector('2dIsing_pinMzLR',parameters,Nx,Ny)
#print(H1)
print("autoHvector--- %s seconds ---" % (time.time() - start_time))
start_time = time.time()
Hk['H2']=krylov_power_vector(Hk['H1'],Hk['H1'],N)
#print(H2)
print("H2--- %s seconds ---" % (time.time() - start_time))
start_time = time.time()
Hk['H3']=krylov_power_vector(Hk['H1'],Hk['H2'],N)
print("H3--- %s seconds ---" % (time.time() - start_time))
start_time = time.time()
Hk['H4']=krylov_power_vector(Hk['H1'],Hk['H3'],N)
print("H4--- %s seconds ---" % (time.time() - start_time))
start_time = time.time()
Hk['H5']=krylov_power_vector(Hk['H1'],Hk['H4'],N)
print("H5--- %s seconds ---" % (time.time() - start_time))

if(writeHk):
    with open('Hkjason.txt', 'w') as outfile:
        json.dump(Hk, outfile)
#    with open('H0.csv', 'w', newline="") as csv_file:  
#        writer = csv.writer(csv_file)
#        for key, value in H0.items():writer.writerow([key, value])
#    with open('H1.csv', 'w', newline="") as csv_file:  
#        writer = csv.writer(csv_file)
#        for key, value in Hk['H1'].items():writer.writerow([key, value])
#    with open('H2.csv', 'w', newline="") as csv_file:  
#        writer = csv.writer(csv_file)
#        for key, value in Hk['H2'].items():writer.writerow([key, value])
#    with open('H3.csv', 'w', newline="") as csv_file:  
#        writer = csv.writer(csv_file)
#        for key, value in Hk['H3'].items():writer.writerow([key, value])
#    with open('H4.csv', 'w', newline="") as csv_file:  
#        writer = csv.writer(csv_file)
#        for key, value in Hk['H4'].items():writer.writerow([key, value])
#    for i in range(m):
#        with open('R'+str(i)+'.csv', 'w', newline="") as csv_file:
#            writer = csv.writer(csv_file)
#        for key, value in R['R'+str(i)].items():writer.writerow([key, value])

'''
with open('H1.csv', 'w', newline="") as csv_file:  
    writer = csv.writer(csv_file)
    for key, value in H1.items():writer.writerow([key, value])
with open('H2.csv', 'w', newline="") as csv_file:  
    writer = csv.writer(csv_file)
    for key, value in H2.items():writer.writerow([key, value])
with open('H3.csv', 'w', newline="") as csv_file:  
    writer = csv.writer(csv_file)
    for key, value in H3.items():writer.writerow([key, value])

with open('H1.csv') as csv_file:
    reader = csv.reader(csv_file)
    read = [(x,float(y)) for (x,y) in reader]
    H1 = dict(read)
with open('H2.csv') as csv_file:
    reader = csv.reader(csv_file)
    read = [(x,float(y)) for (x,y) in reader]
    H2 = dict(read)
with open('H3.csv') as csv_file:
    reader = csv.reader(csv_file)
    read = [(x,float(y)) for (x,y) in reader]
    H3 = dict(read)
'''
start_time = time.time()
gc.collect() 
R=QRgetR(Hk,N)
gc.collect() 
print("QRgetR--- %s seconds ---" % (time.time() - start_time))

m=len(R)
#for i in range(m):
#    print(Hm['H'+str(i)])
#    with open('R'+str(i)+'.csv', 'w', newline="") as csv_file:
#        writer = csv.writer(csv_file)
#        for key, value in R['R'+str(i)].items():writer.writerow([key, value])
with open('Rjason.txt', 'w') as outfile:
    json.dump(R, outfile)
        
        
        
start_time = time.time()
with open('Hkjason.txt') as json_file:
    Hk = json.load(json_file)
    
#with open('H0.csv') as csv_file:
#    reader = csv.reader(csv_file)
#    read = [(x,float(y)) for (x,y) in reader]
#    H0 = dict(read)
#with open('H1.csv') as csv_file:
#    reader = csv.reader(csv_file)
#    read = [(x,float(y)) for (x,y) in reader]
#    Hk['H1'] = dict(read)
#with open('H2.csv') as csv_file:
#    reader = csv.reader(csv_file)
#    read = [(x,float(y)) for (x,y) in reader]
#    Hk['H2'] = dict(read)
#with open('H3.csv') as csv_file:
#    reader = csv.reader(csv_file)
#    read = [(x,float(y)) for (x,y) in reader]
#    Hk['H3'] = dict(read)
#with open('H4.csv') as csv_file:
#    reader = csv.reader(csv_file)
#    read = [(x,float(y)) for (x,y) in reader]
#    Hk['H4'] = dict(read)
gc.collect()    
Hm,listtr=get_Q_Hm(Hk,R,N) 
# listtr={'listtrzH':listtrzH,'listtrzzH'....'listtrxH':...} where listtrzH=[nparray1..m] noted 1st to Hm0 no need

gc.collect()
Hk.clear() 
print("get_Q_Hm--- %s seconds ---" % (time.time() - start_time))
if(writeHm):
    start_time = time.time()
    with open('Hmjason.txt', 'w') as outfile:
        json.dump(Hm, outfile)
    m=len(Hm)
#    start_time = time.time()
#    for i in range(m):
#    print(Hm['H'+str(i)])
#        with open('H'+str(i)+'m.csv', 'w', newline="") as csv_file:
#            writer = csv.writer(csv_file)
#            for key, value in Hm['H'+str(i)].items():writer.writerow([key, value])
    print("Save Hm--- %s seconds ---" % (time.time() - start_time))
'''
with open('C3vsE.csv', 'w', newline="") as csv_file:
    writer = csv.writer(csv_file)
#    d={}
    for E in range(-1000,1000,10):
        Cm=get_Cm(R,E/10)
#        d={Cm[0][0]:E}
        writer.writerow([E,Cm[3][0]])
'''

#for E in range(-10,20):Cm[E*0.1]=get_Cm(R,E*10.0)
#for E in range(3):
start_time = time.time()
m=len(Hm)
print('\nm=')
print(m)
Cm={}
for E in range(-20,-1,1):Cm[E*0.1]=get_Cm(R,E*0.1*N)
#E=-0.0
#Cm[E*0.1]=get_Cm(R,E*0.1*N)
#with open('Cmofe.csv', 'w', newline="") as csv_file:
#    writer = csv.writer(csv_file)
#    for key, value in Cm.items():
#        list=value.tolist()
#        writer.writerow([key, list[0],list[1],list[2],list[3],list[4]])
#print('\nCm[e=0]={}\n'.format(Cm))
print("for 30 get_Cm--- %s seconds ---" % (time.time() - start_time))
 #   if E ==0:Cm[0.4]=get_Cm(R,0.4*N)
 #   if E ==1:Cm[-0.2]=get_Cm(R,-0.2*N)
 #   if E ==2:Cm[-0.8]=get_Cm(R,-0.8*N)
#Cm=get_Cm(R,-22.0)
print('\nwriting result in csv\n')
start_time = time.time()
if(writeMx):

  with open('EvsMx.csv', 'w', newline="") as csv_file:
      writer = csv.writer(csv_file)
      for E,C in Cm.items():
          if(writem):
            for maxm in range(m-1):writer.writerow([E, Mx(C,Hm,N,maxm+1)])
          else:
            writer.writerow([E, Mx(C,Hm,N,m-1)])

if(writeMz):
  with open('EvsMz.csv', 'w', newline="") as csv_file:
      writer = csv.writer(csv_file)
      for E,C in Cm.items():
        if(writem):
          for maxm in range(m-1):
              Mzper,Mzonsite=Mz(C,Hm,N,Nx,Ny,maxm+1)
              writer.writerow([E, Mzper])
              if(writronsiteMz):np.savetxt("MzonsiteE="+str(E)+".csv", Mzonsite, delimiter=",")
        else:
          Mzper,Mzonsite=Mz(C,Hm,N,Nx,Ny,m-1)
          writer.writerow([E, Mzper])
          if(writronsiteMz):np.savetxt("MzonsiteE="+str(E)+".csv", Mzonsite, delimiter=",")
print("for 30 Mx&Mz--- %s seconds ---" % (time.time() - start_time))
      

if(writeCzz17):      
  start_time = time.time()
  with open('EvsCzz.csv', 'w', newline="") as csv_file:
      writer = csv.writer(csv_file)
    
      for E,C in Cm.items():#writer.writerow([E, Czzdr(C,Hm,N,Nx,1)])
              if(writem):
                for maxm in range(m-1):
                  for i in range(7):
                   Czzc,Czz=Czzdr(C,Hm,N,Nx,i+1,maxm+1)
                   writer.writerow([E, Czzc,Czz])
              else:
                for i in range(7):
                  Czzc,Czz=Czzdr(C,Hm,N,Nx,i+1,m-1)
                  writer.writerow([E, Czzc,Czz])
  print("for 30 Czz(1)-(7)--- %s seconds ---" % (time.time() - start_time))
if(writeoneforall):      
  start_time = time.time()
  with open('EvsALL.csv', 'w', newline="") as csv_file:
      writer = csv.writer(csv_file)
    
      for E,C in Cm.items():#writer.writerow([E, Czzdr(C,Hm,N,Nx,1)])
              if(writem):
                for maxm in range(m-1):
                  for i in range(5):
                   sumMx,sumMz,sitesxy,sumCzzc,sumCzz=obsoneforall(C,Hm,N,Nx,maxm+1,i)
                   if(writronsiteMz and i==0):np.savetxt("MzonsiteE="+str(E)+".csv", sitesxy, delimiter=",")
                   writer.writerow([E,sumMx,sumMz,sumCzzc,sumCzz])
              else:
                for i in range(5):
                    sumMx,sumMz,sitesxy,sumCzzc,sumCzz=obsoneforall(C,Hm,N,Nx,m-1,i)
                    if(writronsiteMz and i==0):np.savetxt("MzonsiteE="+str(E)+".csv", sitesxy, delimiter=",")
                    writer.writerow([E,sumMx,sumMz,sumCzzc,sumCzz])
              gc.collect() 
  print("for 30 Czz(1)-(7)--- %s seconds ---" % (time.time() - start_time))
start_time = time.time()
if(writenpobs):
  start_time = time.time()
  with open('Evsnpobs.csv', 'w', newline="") as csv_file:
      writer = csv.writer(csv_file)
      for E,C in Cm.items():
          
          listsumMx,listsumMz,listsumCzzst =obsnp(C,listtr,N,m)
          for i in range(m):writer.writerow([E,listsumMx[i],listsumMz[i],listsumCzzst[i]])
          
  print("for  npobs--- %s seconds ---" % (time.time() - start_time))
if(writeCzzmid):      
  start_time = time.time()
  with open('EvsCzzmid.csv', 'w', newline="") as csv_file:
      writer = csv.writer(csv_file)
    
      for E,C in Cm.items():#writer.writerow([E, Czzdr(C,Hm,N,Nx,1)])
#              print("checking in for E,C in Cm.items()")
#              print(E)
              if(writem):
                   
                   for maxm in range(m-1):
 #                    print("writem=true")
                     Czzdrcnp,Czzdrnp=Czzdrmid(C,Hm,N,Nx,maxm+1)
                     for i in range(len(Czzdrcnp)):
                       writer.writerow([E,Czzdrcnp[i],Czzdrnp[i]])
  #                 if(writronsiteMz and i==0):np.savetxt("MzonsiteE="+str(E)+".csv", sitesxy, delimiter=",")
                   
              else:
                   Czzdrcnp,Czzdrnp=Czzdrmid(C,Hm,N,Nx,m-1)#Czzdrmid(Cm,Hm,N,Nx,maxm)
#                   print("writem=False")
                   for i in range(len(Czzdrcnp)):
                       writer.writerow([E,Czzdrcnp[i],Czzdrnp[i]])
              gc.collect() 
  print("for 30 Czzmid--- %s seconds ---" % (time.time() - start_time))
#print(Mx(Cm,Hm,N))
#print(Mx(Cm,Hm,N))
print('\nresult written in csv\n')
#print(Mz(Cm,Hm,N))

