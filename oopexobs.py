
import numpy as np
from scipy.linalg import solve_triangular
import csv
from joblib import Parallel, delayed
import time
import gc
import os


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
###############input dim#############
Nx=10
Ny=10
########################
N=Nx*Ny
writem=False
writeHk=True
writeHm=True
writronsiteMz=False
writeMz=False
writeMx=False
writeCzz17=False
writeoneforall=False
writeCzzmid=True

Hm={}
allI=''
for p in range(N):allI=allI+'I'


        
        
        
start_time = time.time()
with open('H0m.csv') as csv_file:
    reader = csv.reader(csv_file)
    read = [(x,float(y)) for (x,y) in reader]
    Hm['H0'] = dict(read)
with open('H1m.csv') as csv_file:
    reader = csv.reader(csv_file)
    read = [(x,float(y)) for (x,y) in reader]
    Hm['H1'] = dict(read)
with open('H2m.csv') as csv_file:
    reader = csv.reader(csv_file)
    read = [(x,float(y)) for (x,y) in reader]
    Hm['H2'] = dict(read)
with open('H3m.csv') as csv_file:
    reader = csv.reader(csv_file)
    read = [(x,float(y)) for (x,y) in reader]
    Hm['H3'] = dict(read)
with open('H4m.csv') as csv_file:
    reader = csv.reader(csv_file)
    read = [(x,float(y)) for (x,y) in reader]
    Hm['H4'] = dict(read)
print("loading Hm--- %s seconds ---" % (time.time() - start_time))
m=len(Hm)
R={}
#############input R by hand###############################3
R0={'H0':{allI: 1.0},'H1':{},'H2':{allI: 3646.000000000008},'H3':{allI: -10800.0},'H4':{allI: 39183899.379999794}}
R1={'H1':{allI[1:N]+'x':  60.38211655780212},'H2':{allI[1:N]+'x':  -178.8609047790075},'H3':{allI[1:N]+'x':  648932.1940626296},'H4':{allI[1:N]+'x':  -6284651.5099670775}}
R2={'H2':{allI[1:N]+'y': 5085.134428583784},'H3':{allI[1:N]+'y': -44056.88065405637},'H4':{allI[1:N]+'y':107480748.12699133}}
R3={'H3':{allI[1:N]+'z': 517080.25208282424},'H4':{allI[1:N]+'z':  -8729347.961129703}}
R4={'H4':{allI[2:N]+'xI': 59831907.530641586}}

###############################################################
R['R0']=R0
R['R1']=R1
R['R2']=R2
R['R3']=R3
R['R4']=R4 
'''
for i in range(m):
#    print(Hm['H'+str(i)])
    with open('R'+str(i)+'.csv') as csv_file:
        reader = csv.reader(csv_file)
        read = [(y[3:3+N-1],y[5+N:]) for (x,y) in reader]  
        #R['R'+str(i)] = 
        print(dict(read))
'''
gc.collect()

Cm={}
for E in range(-19,11,2):
  Cm[E*0.1]=get_Cm(R,E*0.1*N)
#Et=-5.0
#Cm[Et*0.1]=get_Cm(R,Et*0.1*N)
gc.collect()
 #   if E ==0:Cm[0.4]=get_Cm(R,0.4*N)
 #   if E ==1:Cm[-0.2]=get_Cm(R,-0.2*N)
 #   if E ==2:Cm[-0.8]=get_Cm(R,-0.8*N)
#Cm=get_Cm(R,-22.0)
print('\nwriting result in csv\n')
if(writeMx):
  start_time = time.time()
  with open('EvsMx.csv', 'w', newline="") as csv_file:
      writer = csv.writer(csv_file)
      for E,C in Cm.items():
          if(writem):
            for maxm in range(m-1):writer.writerow([E, Mx(C,Hm,N,maxm+1)])
          else:
            writer.writerow([E, Mx(C,Hm,N,m-1)])
  print("for 30 Mx--- %s seconds ---" % (time.time() - start_time))
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
  with open('EvsALLobs.csv', 'w', newline="") as csv_file:
      writer = csv.writer(csv_file)
    
      for E,C in Cm.items():#writer.writerow([E, Czzdr(C,Hm,N,Nx,1)])
              if(writem):
                for maxm in range(m-1):
                  for i in range(5):
                   sumMx,sumMz,sitesxy,sumCzzc,sumCzz=obsoneforall(C,Hm,N,Nx,maxm+1,i)
                   if(writronsiteMz and i==0):np.savetxt("MzonsiteE="+str(round(E,1))+"obs.csv", sitesxy, delimiter=",")
                   writer.writerow([E,sumMx,sumMz,sumCzzc,sumCzz])
              else:
                for i in range(5):
                    sumMx,sumMz,sitesxy,sumCzzc,sumCzz=obsoneforall(C,Hm,N,Nx,m-1,i)
                    if(writronsiteMz and i==0):np.savetxt("MzonsiteE="+str(round(E,1))+"obs.csv", sitesxy, delimiter=",")
                    writer.writerow([E,sumMx,sumMz,sumCzzc,sumCzz])
              gc.collect() 
  print("for 30 Czz(1)-(7)--- %s seconds ---" % (time.time() - start_time))
if(writeCzzmid):      
  start_time = time.time()
  with open('EvsCzzmide-1.0to0ex-0.5.csv', 'w', newline="") as csv_file:
      writer = csv.writer(csv_file)
    
      for E,C in Cm.items():#writer.writerow([E, Czzdr(C,Hm,N,Nx,1)])
              if(writem):
                   
                   for maxm in range(m-1):
                     Czzdrcnp,Czzdrnp=Czzdrmid(C,Hm,N,Nx,maxm+1)
                     for i in range(len(Czzdrcnp)):
                       writer.writerow([E,Czzdrcnp[i],Czzdrnp[i]])
  #                 if(writronsiteMz and i==0):np.savetxt("MzonsiteE="+str(E)+".csv", sitesxy, delimiter=",")
                   
              else:
                   Czzdrcnp,Czzdrnp=Czzdrmid(C,Hm,N,Nx,m-1)
                   for i in range(len(Czzdrcnp)):
                       writer.writerow([E,Czzdrcnp[i],Czzdrnp[i]])
              gc.collect() 
  print("for 30 Czzmid--- %s seconds ---" % (time.time() - start_time))
#print(Mx(Cm,Hm,N))
print('\nresult written in csv\n')
print(R)
#print(Mz(Cm,Hm,N))

