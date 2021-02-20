import gc
import math
import cmath
def autoHvector(model,parameters,Nx,Ny): 
    N=Nx*Ny
    result1={}
    result2={}
    result3={}
    result4={}
    result5={}
    c={}
    allI=''
    allx=''
    for p in range(N):
        allI=allI+'I'
        allx=allx+'x'

    if(model=='1dIsing'):
        print("using 1dIsing")
        for i in range(N):##0 1 2 3 N-1
            if(i==0):
                c.update({'zz'+allI[i+2:]:parameters['J']})
                c.update({'x'+allI[i+1:]:parameters['hx']})
                c.update({'z'+allI[i+1:]:parameters['hz']})
            elif(i==(N-2)):
                c.update({allI[:i]+'zz':parameters['J']})
                c.update({allI[:i]+'x'+allI[i+1:]:parameters['hx']})
                c.update({allI[:i]+'z'+allI[i+1:]:parameters['hz']})
            elif(i==(N-1)):
                c.update({allI[:i]+'x':parameters['hx']})
                c.update({allI[:i]+'z':parameters['hz']})
            else:
                c.update({allI[:i]+'zz'+allI[i+2:]:parameters['J']})
                c.update({allI[:i]+'x'+allI[i+1:]:parameters['hx']})
                c.update({allI[:i]+'z'+allI[i+1:]:parameters['hz']})
        for x, y in c.items(): 
            if(abs(y) > 1e-7):result1.update({x:y})
        return result1

    elif(model=='2dIsing'):
        print("using 2dIsing")
        upperconerlist=[]
        for b in range(N-Nx,N,1):upperconerlist.append(b)# a list with sites of upper coner
        for i in range(N):##0 1 2 3 N-1
            
            if((i%Nx)==0):## left edge
                c.update({allI[:i]+'zz'+allI[i+2:]:parameters['J']})## right bond
                c.update({allI[:i]+'x'+allI[i+1:]:parameters['hx']})## on site
                c.update({allI[:i]+'z'+allI[i+1:]:parameters['hz']})## on site
                if(i not in upperconerlist):c.update({allI[:i]+'z'+allI[i+1:i+Nx]+'z'+allI[i+Nx+1:]:parameters['J']}) 
                ## up bond
            elif((i%Nx)==(Nx-2)):##right next edge for inter sites
                c.update({allI[:i]+'zz'+allI[i+2:]:parameters['J']})## right bond
                c.update({allI[:i]+'x'+allI[i+1:]:parameters['hx']})## on site
                c.update({allI[:i]+'z'+allI[i+1:]:parameters['hz']})## on site
                if(i not in upperconerlist):c.update({allI[:i]+'z'+allI[i+1:i+Nx]+'z'+allI[i+Nx+1:]:parameters['J']}) ## up bond
            elif((i%Nx)==(Nx-1)):##right edge for inter sites
              #  c.update({allI[:i]+'zz'+allI[i+2:]:parameters['J']})## right bond
                c.update({allI[:i]+'x'+allI[i+1:]:parameters['hx']})## on site
                c.update({allI[:i]+'z'+allI[i+1:]:parameters['hz']})## on site
                if(i not in upperconerlist):c.update({allI[:i]+'z'+allI[i+1:i+Nx]+'z'+allI[i+Nx+1:]:parameters['J']}) ## up bond
            else:
                c.update({allI[:i]+'zz'+allI[i+2:]:parameters['J']})
                c.update({allI[:i]+'x'+allI[i+1:]:parameters['hx']})
                c.update({allI[:i]+'z'+allI[i+1:]:parameters['hz']})
                if(i not in upperconerlist):c.update({allI[:i]+'z'+allI[i+1:i+Nx]+'z'+allI[i+Nx+1:]:parameters['J']})# up bond
            ##print('\ncurrent sites in 1-d is ',i)
            ##print('\n current vector',c)
        for x, y in c.items():
            if(abs(y) > 1e-7):result2.update({x:y})
        return result2
    elif(model=='2dIsing_pinMz'):
        print("using 2dIsing_pinMz")
        upperconerlist=[]
        connerlistforMz=set()
        for b in range(N-Nx,N,1):upperconerlist.append(b)# a list with sites of upper coner
        
        for b in range(N-Nx,N,1):connerlistforMz.add(b)#upper edge
        for b in range(0,N,Nx):connerlistforMz.add(b)#left edge
        for b in range(Nx-1,N,Nx):connerlistforMz.add(b)#right edge
        for b in range(0,Nx,1):connerlistforMz.add(b)#bottom edge 
#        print(connerlistforMz)
        for i in range(N):##0 1 2 3 N-1
            if(i in connerlistforMz):c.update({allI[:i]+'z'+allI[i+1:]:parameters['hz']}) 
            if((i%Nx)==0):## left edge
                c.update({allI[:i]+'zz'+allI[i+2:]:parameters['J']})## right bond
                c.update({allI[:i]+'x'+allI[i+1:]:parameters['hx']})## on site
#                c.update({allI[:i]+'z'+allI[i+1:]:parameters['hz']})## on site
                if(i not in upperconerlist):c.update({allI[:i]+'z'+allI[i+1:i+Nx]+'z'+allI[i+Nx+1:]:parameters['J']}) 
                ## up bond
            elif((i%Nx)==(Nx-2)):##right next edge for inter sites
                c.update({allI[:i]+'zz'+allI[i+2:]:parameters['J']})## right bond
                c.update({allI[:i]+'x'+allI[i+1:]:parameters['hx']})## on site
#                 c.update({allI[:i]+'z'+allI[i+1:]:parameters['hz']})## on site
                if(i not in upperconerlist):c.update({allI[:i]+'z'+allI[i+1:i+Nx]+'z'+allI[i+Nx+1:]:parameters['J']}) ## up bond
            elif((i%Nx)==(Nx-1)):##right edge for inter sites
              #  c.update({allI[:i]+'zz'+allI[i+2:]:parameters['J']})## right bond
                c.update({allI[:i]+'x'+allI[i+1:]:parameters['hx']})## on site
#                c.update({allI[:i]+'z'+allI[i+1:]:parameters['hz']})## on site
                if(i not in upperconerlist):c.update({allI[:i]+'z'+allI[i+1:i+Nx]+'z'+allI[i+Nx+1:]:parameters['J']}) ## up bond
            else:
                c.update({allI[:i]+'zz'+allI[i+2:]:parameters['J']})
                c.update({allI[:i]+'x'+allI[i+1:]:parameters['hx']})
#                 c.update({allI[:i]+'z'+allI[i+1:]:parameters['hz']})
                if(i not in upperconerlist):c.update({allI[:i]+'z'+allI[i+1:i+Nx]+'z'+allI[i+Nx+1:]:parameters['J']})# up bond
            ##print('\ncurrent sites in 1-d is ',i)
            ##print('\n current vector',c)
        for x, y in c.items():
            if(abs(y) > 1e-7):result4.update({x:y})
        return result4
    elif(model=='2dIsing_pinMzLR'):
        print("using 2dIsing_pinMzLR")
        upperconerlist=[]
        connerlistforMz=set()
        for b in range(N-Nx,N,1):upperconerlist.append(b)# a list with sites of upper coner
        ## just ignore upper and bottom edges
#        for b in range(N-Nx,N,1):connerlistforMz.add(b)#upper edge
        for b in range(0,N,Nx):connerlistforMz.add(b)#left edge
        for b in range(Nx-1,N,Nx):connerlistforMz.add(b)#right edge
#        for b in range(0,Nx,1):connerlistforMz.add(b)#bottom edge 
#        print(connerlistforMz)
        for i in range(N):##0 1 2 3 N-1
            if(i in connerlistforMz):c.update({allI[:i]+'z'+allI[i+1:]:parameters['hz']}) 
            if((i%Nx)==0):## left edge
                c.update({allI[:i]+'zz'+allI[i+2:]:parameters['J']})## right bond
                c.update({allI[:i]+'x'+allI[i+1:]:parameters['hx']})## on site
#                c.update({allI[:i]+'z'+allI[i+1:]:parameters['hz']})## on site
                if(i not in upperconerlist):c.update({allI[:i]+'z'+allI[i+1:i+Nx]+'z'+allI[i+Nx+1:]:parameters['J']}) 
                ## up bond
            elif((i%Nx)==(Nx-2)):##right next edge for inter sites
                c.update({allI[:i]+'zz'+allI[i+2:]:parameters['J']})## right bond
                c.update({allI[:i]+'x'+allI[i+1:]:parameters['hx']})## on site
#                 c.update({allI[:i]+'z'+allI[i+1:]:parameters['hz']})## on site
                if(i not in upperconerlist):c.update({allI[:i]+'z'+allI[i+1:i+Nx]+'z'+allI[i+Nx+1:]:parameters['J']}) ## up bond
            elif((i%Nx)==(Nx-1)):##right edge for inter sites
              #  c.update({allI[:i]+'zz'+allI[i+2:]:parameters['J']})## right bond
                c.update({allI[:i]+'x'+allI[i+1:]:parameters['hx']})## on site
#                c.update({allI[:i]+'z'+allI[i+1:]:parameters['hz']})## on site
                if(i not in upperconerlist):c.update({allI[:i]+'z'+allI[i+1:i+Nx]+'z'+allI[i+Nx+1:]:parameters['J']}) ## up bond
            else:
                c.update({allI[:i]+'zz'+allI[i+2:]:parameters['J']})
                c.update({allI[:i]+'x'+allI[i+1:]:parameters['hx']})
#                 c.update({allI[:i]+'z'+allI[i+1:]:parameters['hz']})
                if(i not in upperconerlist):c.update({allI[:i]+'z'+allI[i+1:i+Nx]+'z'+allI[i+Nx+1:]:parameters['J']})# up bond
            ##print('\ncurrent sites in 1-d is ',i)
            ##print('\n current vector',c)
        for x, y in c.items():
            if(abs(y) > 1e-7):result5.update({x:y})
        return result5
        
    elif(model=='2dIsing_projectout_S=-1'):
        print("using 2dIsing_projectout_S=-1")
        upperconerlist=[]
        for b in range(N-Nx,N,1):upperconerlist.append(b)# a list with sites of upper coner
        for i in range(N):##0 1 2 3 N-1

            if((i%Nx)==0):## left edge
                c.update({allI[:i]+'zz'+allI[i+2:]:0.5*parameters['J']})## right bond
                c.update({allx[:i]+'yy'+allx[i+2:]:(-0.5)*parameters['J']})## right bond
                c.update({allI[:i]+'x'+allI[i+1:]:0.5*parameters['hx']})## on site
                c.update({allx[:i]+'I'+allx[i+1:]:0.5*parameters['hx']})## on site
#                 c.update({allI[:i]+'z'+allI[i+1:]:0.5*parameters['hz']})## on site
#                 c.update({allx[:i]+'y'+allx[i+1:]:0.5*parameters['hz']})## on site
                if(i not in upperconerlist):
                    c.update({allI[:i]+'z'+allI[i+1:i+Nx]+'z'+allI[i+Nx+1:]:0.5*parameters['J']})
                    c.update({allx[:i]+'y'+allx[i+1:i+Nx]+'y'+allx[i+Nx+1:]:(-0.5)*parameters['J']}) 
                ## up bond
            elif((i%Nx)==(Nx-2)):##right next edge for inter sites
                c.update({allI[:i]+'zz'+allI[i+2:]:0.5*parameters['J']})## right bond
                c.update({allx[:i]+'yy'+allx[i+2:]:(-0.5)*parameters['J']})## right bond
                c.update({allI[:i]+'x'+allI[i+1:]:0.5*parameters['hx']})## on site
                c.update({allx[:i]+'I'+allx[i+1:]:0.5*parameters['hx']})## on site
#                 c.update({allI[:i]+'z'+allI[i+1:]:0.5*parameters['hz']})## on site
#                 c.update({allx[:i]+'y'+allx[i+1:]:0.5*parameters['hz']})## on site
                if(i not in upperconerlist):
                    c.update({allI[:i]+'z'+allI[i+1:i+Nx]+'z'+allI[i+Nx+1:]:0.5*parameters['J']}) ## up bond
                    c.update({allx[:i]+'y'+allx[i+1:i+Nx]+'y'+allx[i+Nx+1:]:(-0.5)*parameters['J']}) 
            elif((i%Nx)==(Nx-1)):##right edge for inter sites
              #  c.update({allI[:i]+'zz'+allI[i+2:]:parameters['J']})## right bond
                c.update({allI[:i]+'x'+allI[i+1:]:0.5*parameters['hx']})## on site
                c.update({allx[:i]+'I'+allx[i+1:]:0.5*parameters['hx']})## on site
#                 c.update({allI[:i]+'z'+allI[i+1:]:0.5*parameters['hz']})## on site
#                 c.update({allx[:i]+'y'+allx[i+1:]:0.5*parameters['hz']})## on site
                if(i not in upperconerlist):
                    c.update({allI[:i]+'z'+allI[i+1:i+Nx]+'z'+allI[i+Nx+1:]:0.5*parameters['J']}) ## up bond
                    c.update({allx[:i]+'y'+allx[i+1:i+Nx]+'y'+allx[i+Nx+1:]:(-0.5)*parameters['J']}) 
            else:
                c.update({allI[:i]+'zz'+allI[i+2:]:0.5*parameters['J']})
                c.update({allx[:i]+'yy'+allx[i+2:]:(-0.5)*parameters['J']})
                c.update({allI[:i]+'x'+allI[i+1:]:0.5*parameters['hx']})
                c.update({allx[:i]+'I'+allx[i+1:]:0.5*parameters['hx']})
#                 c.update({allI[:i]+'z'+allI[i+1:]:0.5*parameters['hz']})
#                 c.update({allx[:i]+'y'+allx[i+1:]:0.5*parameters['hz']})
                if(i not in upperconerlist):
                    c.update({allI[:i]+'z'+allI[i+1:i+Nx]+'z'+allI[i+Nx+1:]:0.5*parameters['J']})# up bond
                    c.update({allx[:i]+'y'+allx[i+1:i+Nx]+'y'+allx[i+Nx+1:]:(-0.5)*parameters['J']})# up bond
                    
            ##print('\ncurrent sites in 1-d is ',i)
            ##print('\n current vector',c)
        
        for x, y in c.items():
            if(abs(y) > 1e-7):result3.update({x:y})
        return result3
    else:
        print("no such model has been found and return empty vector")
        return {}
    #else: print("no such model and return empty vector")
    #return result
    
def newkey(i,j):
    sign=1.0
    k=''
    for q in range(len(i)):
        shoot=i[q]
        target=j[q]
        if(shoot == target): output='I'
        elif(shoot == 'I'): output=target
        elif(target == 'I'): output=shoot
        elif(shoot=='x' and target=='y'): output='z';sign=sign*(1.0j)
        elif(shoot=='y' and target=='x'): output='z';sign=sign*(-1.0j)
        elif(shoot=='z' and target=='x'): output='y';sign=sign*(1.0j)
        elif(shoot=='x' and target=='z'): output='y';sign=sign*(-1.0j)
        elif(shoot=='y' and target=='z'): output='x';sign=sign*(1.0j)
        elif(shoot=='z' and target=='y'): output='x';sign=sign*(-1.0j)
        k=k+output
    return k,sign


def krylov_power_vector(a,b,n):# Parallel(n_jobs=2)(delayed(sqrt)(i ** 2) for i in range(10))
    allI=''
    for p in range(n):allI=allI+'I'
    c={}
    result={}
    sumallI=0
    for i in a.keys():
        for j in b.keys():
            if (i == j):
                sumallI=sumallI+a[i]*b[i]
            else:
                k,sign=newkey(i,j)
                if k in c:
                    value=c[k]+a[i]*b[j]*sign
                    if(abs(value)>1e-7):c.update( { k:value} )
                    # and abs(value.imag)>1e-7
                    else:del c[k]
                else:
                    value=a[i]*b[j]*sign
                    if(abs(value)>1e-7):c.update( { k: value})
                    #and abs(value.imag)>1e-7
    
    c.update( {allI:sumallI} )

    for x, y in c.items(): 
        if(abs(y.imag) > 1e-7):
                print('warning:complex element in Ham')
                result.update({x:y.real})
        else:result.update({x:y.real}) 
    return result

