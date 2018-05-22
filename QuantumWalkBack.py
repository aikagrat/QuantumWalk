import math
import numpy as np
import cmath

#define final state
k=3 #number of sites

final = np.zeros((2*k,1),dtype=complex) #2k number of total sites (up and down)
#final[0][0] = 1. #spin up 
#final[1][0] = 1. #spin down
#final[2*k-1][0]=1.

#final = np.ones((2*k,1),dtype=complex)/math.sqrt(2*k)
#final[2*k-1][0] = -final[2*k-1][0]
#final[3][0] = final[2*k-2][0] = 0. #| 1 down > , |3 up >
#final[2*k-3][0] = -final[2*k-3][0]

initial = np.zeros((2*k,1),dtype=complex)
initial[0][0]=initial[1][0]=1./math.sqrt(2)

u1 = np.ones((2,1),dtype=complex)/math.sqrt(2)
u2 = np.ones((2,1),dtype=complex)/math.sqrt(2)

final[0][0] = u1[0][0]
final[2][0]= u1[1][0]
final[2*k-3][0]= u2[0][0]
final[2*k-1][0]= -u2[1][0]

Final = final


print (Final)


#definition of invS
invS = np.zeros((2*k,2*k),dtype=complex)

for i in range (0,2*k,2) :
    invS[0+i][0+i] = 1.
    if (i+3)< 2*k :
        invS[1+i][3+i] = 1. #S-1
        #matrixS[3+i][1+i] = 1.
        


n=2 #number os steps
listSt = []
listC = []

listSt.append (final)

for j in range (n,1,-1) : 
    
    #definition of C
    v1 = np.array([[final[2][0],final[3][0]]])
    #v1=u1
    #vn = np.array([[final[2*k-3][0],final[2*k-1][0]]])

    matrixC = np.zeros((2*k,2*k),dtype=complex)
    r1 = cmath.polar(final[2][0])
    r2 = cmath.polar(final[3][0])
    x=r2[0]/r1[0]
    theta = math.atan(x)
    
    phase =np.angle( final[2][0]/final[3][0])
    #print ("phase",phase)
    
    sum_values= ( phase - math.pi)/2. #sum ksi + z
    z=0.5 #freedom in z --> freedom in alpha
    ksi= sum_values - z

    C = np.array([[cmath.exp(ksi*1j)*math.cos(theta),cmath.exp(z*1j)*math.sin(theta)],[-cmath.exp(-z*1j)*math.sin(theta),cmath.exp(-ksi*1j)*math.cos(theta)]])
    print (C)
    
    invC = np.linalg.inv(C)
    for i in range (0,2*k,2):
        matrixC[0+i][0+i] = invC[0][0]
        matrixC[1+i][1+i] = invC[1][1]
        matrixC[0+i][1+i] = invC[0][1]          
        matrixC[1+i][0+i] = invC[1][0]
        
    listC.append (matrixC)    
    
    
    prevstate = np.dot(invS,final)
    #print ("1",  prevstate )
    prevstate = np.dot(matrixC,prevstate)
    #print ("2",  prevstate )
    print ("after step j=",j," ",prevstate)
    listSt.append (prevstate)
    final = prevstate
    
C1 = np.dot(final,initial.transpose())
listC.append (C1)

    



