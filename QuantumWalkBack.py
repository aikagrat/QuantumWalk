import math
import numpy as np
import cmath

#define final state
k=3 #number of sites

final = np.zeros((2*k,1),dtype=complex) #2k number of total sites (up and down)
#final[0][0] = 1. #spin up 
#final[1][0] = 1. #spin down
final[2*k-1][0]=1.

final = np.ones((2*k,1),dtype=complex)/math.sqrt(2*k)
final[2*k-1][0] = -final[2*k-1][0]
#final[2*k-3][0] = -final[2*k-3][0]

initial = np.zeros((2*k,1),dtype=complex)
initial[0][0]=initial[1][0]=1./math.sqrt(2)

print (final)


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
    #vn = np.array([[final[6][0],final[7][0]]])

    matrixC = np.zeros((2*k,2*k),dtype=complex)
    r1 = cmath.polar(final[2][0])
    r2 = cmath.polar(final[3][0])
    x=r2[0]/r1[0]
    theta = math.atan(x)
    
    phase =np.angle( final[2][0]/final[3][0])
    print ("phase",phase)
    
    sum_values= ( phase - math.pi)/2. #sum ksi + z
    z=2. #freedom in z --> freedom in alpha
    ksi= sum_values - z

    C = np.array([[cmath.rect(1,ksi)*math.cos(theta),cmath.rect(1,z)*math.sin(theta)],[-cmath.rect(1,-z)*math.sin(theta),cmath.rect(1,-ksi)*math.cos(theta)]])
    #print (C_n)
    
    invC = np.linalg.inv(C)
    for i in range (0,2*k,2):
        matrixC[0+i][0+i] = invC[0][0]
        matrixC[1+i][1+i] = invC[1][1]
        matrixC[0+i][1+i] = invC[0][1]          
        matrixC[1+i][0+i] = invC[1][0]
        
    listC.append (C)    
    
    
    prevstate = np.dot(invS,final)
    #print ("1",  prevstate )
    prevstate = np.dot(matrixC,prevstate)
    #print ("2",  prevstate )
    print ("after step j=",j," ",prevstate)
    listSt.append (prevstate)
    final = prevstate
    
C1 = np.dot(final,initial.transpose())
listC.append (C1)

    



