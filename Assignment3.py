# -*- coding: utf-8 -*-

'''
Assignment 3. 


Instructions:
    
 A simple polygon is closed polygonal chain of line segments that do not 
cross each other. That is, it consists of finitely many line segments, each 
line segment endpoint is shared by two segments, and the segments do not 
otherwise intersect. In other words, a simple polygon is a polygon whose 
sides do not cross.

1. Complete all the functions given at the end of the program. 

2. Initialize the first two variables in your program with your name and ID

3. Name the file as shown in the following example 
   Assuming your ID is 2020ABA40001
   BITS_2020ABA40001.py
   
4. Upload the file,

5. Avoid Plagiarism. This may have serious consequences. 

Deadline 25th of April 5:00 PM. 


'''

STUDENT_NAME='SOMYA UPADHYAY'
STUDENT_ID='2017B5A40962P'

import numpy as np
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon
import matplotlib.pyplot as plt


fig, ax = plt.subplots(figsize=(6,6))
VL=np.array([[0.2,0.5],[1,0.2],[2,.4],[1.5,3],[0,2]])

RL=np.random.rand(10,2)

polygon = Polygon(VL,True)
polygon2 = Polygon(RL,True)
patches=[]
patches.append(polygon)

p=PatchCollection(patches,color='green',alpha=0.5)
ax.add_collection(p)
plt.xlim([np.min(VL[:,0]),np.max(VL[:,0])])
plt.ylim([np.min(VL[:,1]),np.max(VL[:,1])])
plt.show()


########################################################################

def check_simplepolygon(VL):
    '''
    Parameters
    ----------
    VL : An sequence of vertices in the form of a 2D array

    Returns
    -------
    True/False depending on whether VL is a simple polygon or not.

    '''

    x=[]
    y=[]
    for i,j in VL:
        x.append (i)
        y.append(j)
    l=len(VL)
    for ii in range (0,l-1):
        for jj in range (ii+1,l):
            i=ii%(l-1)
            j=jj%(l-1)
            X=[[x[i+1]-x[i] ,x[j]-x[j+1]], 
               [y[i+1]-y[i], y[j]-y[j+1]]]
            C=[[x[j]-x[i]               ], [y[j]-y[i]               ]]
            Xinv=np.linalg.inv(X)
            Xnorm=np.linalg.norm(X)
            mu=1/Xnorm*np.dot(Xinv,C)
            mu.reshape([2])
            if ((mu[0]<1 and mu[0]>0) and (mu[1]>0 and mu[1]<1)):
                return False
    return True

def check_convexity(VL):
    '''
    Parameters
    ----------
    VL : An sequence of vertices in the form of a 2D array

    Returns
    -------
    True/False depending on whether VL forms a boundary of a convexy polygon.

    '''
    l=len(VL)
    a=VL[0][0]*VL[1][1]-VL[1][0]*VL[0][1]
    for i in range (1,l-1):
        a=a*(VL[i][0]*VL[i+1][1]-VL[i+1][0]*VL[i][1])
        if a<=0:
            return False
    return True

def point_membership(P,VL):
    '''
    Parameters
    ----------
    P : a 2D point example, P = np.array([1,2])
    VL : An sequence of vertices in the form of a 2D array

    Returns
    -------
    Should an integer type 
    1 if the P is inside the boundaries defined by VL
    0 if the P is outside the boundaries defined by VL
    -1 if the P is on the boundary defined by VL

    '''
     def cro(V1,V2):
        v1=(V1[0]-V1[1])
        v2=(V2[0]-V2[1])
        z= v1[0]*v2[1] - v2[0]*v1[1]
        return z
    O=np.array([0,0])#one point outside
    x=[]
    y=[]
    count=0 #counts intersection
    
    op=-2
    for i,j in VL:
        x.append (i)
        y.append(j)
    l=len(VL)
    for i in range (0,l):
        X=[[x[(i+1)%l]-x[i] ,
            -P[0]+O[0]], 
        [y[(i+1)%l]-y[i], 
         -P[1]+O[1]]]
        C=np.array([O[0]-x[i] , O[1]-y[i] ])
        Xnorm=np.linalg.det(X)
        if (Xnorm==0):
            if(P[0]==O[0] and P[1]==O[1]): 
                op= 0
            else:
                m=(P[0]-x[i])/(x[(i+1)%l]-x[i])
                if((m<=1 and m>=0)):
                    return -1
            continue              
        
        Xinv=np.linalg.inv(X)
        #print('X\'',Xinv*Xnorm)
        mu=np.dot(Xinv,C)
        mu.reshape([2])
        I=np.array([O[0]*(1-mu[1])+P[0]*mu[1],O[1]*(1-mu[1])+P[1]*mu[1]]).T
        #print('I=',I,'mu1=',mu[1],'mu0=',mu[0])
        if ((mu[0]<1 and mu[0]>0) and (mu[1]>0 and mu[1]<1)):#counting intersection lines
            count=count+1
        if ((mu[1]==1 or mu[1]==0) and (mu[0]<=1 and mu[0]>=0)):#boundary  case
            return -1
        if ((mu[0]==1 or mu[0]==0)):# and (mu[1]>=0 and mu[1]<=1) ):#intersecting a vertex
            L=cro([P,I],[np.array([x[i],y[i]]),I])*cro([P,I],[np.array(x[(i+1)%l],y[(i+1)%l]),I])
            #print('L=',L)
            if (L < 0):
                count=count+1
            if (I[0]==P[0] and I[1]==P[1]):
                op= -1
    if (count%2==1 and op!=-1):
        op= 1
    if (count%2==0 and op!=-1):
        op= 0
    #print( count)
    return op

def find_intersection(VL1,VL2):
    '''
    Parameters
    ----------
    VL1 : a 2D array, shape=(N rows, 2 columns)
        A sequence of vertices, which form the boundary of a solid 1
    VL2 : a 2D array, shape=(N rows, 2 columns)
        A sequence of vertices, which form the boundary of a solid 2

    Returns
    -------
    VL_int : A sequence of vertices of the boundary of the intersection 
     of the two solids 1 and 2.

    '''
    def intersept(A1,A2,B1,B2):
        m1=(A1[1]-A2[1])/(A1[0]-A2[0])
        m2=(B1[1]-B2[1])/(B1[0]-B2[0])
        X=[ [A2[0]-A1[0] ,  -B2[0]+B1[0]], 
            [A2[1]-A1[1] ,  -B2[1]+B1[1]] ]
        C=np.array([-A1[0]+B1[0] , -A1[1]+B1[1] ])
        Xnorm=np.linalg.det(X)
        if (Xnorm==0):  
            if(m1==m2):
                if((B1[1]-A1[1])*(A2[0]-A1[0])<(A2[1]-A1[1])*(A2[0]-A1[0])):
                    return B1
                elif((B2[1]-A1[1])*(A2[0]-A1[0])<(A2[1]-A1[1])*(A2[0]-A1[0])):
                    return B2
            else:
                return [-1]    
        
        Xinv=np.linalg.inv(X)
        #print('X\'',Xinv*Xnorm)
        mu=np.dot(Xinv,C)
        mu.reshape([2])
        I=[B1[0]*(1-mu[1])+B2[0]*mu[1],B1[1]*(1-mu[1])+B2[1]*mu[1]]
        if(mu[0]> 0 and mu[0] < 0.9999999999999999 and mu[1] > 0 and mu[1] < 0.9999999999999999 ):
            return I
        else:
            return[-1]
    
        
    
    V1=[]
    V2=[]
    VL_int=[]
    l1=len(VL1)
    l2=len(VL2)
    for i in range(l1):
        ii=(i+1)%l1
        V1.append(VL1[i])
        V3=[]
        for j in range (l2):
            jj=(j+1)%l2
            I=intersept(VL1[i],VL1[ii],VL2[j],VL2[jj])
            if (len(I)==2):
                V3.insert(0, I)
        if (len(V3)>0):
            for k in range (len(V3)):
                V1.append(np.array(V3[k]))
   
    for i in range(l2):
        ii=(i+1)%l2
        V2.append(VL2[i])
        V3=[]
        for j in range (l1):
            jj=(j+1)%l1
            I=intersept(VL2[i],VL2[ii],VL1[j],VL1[jj])
            if (len(I)==2):
                V3.insert(0,I)
        if (len(V3)>0):
            for k in range (len(V3)):
                V2.append(np.array(V3[k]))
    
    for i in range(len(V2)):
        if(point_membership(V2[i],V1)!=0):
            VL_int.append(V2[i])
    if(len(VL_int)==0):
        for i in range(len(V1)):
            if(point_membership(V1[i],V2)!=0):
                VL_int.append(V1[i])
    return np.array(VL_int)
    
def find_union(VL1,VL2):
    '''
    Parameters
    ----------
    VL1 : a 2D array, shape=(N rows, 2 columns)
        A sequence of vertices, which form the boundary of a solid 1
    VL2 : a 2D array, shape=(N rows, 2 columns)
        A sequence of vertices, which form the boundary of a solid 2

    Returns
    -------
    VL_int : A sequence of vertices of the boundary of the union 
     of the two solids 1 and 2.

    '''
    def intersept(A1,A2,B1,B2):
        m1=(A1[1]-A2[1])/(A1[0]-A2[0])
        m2=(B1[1]-B2[1])/(B1[0]-B2[0])
        X=[ [A2[0]-A1[0] ,  -B2[0]+B1[0]], 
            [A2[1]-A1[1] ,  -B2[1]+B1[1]] ]
        C=np.array([-A1[0]+B1[0] , -A1[1]+B1[1] ])
        Xnorm=np.linalg.det(X)
        if (Xnorm==0):  
            if(m1==m2):
                if((B1[1]-A1[1])*(A2[0]-A1[0])<(A2[1]-A1[1])*(A2[0]-A1[0])):
                    return B1
                elif((B2[1]-A1[1])*(A2[0]-A1[0])<(A2[1]-A1[1])*(A2[0]-A1[0])):
                    return B2
            else:
                return [-1]    
        
        Xinv=np.linalg.inv(X)
        #print('X\'',Xinv*Xnorm)
        mu=np.dot(Xinv,C)
        mu.reshape([2])
        I=[B1[0]*(1-mu[1])+B2[0]*mu[1],B1[1]*(1-mu[1])+B2[1]*mu[1]]
        if(mu[0]> 0 and mu[0] < 0.9999999999999999 and mu[1] > 0 and mu[1] < 0.9999999999999999 ):
            return I
        else:
            return[-1]
    
    def compare(a,b):
        ep=0.00000001
        if(abs(a[0]-b[0])<=ep and abs(a[1]-b[1])<=ep):
            return True
        else:
            return False
    
    def cw(a,b,c):
        x1=-(b[0]-a[0])
        y1=-(b[1]-a[1])
        x2=(c[0]-a[0])
        y2=(c[1]-a[1])
        return(x1*y2-x2*y1)
        
    V1=[]
    V2=[]
    VL_int=[]
    l1=len(VL1)
    l2=len(VL2)    
    for i in range(l1):
        ii=(i+1)%l1
        V1.append(VL1[i])
        V3=[]
        for j in range (l2):
            jj=(j+1)%l2
            I=intersept(VL1[i],VL1[ii],VL2[j],VL2[jj])
            if (len(I)==2):
                V3.insert(0, I)
        if (len(V3)>0):
            for k in range (len(V3)):
                V1.append(np.array(V3[k]))
    for i in range(l2):
        ii=(i+1)%l2
        V2.append(VL2[i])
        V3=[]
        for j in range (l1):
            jj=(j+1)%l1
            I=intersept(VL2[i],VL2[ii],VL1[j],VL1[jj])
            if (len(I)==2):
                V3.insert(0,I)
        if (len(V3)>0):
            for k in range (len(V3)):
                V2.append(np.array(V3[k]))
    ccw1=0
    ccw2=0
    for i in range (l1-2):
        ccw1=ccw1 + cw(V1[i+1],V1[0],V1[i+2])
    for i in range (l2-2):
        ccw2=ccw2 + cw(V2[i+1],V2[0],V2[i+2])                   
    if (ccw1<0):
        V1.reverse()
    if (ccw2>0):
        V2.reverse()  
    
    ###segmengts in V1 and V2###
    VL_i=[]
    VL_1=[]
    i=-1
    while(i+1<(len(V1))):
        i=i+1
        ip=(i+1)%len(V1)
        while(i<len(V1)and point_membership((V1[i%len(V1)]+V1[ip])/2,V2)==0):
            ip=(i+1)%len(V1)
            VL_i.append(V1[i])
            i=i+1
        VL_i.append(V1[i%len(V1)])
        if (len(VL_i)>2):
            VL_1.append(VL_i)
        VL_i=[]
    VL_2 =[]
    VL_i =[]
    i=-1
    while(i+1<(len(V2))):
        i=i+1
        ip=(i+1)%len(V2)
        while(i<len(V2)and point_membership((V2[i%len(V2)]+V2[ip])/2,V1)==0):
            ip=(i+1)%len(V2)
            VL_i.append(V2[i])
            i=i+1
        VL_i.append(V2[i%len(V2)])
        if (len(VL_i)>2):
            VL_2.append(VL_i)
        VL_i=[]
    
    i=0
    j=0
    #print('VL_1=',np.array(VL_1))
    #print('VL_2=',np.array(VL_2))
    ##join segments###
    while(i<len(VL_1)):
        for k in VL_1[i]:
            VL_int.append(k)    
        for k in VL_2:
            if ( compare(k[0],VL_1[i][-1])):
                for j in k[1:]:
                     VL_int.append(j)
                
            else:
                k.reverse()
                if( compare(k[0],VL_1[i][-1])):
                    for j in k[1:]:
                        VL_int.append(j)
        i=i+1       
    return np.array(VL_int)
    
def find_difference(VL1, VL2):
    '''
    Parameters
    ----------
    VL1 : a 2D array, shape=(N rows, 2 columns)
        A sequence of vertices, which form the boundary of a solid 1
    VL2 : a 2D array, shape=(N rows, 2 columns)
        A sequence of vertices, which form the boundary of a solid 2

    Returns
    -------
    VL_int : A sequence of vertices of the boundary of the difference  
     of the two solids 1 and 2.
     S1-S2.
    '''
    def intersept(A1,A2,B1,B2):
        m1=(A1[1]-A2[1])/(A1[0]-A2[0])
        m2=(B1[1]-B2[1])/(B1[0]-B2[0])
        X=[ [A2[0]-A1[0] ,  -B2[0]+B1[0]], 
            [A2[1]-A1[1] ,  -B2[1]+B1[1]] ]
        C=np.array([-A1[0]+B1[0] , -A1[1]+B1[1] ])
        Xnorm=np.linalg.det(X)
        if (Xnorm==0):  
            if(m1==m2):
                if((B1[1]-A1[1])*(A2[0]-A1[0])<(A2[1]-A1[1])*(A2[0]-A1[0])):
                    return B1
                elif((B2[1]-A1[1])*(A2[0]-A1[0])<(A2[1]-A1[1])*(A2[0]-A1[0])):
                    return B2
            else:
                return [-1]    
        
        Xinv=np.linalg.inv(X)
        #print('X\'',Xinv*Xnorm)
        mu=np.dot(Xinv,C)
        mu.reshape([2])
        I=[B1[0]*(1-mu[1])+B2[0]*mu[1],B1[1]*(1-mu[1])+B2[1]*mu[1]]
        if(mu[0]> 0 and mu[0] < 0.9999999999999999 and mu[1] > 0 and mu[1] < 0.9999999999999999 ):
            return I
        else:
            return[-1]
    
    def compare(a,b):
        ep=0.00000001
        if(abs(a[0]-b[0])<=ep and abs(a[1]-b[1])<=ep):
            return True
        else:
            return False
    
    def cw(a,b,c):
        x1=-(b[0]-a[0])
        y1=-(b[1]-a[1])
        x2=(c[0]-a[0])
        y2=(c[1]-a[1])
        return(x1*y2-x2*y1)
        
    V1=[]
    V2=[]
    VL_int=[]
    l1=len(VL1)
    l2=len(VL2)    
    for i in range(l1):
        ii=(i+1)%l1
        V1.append(VL1[i])
        V3=[]
        for j in range (l2):
            jj=(j+1)%l2
            I=intersept(VL1[i],VL1[ii],VL2[j],VL2[jj])
            if (len(I)==2):
                V3.insert(0, I)
        if (len(V3)>0):
            for k in range (len(V3)):
                V1.append(np.array(V3[k]))
    for i in range(l2):
        ii=(i+1)%l2
        V2.append(VL2[i])
        V3=[]
        for j in range (l1):
            jj=(j+1)%l1
            I=intersept(VL2[i],VL2[ii],VL1[j],VL1[jj])
            if (len(I)==2):
                V3.insert(0,I)
        if (len(V3)>0):
            for k in range (len(V3)):
                V2.append(np.array(V3[k]))
    ccw1=0
    ccw2=0
    for i in range (l1-2):
        ccw1=ccw1 + cw(V1[i+1],V1[0],V1[i+2])
    for i in range (l2-2):
        ccw2=ccw2 + cw(V2[i+1],V2[0],V2[i+2])                   
    if (ccw1<0):
        V1.reverse()
    if (ccw2>0):
        V2.reverse()  
    
    ##########################################comparing V1 and V2#####################################
    VL_i=[]
    VL_1=[]
    i=-1
    while(i+1<(len(V1))):
        i=i+1
        ip=(i+1)%len(V1)
        while(i<len(V1)and point_membership((V1[i%len(V1)]+V1[ip])/2,V2)==0):
            ip=(i+1)%len(V1)
            VL_i.append(V1[i])
            i=i+1
        VL_i.append(V1[i%len(V1)])
        if (len(VL_i)>2):
            VL_1.append(VL_i)
        VL_i=[]
    if (point_membership(VL_1[0][0],V2)==0):
        for j in VL_1[0][1:]:
            VL_1[-1].append(j)
        VL_1=VL_1[1:]
    VL_2 =[]
    VL_i =[]
    i=-1
    while(i+1<(len(V2))):
        i=i+1
        ip=(i+1)%len(V2)
        while(i<len(V2)and point_membership((V2[i%len(V2)]+V2[ip])/2,V1)==1):
            ip=(i+1)%len(V2)
            VL_i.append(V2[i])
            i=i+1
        VL_i.append(V2[i%len(V2)])
        if (len(VL_i)>1):
            VL_2.append(VL_i)
        VL_i=[]
    if (point_membership(VL_2[0][0],V1)==1):
        for j in VL_2[0][1:]:
            VL_2[-1].append(j)
        VL_2=VL_2[1:]
    print('VL_1=',np.array(VL_1))
    print('VL_2=',np.array(VL_2))
    i=-1
    j=0
   ############################################## merge ################################################
    while(i<len(VL_1)):
        for k in VL_1[i]:
            VL_int.append(k)
        
        for k in VL_2:
            if ( compare(k[0],VL_1[i][-1]) and compare(k[-1],VL_1[i][0])):
                for j in k[1:]:
                     VL_int.append(j)
                
            else:
                k.reverse()
                if( compare(k[0],VL_1[i][-1]) and compare(k[-1],VL_1[i][0])):
                    for j in k[1:]:
                        VL_int.append(j)
                        
        i=i+1       
    return np.array(VL_int)
    
