#!/usr/bin/env python
# coding: utf-8

# In[3]:


import numpy as np
import matplotlib.pyplot as plt
import math
def rowexchange(a):
    """"Checks for Non zero elements along the diagonal does row exchanges accordingly"""
    i=1
    while a[0][0]==0:
        if i>=a.shape[0]:
            print("Matrix is singular")
            break
        c=np.copy(a[0])
        d=np.copy(a[i])
        a[0]=d
        a[i]=c
        i+=1
    d=elemination(a)
    return d
def elemination(a):
    """ Once rows are exchanged, this function produces an upper triangular matrix"""
    for i in range(1,a.shape[0]):
        b=a[i][0]/a[0][0]
        c=a[i]-a[0]*b
        a[i]=c
    return a
def gausselimination(a):
    """ This function calls the above functions and runs a loop to get row echelon form"""
    for i in range(a.shape[0]):
        b=np.copy(a[i::,i::])
        a[i::,i::]=rowexchange(b)
    x=backsubstitution(a)
    return x
def backsubstitution(a):
    """" The Reduced Row echelon form is then back substituted to get the solution."""
    b=a[:,-1]
    n = b.size
    x = np.zeros_like(b)
    if a[n-1, n-1] == 0:
        raise ValueError
    x[n-1] = b[n-1]/a[n-1, n-1]
    C = np.zeros((n,n))
    for i in range(n-2, -1, -1):
        bb = 0
        for j in range (i+1, n):
            bb += a[i, j]*x[j]
        C[i, i] = b[i] - bb
        x[i] = C[i, i]/a[i, i]
    return x
###############################################################################################################
a=float(input("Enter coefficient a: "))
b=float(input("Enter coefficient b: "))
c=float(input("Enter coefficient c: "))
n=int(input("Enter grid size: "))
lb=float(input("Enter lower boundary condition: "))
ub=float(input("Enter upper boundary condition: "))
dx=16.5/n
C=np.arange(0,16.50000001,dx)
###############################################################################################################
def fdmsolve(n,a,b,c,lb,ub):
    dx=16.5/n
    #C=np.array([0, 1.5, 3, 4.5, 6, 7.5, 9, 10.5, 12, 13.5, 15, 16.5])
    C=np.arange(0,16.50000001,dx)
    A=np.identity(len(C))
    for i in range(1,n):
        A[i,i-1]=a/(dx**2)
        A[i,i]=(c*(dx**2)-b*dx-2*a)/(dx**2)
        A[i,i+1]=(a+b*dx)/(dx**2)
    #print(A)
    B=np.sin(C*2)
    B[0]=lb
    B[-1]=ub
    A=np.concatenate((A,B[:,None]),axis=1)
    A=A.astype(np.float)
    X=gausselimination(A)
    return X
def analytical_sol(n):
    dx=16.5/n
    C=np.arange(0,16.50000001,dx)
    G=(5*C/66)-(np.sin(2*C)/4)
    return G
def rmserror(X,n):
    G=analytical_sol(n)
    rms=math.sqrt(np.sum(np.square(G-X))/(n+1))
    return rms
G=analytical_sol(n)
X=fdmsolve(n,a,b,c,lb,ub)
plt.subplot(211)
plt.plot(C,X)
plt.plot(C,G)
plt.xlabel("Values of the nodes\n")
plt.ylabel("Function value")
plt.legend(["FDM Solution","Analytical Solution"])
plt.show()
print("RMS error is {}".format(rmserror(X,n)))
rmse=[]
N=np.arange(10,310,10)
for i in N:
    rmse.append(rmserror(fdmsolve(i,a,b,c,lb,ub),i))
plt.subplot(212)
plt.plot(N,rmse)
plt.xlabel("No of grid spacings")
plt.ylabel("RMS error")
plt.title("RMS error as function of no of grid spaces")
plt.show()



    



# In[ ]:





# In[ ]:





# In[ ]:




