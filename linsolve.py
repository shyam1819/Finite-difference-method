#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
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
    print(a)
    X=backsubstitution(a)
    return X
def backsubstitution(a):
    """" The Reduced Row echelon form is then back substituted to get the solution. this will do the back sunstitution and prints the answer"""
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
    print("The solution is: ",x)
        
#a=np.array([[6,12,1,8],[4,7,3,9],[2,2,-2,10]])
#a=np.array([[0,0,10,10],[11,0,1,20],[0,10,1,30]])    # [1.7272, 2.9, 1]
#a=np.array([[-3,2,-5,-14],[2,-3,4,10],[1,1,1,4]])     # [3, 0, 1] infinitely many solutions
#a=np.array([[1,3,-1,-4],[4,-1,2,3],[2,-1,-3,1]])           # [0.27272727 -1.32727272 0.29090909]
#a=np.array([[1,3,4],[3,-1,7]])
#a=np.array([[0,2,1,4],[1,1,2,6],[2,1,1,7]])
#a=np.array([[2,1,3],[4,6,9]])                             #[1.125 0.75]
print("This codes takes Augmented matrix as input [A|B] where AX=B")
n=int(input("Enter no of linear equations: "))
a=np.zeros((n,n+1))
a=a.astype(np.float)
for i in range(n):
    for j in range(n+1):
        a[i,j]=float(input("Enter value of a[{}][{}]: ".format(i,j)))
print("############################################")
print("The entered matrix is: ")
print(a)
print("Row echelon form of the given matrix is: ")
gausselimination(a)
print("############################################")


# In[ ]:





# In[ ]:





# In[ ]:




