================================== OPPGAVE 2 ===================================
################################################################################
import sys
import os
import numpy as np
import numpy.linalg as nl
import scipy.linalg as sl
import scipy.sparse as sp
import matplotlib.pyplot as plt
from scipy.sparse import diags
################################################################################


def Eigenvalue_Program(program,algorithm,matrix):
    return 0

def Matrix_Generator(matrix,n):
    A = np.zeros((n, n)) #lager 0matrise
    forskjell = -int(len(matrix)/2)
    for i in matrix:
        A += i*np.eye(n,n,forskjell)
        forskjell+=1
    return A


def Run_Simulation(matrix,algorithm):   	# tar inn et tall matrix og en string algorithm
    if(algorithm == "Power_Eig"):       
        powerfile = open("Power_1.txt","w") 	#åpner fil for å legge til data
        for i in range(10,matrix):
            x = [1/matrix**(1/2)]*matrix    	#lager en vektor med lengde matrix med hvert element 1/matrix^1/2
            A = Matrix_Generator([4,11,4],matrix)
            egen, it = Power_Eig(A, x)
            powerfile.write(str(it)+"\t"+str(egen)+"\n")
        powerfile.close()
        return Power_Eig(matrix, x)
    else:
        return QR_Eig(matrix,len(matrix))




def Power_Eig(A, x):  # modified power method with normalization, p 399
    k_max = 100       # steps
    for k in range(k_max):
        y = np.dot(A, x)
        y_abs = nl.norm(y, np.inf)
        r = y[1] / x[1] # y |-> y[1] is arbitrary linear functional
        x = y / y_abs   # normalization
        print(k, x, r)
    print(0, x)


def QR_Eig(A,n):
		spectrum = [0]*n
    N = 0
    tol = 10**-14
    A = Hessenberg(A,n)
    for i in range(n,1):
        spectrum, A, t = QR_Shift(A,i,tol)
        N += t
    return spectrum,N


def Hessenberg(A,n): #fungere ikkje
    z = []
    e = [0]*n
    u = []
    e[0] = 1
    for k in range(1, n - 2):
        for i in range(k+1,n):
            z.append(A[i][k])
        for i in range(0,n):
            e.append((n-k))
        u = z + np.sign(z[0])*nl.norm(z)*e
        u /= nl.norm(u)
        for i in range(k + 1, n):
            for j in range(k,n):
                A[i][j] += -2*u*(np.transpose(u)*A[i][j])
        for i in range(1, n):
            for j in range(k+1,n):
                A[i][j] += -2*(A[i][j]*u)*np.transpose(u)
    return A


def QR_Shift(A,m,tol):
	la = A[m-1][m-1]
    t = 0
    e = 1 
    identitetsmatrise = np.identity(m)
    if m > 1:
        while e > tol:
            t += 1
            Q, R = sl.qr(A-la*identitetsmatrise)
            A = np.dot(R,Q) + la*identitetsmatrise
            la = A[m-1][m-1]
            e = A[m-1][m-2]
    return la,A,t


def Plot_Iterations(algorithm, matrices):  # takes in text document with x and y in colums. plot them
    x, y = [], []
    for line in open('power_1.txt', 'r'):
        values = [float(s) for s in line.split()]
        x.append(values[0])
        y.append(values[1])

    plt.plot(x, y)
    plt.show()
    return 0

def Plot_Iterations(algorithm, matrices):  # takes in text document with x and y in colums. plot them
    dim, r_A, r_B, r_C = [], [], [], []
    for line in open('power_1.txt', 'r'):
        values = [float(s) for s in line.split()]
        dim.append(values[0])
        r_A.append(values[1])
        r_B.append(values[2])
        r_C.append(values[3])

    fig, ax1 = plt.subplots()

    fig.suptitle(r'Largest eigenvalues of $A$, $B$, and $C$', fontsize=16)

    ax1.plot(dim, r_A, label=r'Numerical $r_{A}(n)$', lw=2)
    ax1.plot(dim, r_B, label=r'Numerical $r_{B}(n)$', lw=2)
    ax1.plot(dim, r_C, label=r'Numerical $r_{C}(n)$', lw=2)
    ax1.set_xlabel(r'$n $', fontsize=20)
    ax1.set_ylabel(r'$r(n) $', fontsize=20)
    ax1.legend(loc='best')
    ax1.grid()
    plt.show()
    return 0

#========================= Anders' forsøk på ting ========================= 

def QR_Eig(A, n):
    L = np.zeros(n)  # vector with zeroes
    N = 0            # number of steps in QR_shift
    tol = 10 ** -14  # tolerance
    A = Hessenberg(A, n)
    for i in range(n, 1):
        L[i], A, t = QR_Shift(np.resize(A, (i, i)), i, tol)
        N += t
    return L, N


def Hessenberg(A, n):  # Returns the Hessenberg form of a matrix
    z = []          #
    e = [0] * n     # first basis vector
    e[0] = 1
    for k in range(1, n - 2):
        for i in range(k + 1, n):
            z.append(A[i][k])
        np.resize(e, n-k)       # makes it same size as z
        u = z + (np.sign(z[0]) * nl.norm(z)) * e
        u = nl.norm(u)
        for i in range(k + 1, n):
            for j in range(k, n):
                A[i][j] += -2 * u * (np.transpose(u) * A[i][j])
        for i in range(1, n):
            for j in range(k + 1, n):
                A[i][j] += -2 * (A[i][j] * u) * np.transpose(u)
    return A


def QR_Shift(A, m, tol):
    la = A[m - 1][m - 1]  # lambda
    t = 0                 # number of iterations
    e = 1
    I = np.identity(m)    # mxm identity matrix
    if m > 1:
        while e > tol:
            t += 1
            Q, R, P = sl.qr(A - la*I, pivoting=False)
            A = np.dot(R, Q) + la * I
            la = A[m - 1][m - 1]
            e = A[m - 1][m - 2]
    return la, A, t
