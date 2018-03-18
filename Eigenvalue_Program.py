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


def Run_Simulation(matrix,algorithm):
    if(algorithm == "Power_Eig"):
        x = []
        for i in range(len(matrix)):
            x.append(1 / (len(matrix)) ** (1 / 2))
        return Power_Eig(matrix, x)
    else:
        return QR_Eig(matrix,len(matrix))




def Power_Eig(A,x):
    def Power_Eig(A,x):
    r=0
    
    for it in range(100):
        x_ny = np.dot(A, x)
        x_ny_abs = nl.norm(x_ny, np.inf)
        r = x_ny[1]/x[1]
        x = x_ny / x_ny_abs
    return r, it


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
		la = A[m][m]
    t = 0
    e = 1 
    identitetsmatrise = np.identity(m)
    if m > 1:
        while e > tol:
            t += 1
            Q, R = sl.qr(A-la*identitetsmatrise)
            A = np.dot(R,Q) + la*identitetsmatrise
            la = A[m][m]
            e = A[m][m-1]
    return la,A,t


def Plot_Iterations(algorithm,matrices):
    return 0
