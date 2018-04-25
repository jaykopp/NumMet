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
    if(algorithm == "Power_Eig"):       # og en string algorithm
        powerfile = open("Power_1.txt","w") #åpner fil for å legge til data
        for i in range(10,matrix):
            x = [1/i**(1/2)]*i   #lager en vektor med lengde matrix med hvert element 1/matrix^1/2
            A = Matrix_Generator([4,11,4],i) # denne er bare her for å teste
            B = Matrix_Generator([2,-7,20,-7,2],i) # vet ikke om jeg skal beholde dem her eller gjøre det på en bedre måte
            C = Matrix_Generator([6,-3,-7,19,-7,-3,6],i) # todo er å gjøre denne bedre
            egen1, it1 = Power_Eig(A, x)
            egen2, it2 = Power_Eig(B, x)
            egen3, it3 = Power_Eig(C, x)
            powerfile.write(str(i)+"\t"+str(it1)+"\t"+str(it2)+"\t"+str(it3)+"\n") # legger til i en fil
        powerfile.close()
        return 0
    elif(algorithm == "QR_Eig"):
        QRfile = open("QR_1.txt", "w")
        for i in range(10,matrix):
            B = Matrix_Generator([2,-7,20,-7,2],i) # vet ikke om jeg skal beholde dem her eller gjøre det på en bedre måte
            C = Matrix_Generator([6,-3,-7,19,-7,-3,6],i) # todo er å gjøre denne bedre
            eig1, it1 = QR_Eig(B, i)
            eig2, it2 = QR_Eig(C, i)
            QRfile.write(i)+"\t"+str(it1)+"\t"+str(i)+"\t"+str(it2)+"\n")
        QRfile.close()
        return 0
    else: return 0


def Power_Eig(A,x): 		#Jacobs Power_Eig
    r=0
    for it in range(10000):
        x_ny = np.dot(A, x)
        x_ny_abs = nl.norm(x_ny, np.inf)
        r = x_ny[1]/x[1]
        x = x_ny / x_ny_abs
    return r, it


def Power_Eig(A,x):
    r = 0
    it = 0
    
    while nl.norm(np.dot(A,x) - r*x)- > 10**-14:
        it += 1
        y = np.dot(A, x)
        y_abs = nl.norm(y, np.inf)
        r = y[1] / x[1] 	# y |-> y[1] is arbitrary linear functional
        x = y / y_abs   	# normalization
    return r, it



def QR_Eig(A, n):
    L = np.zeros(n)         # vector with zeroes
    N = 0                   # number of steps in QR_shift
    tol = 10**-14           # tolerance
    A = Hessenberg(np.array(A), n)
    for i in range(n-1, 0, -1):
        L[i-1], A, t = QR_Shift(np.resize(A, (i, i)), i, tol)
        N += t
    return L, N


def Hessenberg(A, n):  # Returns the Hessenberg form of a matrix
    for k in range(1, n - 2):
        z = A[k + 1:n, k]
        e = np.array([0] * (n-k))
        e[0] = 1
        u = z + (np.sign(z[0]) * z.dot(z)) * e
        u = np.asarray(u / u.dot(u))
        A[k + 1:n, k:n] = A[k + 1:n, k:n] - 2 * np.outer(u, np.dot(np.transpose(u), A[k + 1:n, k:n]))
        A[1:n, k + 1:n] = A[1:n, k + 1:n] - 2 * np.outer(np.dot(A[1:n, k + 1:n], u), np.transpose(u))
    return A


def QR_Shift(A, m, tol):
    la = A[m - 1][m - 1]  # lambda
    t = 0                 # number of iterations
    e = 1
    I = np.identity(m)    # mxm identity matrix
    if m > 1:
        while e > tol:
            t += 1
            Q, R = sl.qr(A - (la)*I)
            A = np.dot(R, Q) + la * I
            la = A[m - 1][m - 1]
            e = A[m - 1][m - 2]
    return la, A, t


def Plot_Iterations(algorithm, matrices):  # takes in text document with x and y in colums. plot them
    if algorithm == "Power_Eig":
        dim, r_A, r_B, r_C = [], [], [], []
        for line in open('power_1.txt', 'r'):
            values = [float(s) for s in line.split()]
            dim.append(values[0])
            r_A.append(values[1])
            r_B.append(values[2])
            r_C.append(values[3])

        fig, ax1 = plt.subplots()

        fig.suptitle(r'Eigenvalues of $A$, $B$, and $C$', fontsize=16)

        ax1.plot(dim, r_A, label=r'Numerical $r_{A}(n)$', lw=2)
        ax1.plot(dim, r_B, label=r'Numerical $r_{B}(n)$', lw=2)
        ax1.plot(dim, r_C, label=r'Numerical $r_{C}(n)$', lw=2)
        ax1.set_xlabel(r'$n$', fontsize=20)
        ax1.set_ylabel(r'$r(n)$', fontsize=20)
        ax1.legend(loc='best')
        ax1.grid()
        plt.show()
	return 0
    elif algorithm == "QR_Eig":
        N_B, L0_B, N_C, L0_C = [], [], [], []
        for line in open('QR_1.txt', 'r'):
            values = [float(s) for s in line.split()]
            N_B.append(values[0])
            L0_B.append(values[1])
            N_C.append(values[2])
            L0_C.append(values[3])

        fig, ax1 = plt.subplots()

        fig.suptitle(r'$B$, and $C$', fontsize=16)

        ax1.plot(N_B, L0_B, label=r'Numerical $L_{B}(N)$', lw=2)
        ax1.plot(N_C, L0_C, label=r'Numerical $L_{C}(N)$', lw=2)
        ax1.set_xlabel(r'$ N $', fontsize=20)
        ax1.set_ylabel(r'$ lambdas $', fontsize=20)
        ax1.legend(loc='best')
        ax1.grid()
        plt.show()
	return 0
	
    return 0
