#================================== OPPGAVE 2 ==================================
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
    A = []    
    if matrix == 1:
        A = [4,11,4]
    elif matrix == 2:
        A = [2,-7,20,-7,2]
    elif matrix == 3:
        A = [6,-3,-7,19,-7,-3,6]
    
    
    if(algorithm == "Power_Eig"):           # og en string algorithm
        powerfile = open("Power_"+str(matrix)+".txt","w") #åpner fil for å legge til data
        for i in range(10,201):
            x = [1/i**(1/2)]*i                          #lager en vektor med lengde matrix med hvert element 1/matrix^1/2
            B = Matrix_Generator(A,i)                   # denne er bare her for å teste
            egen, it = Power_Eig(B, x)
            powerfile.write(str(i)+"\t"+str(it)+"\n") # legger til i en fil
        powerfile.close()
        return 0
    elif(algorithm == "QR_Eig"):
        QRfile = open("QR_"+str(matrix)+".txt", "w")
        for i in range(10,201):
            B = Matrix_Generator(A,i) 
            eig, it = QR_Eig(B, i)
            QRfile.write(str(i)+"\t"+str(it)+"\n")
        QRfile.close()
        return 0
    else: return 0


def Power_Eig(A,x):
    r = 0
    it = 0
    err = 1
    while err > 10**-14:
        it += 1
        y = np.dot(A, x)
        r = y[0] / x[0]
        y = y/nl.norm(y, np.inf)
        err = nl.norm(x - y,np.inf)                 # y |-> y[1] is arbitrary linear functional
        x = y                                       
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
        z = A[k + 1:n+1, k]
        e = np.array([0] * (n-k-1))
        e[0] = 1
        u = z + (np.sign(z[0]) * np.sqrt(z.dot(z))) * e
        u = np.asarray(u / np.sqrt(u.dot(u)))
        A[k + 1:n+1, k:n+1] = A[k + 1:n+1, k:n+1] - 2 * np.outer(u, np.dot(np.transpose(u), A[k + 1:n+1, k:n+1]))
        A[1:n+1, k + 1:n+1] = A[1:n+1, k + 1:n+1] - 2 * np.outer(np.dot(A[1:n+1, k + 1:n+1], u), np.transpose(u))
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
        plt.savefig("Power_plot.png",transparent = True)	
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
        plt.savefig("QR_plot.png", transparent = True)
        plt.show()
	return 0
	
    else:
	return 0

def Plot_Iterations(algorithm, matricies): # matricies = [1], [2], [3], [1, 2], [1, 3], [2, 3], [1,2,3]
    data = [] # eller [[]*(len(matricies))] ?
    t = 0
    fig, ax1 = plt.subplots()
    if algorithm == "Power_Eig":
        fig.suptitle(r'Iterations over matrices of size $n$', fontsize=16)
        for i in matricies:
            size = []
            iter = []
            for line in open("Power_"+str(i)+".txt", "r"):
                d = [float(s) for s in line.split()]
                size.append(d[0])
                iter.append(d[1])
            data.append(size)
            data.append(iter)
            
            ax1.plot(data[t], data[t+1]) 
            t += 2
        """
        for i in len(matricies):
            if matricies[i] == 1:
                matricies[i] = "A"
            elif matricies[i] == 2:
                matricies[i] = "B"
            elif matricies[i] == 3:
                matricies[i] = "C"
        """
        ax1.set_xlabel(r'$size of matrix$', fontsize=20)
        ax1.set_ylabel(r'$iterations$', fontsize=20)
        ax1.legend(loc='best')
        ax1.grid()
        plt.savefig("Power_plot.png", transparent=True)
        plt.show()
        return 0
    
    elif algorithm == "QR_Eig":
        fig.suptitle(r'Iterations over matrices of size $n$', fontsize=16)
        for i in matricies:
            size = []
            iter = []
            for line in open("QR_" + str(i) + ".txt", "r"):
                d = [float(s) for s in line.split()]
                size.append(d[0])
                iter.append(d[1])
            data.append(size)
            data.append(iter)

            ax1.plot(data[t], data[t+1])
            t += 2

        ax1.set_xlabel(r'$size of matrix$', fontsize=20)
        ax1.set_ylabel(r'$iterations$', fontsize=20)
        ax1.legend(loc='best')
        ax1.grid()
        plt.savefig("QR_plot.png", transparent=True)
        plt.show()
