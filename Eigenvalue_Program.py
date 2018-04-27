# ================================== OPPGAVE 2 ==================================
################################################################################
import sys
import os
import numpy as np
import numpy.linalg as nl
import scipy.linalg as sl
import scipy.sparse as sp
import matplotlib.pyplot as plt


################################################################################


def Eigenvalue_Program(program, algorithm,
                       matrix):  # tar inn hvilket program og hvilken algoritme, samt en 1 til 3 lang array med tall fra 1 til tre, som har en tilhorende matrise !, B eller C
    if program == "Run_Simulation":
        for i in matrix:  # itererer gjennom array-en og kjorer simulering pa alle tall med i den
            Run_Simulation(i, algorithm)
        return 0
    elif program == "Plot_Iterations":  # plotter matrisene tilsvarende tallene
        Plot_Iterations(algorithm, matrix)
        return 0
    else:
        return 0


def Matrix_Generator(matrix, n):  # tar inn en liste matrix og et tall n
    A = np.zeros((n, n))  # fyller A med 0
    forskjell = -int(len(matrix) / 2)
    for i in matrix:
        A += i * np.eye(n, n, forskjell)  # legger til i ganger forskjell pa diagonalen
        forskjell += 1
    return A


def Run_Simulation(matrix, algorithm):  # tar inn et tall matrix tilsvarende en matrise og en string algorithm
    A = []
    if matrix == 1:  # legger matrise til riktig tall
        A = [4, 11, 4]
    elif matrix == 2:
        A = [2, -7, 20, -7, 2]
    elif matrix == 3:
        A = [6, -3, -7, 19, -7, -3, 6]

    if (algorithm == "Power_Eig"):  # sjekker hvilken algoritme man skal bruke
        powerfile = open("Power_" + str(matrix) + ".txt", "w")  # apner fil for a legge til data
        for i in range(10, 201):  # itererer fra 10 til 2000
            x = [1 / i ** (1 / 2)] * i  # lager en vektor med lengde i med hvert element 1/matrix^1/2
            B = Matrix_Generator(A, i)  # Genererer en matrise av type "matrix" med dimensjon ixi
            egen, it = Power_Eig(B, x)  # kjorer algoritmen pa B
            powerfile.write(str(i) + "\t" + str(it) + "\n")  # legger til i en fil
        powerfile.close()  # lukker filen
        return 0
    elif (algorithm == "QR_Eig"):  # tilsvarende her for QR-filen
        QRfile = open("QR_" + str(matrix) + ".txt", "w")
        for i in range(10, 201):
            B = Matrix_Generator(A, i)
            eig, it = QR_Eig(B, i)
            QRfile.write(str(i) + "\t" + str(it) + "\n")
        QRfile.close()
        return 0
    else:
        return 0


def Power_Eig(A, x):
    r = 0  # egenverdi, ikke bestemt enna
    it = 0  # iterasjoner
    err = 1  # error
    while err > 10 ** -14:
        it += 1
        y = np.dot(A, x)  # Prikkprodukt mellom matrise A og vektor x
        r = y[0] / x[0]  # setter egenverdi
        y = y / nl.norm(y, np.inf)  # normaliserer
        err = nl.norm(x - y, np.inf)  # y |-> y[1] is arbitrary linear functional
        x = y
    return r, it


def QR_Eig(A, n):
    L = np.zeros(n)  # vector with zeroes
    N = 0  # number of steps in QR_shift
    tol = 10 ** -14  # tolerance
    A = Hessenberg(np.array(A), n)
    for i in range(n - 1, 0, -1):
        L[i - 1], A, t = QR_Shift(np.resize(A, (i, i)), i, tol)
        N += t
    return L, N


def Hessenberg(A, n):  # returnerer hessenbergformen til en matrisen A
    for k in range(1, n - 2):
        z = A[k + 1:n + 1, k]
        e = np.array([0] * (n - k - 1))
        e[0] = 1
        u = z + (np.sign(z[0]) * np.sqrt(z.dot(z))) * e
        u = np.asarray(u / np.sqrt(u.dot(u)))
        A[k + 1:n + 1, k:n + 1] = A[k + 1:n + 1, k:n + 1] - 2 * np.outer(u, np.dot(np.transpose(u),
                                                                                   A[k + 1:n + 1, k:n + 1]))
        A[1:n + 1, k + 1:n + 1] = A[1:n + 1, k + 1:n + 1] - 2 * np.outer(np.dot(A[1:n + 1, k + 1:n + 1], u),
                                                                         np.transpose(u))
    return A


def QR_Shift(A, m, tol):  # QR-shifter matrise A
    la = A[m - 1][m - 1]
    t = 0
    e = 1
    I = np.identity(m)
    if m > 1:
        while e > tol:
            t += 1
            Q, R = sl.qr(A - (la) * I)
            A = np.dot(R, Q) + la * I
            la = A[m - 1][m - 1]
            e = A[m - 1][m - 2]
    return la, A, t


def Plot_Iterations(algorithm, matricies):  # matricies = [1], [2], [3], [1, 2], [1, 3], [2, 3], [1,2,3]
    data = []  # data far vekslende kolonner med input og output av enten Power eller QR
    t = 0
    fig, ax1 = plt.subplots()
    if algorithm == "Power_Eig":
        for i in matricies:
            size = []
            iter = []
            for line in open("Power_" + str(i) + ".txt", "r"):
                d = [float(s) for s in line.split()]
                size.append(d[0])
                iter.append(d[1])
            data.append(size)
            data.append(iter)

            ax1.plot(data[t], data[t + 1])
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
        ax1.set_xlabel(r'size of matrix', fontsize=10)
        ax1.set_ylabel(r'iterations', fontsize=10)
        ax1.legend(loc='best')
        plt.savefig("Power_plot.png", transparent=True)
        plt.show()
        return 0

    elif algorithm == "QR_Eig":
        for i in matricies:
            size = []
            iter = []
            for line in open("QR_" + str(i) + ".txt", "r"):
                d = [float(s) for s in line.split()]
                size.append(d[0])
                iter.append(d[1])
            data.append(size)
            data.append(iter)

            ax1.plot(data[t], data[t + 1])
            t += 2

        ax1.set_xlabel(r'size of matrix', fontsize=10)
        ax1.set_ylabel(r'iterations$', fontsize=10)
        ax1.legend(loc='best')
        plt.savefig("QR_plot.png", transparent=True)
        plt.show()
        return 0


#Eigenvalue_Program("Run_Simulation", "Power_Eig", [1, 2, 3])
#Eigenvalue_Program("Run_Simulation", "QR_Eig", [2, 3])
Eigenvalue_Program("Plot_Iterations", "Power_Eig", [1, 2, 3])
Eigenvalue_Program("Plot_Iterations", "QR_Eig", [2, 3])
