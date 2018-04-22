from __future__ import division
################################################################################
import sys
import os
import csv
import numpy as np
import xml.etree.ElementTree as et
import matplotlib.pyplot as plt
from math import pi, sqrt
################################################################################

#heko
def Interpolation_Program(XMLFILE, PROGRAM, METHOD, GRID, n):
    #if PROGRAM == "Evaluation":
    f0, f1, I = XML_Extraction(XMLFILE)
    PREFIX = XMLFILE.split(".")[0]
    if(PROGRAM == "Evaluation"):
        X = Partition(GRID, I, n+1)
        Compute_Points(PREFIX, METHOD, GRID, f0, f1, I, X, n)
    else:
        DATA, T, PATH = Collect_Data(PREFIX,METHOD,GRID,4)
        if(PROGRAM == "Error"):
            Plot_Error(PATH, DATA, f0, T, I)
    return 0

def XML_Extraction(XMLFILE):
    file = et.parse(XMLFILE)
    root = file.getroot()
    f0 = lambda x : eval(root[0][0].text)
    f1 = lambda x : eval(root[0][1].text)
    I = [int(root[1][0].text), int(root[1][1].text)]
    return f0, f1, I

def Partition(GRID, I, n):
        partition = []
        if GRID == "Uniform":
            for i in range(0, n):
                partition.append(i/(n-1)*(I[1]-I[0]) + I[0]) #Alternative: 1/n (is without last point)
        else:
            for i in range(1, n+1):
                partition.append((1/2)*(I[0]+I[1]) + (1/2)*(I[1]-I[0]) * np.cos(((2*i-1) / (2*n)) * pi))
            '''for i in range(1, n):
                partition.append((1/2)*(I[0]+I[1]) + (1/2)*(I[1]-I[0]) * np.cos(((2*i+1) / (2*n+2)) * pi))'''
        return partition


def Compute_Points(PREFIX, METHOD, GRID, f0, f1, I, X, n):
    fname = "oppg3/"+PREFIX+"/"+METHOD+"_"+GRID+"/"+PREFIX+"_"+METHOD+"_"+GRID+"_"+str(n)+".txt"
    file = open(fname, "w")
    if METHOD == "Lagrange":
        a = Lagrange_Newton_Coefficients(f0, X, n)
        t = I[0]
        while t < I[1]:
            file.write(str(Lagrange_Newton_Evaluation(X, a, n, t)))
            file.write(";")
            t += 0.01
    else:
        a, Z = Hermite_Newton_Coefficients(f0, f1, X, n)
        t = I[0]
        notForste = False
        while t < I[1]:
            if(notForste):
                file.write(";")
            else:
                notForste = True
            file.write(str(Hermite_Newton_Evaluation(Z,a,n,t)))
            t += 0.01
    file.close()
    return

def Lagrange_Newton_Coefficients(f0, X, n):
    a = [0]*n
    for i in range(n):
        a[i] = f0(X[i])
    for j in range(1, n):
        for i in range(n - 1, j - 1, -1):
            print(X[i] - X[i-j])
            a[i] = (a[i] - a[i-1]) / (X[i] - X[i-j]);
    return a


def Lagrange_Newton_Evaluation(X, a, n, t):
    temp = a[n-1]
    for i in range(n - 1, -1, -1):
        temp = temp*(t - X[i]) + a[i]
    return temp


def Hermite_Newton_Coefficients(f0,f1,X,n):
    k = 2*n + 1
    Z = [0] * (k+1)
    Q = [[0]*(k+1) for i in range(k+1)]
    for i in range(n):
        Z[2*i] = X[i]
        Z[2*i + 1] = X[i]
        Q[2*i][0] = f0(X[i])
        Q[2*i + 1][0] = f0(X[i])
        Q[2*i][1] = f1(X[i])
        if i != 0:
            Q[2*i][1] = (Q[2*i][0] - Q[2*i - 1][0])/(Z[2*i] - Z[2*i-1])
    for i in range(2, k+1):
        for j in range(2, i+1):
            if(Z[i] == Z[i-j]):
                Q[i][j] = f1(Z[i])
            else:
                Q[i][j] = (Q[i][j-1] - Q[i-1][j-1]) / (Z[i] - Z[i-j])
    return Z, Q

    #https://docs.scipy.org/doc/numpy-1.14.0/reference/generated/numpy.polynomial.hermite.hermfit.html
    return a,Z


def Hermite_Newton_Evaluation(Q,a,n,t):
    k = 2 * n + 1
    s = Q[k][k] * (t - a[k-1])
    for i in range(2, k + 1):
        j = k - i + 1
        s = (s+Q[j][j]) * (t - a[j-1])
    s = s + Q[0][0]
    return s


def Collect_Data(PREFIX,METHOD,GRID,I):
    DATA = []
    PATH = "oppg3/"+PREFIX+"/"+METHOD+"_"+GRID+"/"
    for i in range(2, 21):
        fname = PATH + PREFIX+"_"+METHOD+"_"+GRID+"_"+str(i)+".txt"

        with open(fname) as file:
            points = file.readline().split(";")
            DATA.append(points)
    T = 0
    return DATA, T, PATH


def Plot_Error(PATH, DATA, f0, T, I):
    x = I[0]
    yPlot = []
    xPlot = []
    for j in range(len(DATA)):
        s = 0
        for i in range(len(DATA[j])):
            if(DATA[j][i] == ""):
                break
            x = i * 0.01 + I[0]
            s += (f0(x) - float(DATA[j][i]))**2
        s /= len(DATA[j])
        s = sqrt(s)
        yPlot.append(s)
        xPlot.append(j)
    plt.plot(xPlot, yPlot)
    print(yPlot)
    plt.show()
    return 0

def Plot_Polynomials(PATH,DATA,f0,X):
    return 0
for i in range(2, 21):
    Interpolation_Program("f1.xml", "Evaluation", "Lagrange", "Cheboshev", i)
Interpolation_Program("f1.xml", "Error", "Lagrange", "Cheboshev", 4)
#Interpolation_Program("f1.xml", "Evaluation", "Hermite", "Uniform", 4)