from __future__ import division
################################################################################
import sys
import os
import csv
import numpy as np
import xml.etree.ElementTree as et
import matplotlib.pyplot as plt
import math
################################################################################

##
## Main program to perform interpolation
## n is the degree of the polynom
##
def Interpolation_Program(XMLFILE, PROGRAM, METHOD, GRID, n):
    
    f0, f1, I = XML_Extraction(XMLFILE)
    
    n += 1 ## Add 1 to degree to include extra point in partition

    PREFIX = XMLFILE.split(".")[0]

    if(PROGRAM == "Evaluation"):
        X = Partition(GRID, I, n)
        Compute_Points(PREFIX, METHOD, GRID, f0, f1, I, X, n)
    else:
        DATA, T, PATH = Collect_Data(PREFIX,METHOD,GRID,4)
        if(PROGRAM == "Error"):
            Plot_Error(PATH, DATA, f0, T, I)
        else:
            Plot_Polynomials(PATH,DATA,f0,I)
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
        e = I[1]
        s = I[0]
        if GRID == "Uniform":
            for i in range(0, n):
                partition.append(i/(n-1)*(e-s) + s)
        else:
            for i in range(1, n):
                partition.append((1/2)*(s+e) + (1/2)*(e-s) * np.cos(((2*i-1) / (2*n)) * math.pi))
        print(partition)
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
    n = len(X)
    a = [0]*n
    for i in range(n):
        a[i] = f0(X[i])
    for j in range(1, n):
        for i in range(n-1, j - 1, -1):
            a[i] = (a[i] - a[i-1]) / (X[i] - X[i-j]);
    return a


def Lagrange_Newton_Evaluation(X, a, n, t):
    n = len(X)
    temp = a[n-1]
    for i in range(n - 2, -1, -1):
        temp = temp*(t - X[i]) + a[i]
    #print(temp)
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
    for i in range(2, 7):
        fname = PATH + PREFIX+"_"+METHOD+"_"+GRID+"_"+str(i)+".txt"

        with open(fname) as file:
            points = file.readline().split(";")
            DATA.append(points[:-1])# Remove final data point (it is empty)
    T = 8
    return DATA, T, PATH


def Plot_Error(PATH, DATA, f0, T, I):
    plt.figure()
    plt.xlabel('Plynomial Degree')
    plt.ylabel('Error')
    plt.title('Error for polynomials in ' + PATH)

    yPlot = []
    xPlot = []
    x = I[0] # start position

    for j in range(len(DATA)):
        s = 0 # sumation
        for i in range(len(DATA[j])):
            x = i * 0.01 + I[0]
            s += (f0(x) - float(DATA[j][i]))**2
        s /= len(DATA[j])    # Divide sum by lengt to get average
        error = math.sqrt(s) # Square s to get square error
        yPlot.append(error)  # Add error to list for graph
        xPlot.append(j + 2)  # Add degree of polynom (starts at 2 not 0, so add 2)
    plt.plot(xPlot, yPlot) 
    plt.show()
    return 0

def Plot_Polynomials(PATH,DATA,f0,I):
    plt.subplots()
    plt.xlabel('X values')
    plt.ylabel('Y Values')
    plt.title('Visualization for ' + PATH)

    xPlot = []
    yPlot = []
    zPlot = []

    for i in range(len(DATA)):
        # Loop through datasets (each different polynom)
        yPlot = [] 
        for j in range(len(DATA[i])):

            if i == 0:
                # Though first loop, add target value function and x values
                x = j * 0.01 + I[0]
                xPlot.append(x)
                zPlot.append(f0(x))

            ## Append value at point
            yPlot.append(float(DATA[i][j]))
        # Plot the dataset
        plt.plot(xPlot, yPlot)
    # Plot the target function
    plt.plot(xPlot, zPlot)
    plt.show()
    return 0
for i in range(1, 19):
    Interpolation_Program("f1.xml", "Evaluation", "Hermite", "Uniform", i)
Interpolation_Program("f1.xml", "Error", "Hermite", "Uniform", 4)
Interpolation_Program("f1.xml", "Visualization", "Hermite", "Uniform", 4)

#Interpolation_Program("f1.xml", "Evaluation", "Hermite", "Uniform", 4)