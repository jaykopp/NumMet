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
        a = I[0]
        b = I[1]
        if GRID == "Uniform":
            for i in range(0, n):
                partition.append(i/(n-1)*(b-a) + a)
        else:
            for i in range(0, n):
                partition.append((1/2)*(a+b) + (1/2)*(b-a) * np.cos(((2*i+1) / (2*n)) * math.pi))
        return partition


def Compute_Points(PREFIX, METHOD, GRID, f0, f1, I, X, n):
    fname = "oppg3/"+PREFIX+"/"+METHOD+"_"+GRID+"/"+PREFIX+"_"+METHOD+"_"+GRID+"_"+str(n-1)+".txt"
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
    return temp


def Hermite_Newton_Coefficients(f0,f1,X,n):
    a = np.zeros(2*(n))
    Z = np.zeros(2*(n))
    for i in range(0,n):
        Z[2*i] = X[i]
        Z[2*i+1] = X[i] #X =  [3, 4, 5]
                        #Z = [3, 3, 4, 4, 5, 5]
    for i in range(2*n):
        a[i] = f0(Z[i]) #a = [f(3), f(3), f(4), f(4), f(5), f(5)]

    for j in range(1,2*(n)):
        for i in range(2*(n)-1, j-1, -1):
            if math.fabs(Z[i]-Z[i-j]) < 10**(-10): #Almost identical
                a[i] = (a[i]-a[i-1])/f1(Z[j-1])
            else:
                a[i] = (a[i]-a[i-1])/(Z[i]-Z[i-j])
    print(a, Z)
    return a,Z


def Hermite_Newton_Evaluation(Z,a,n,t):
    H = a[2*n-1]
    for i in range(2*n-1,-1,-1):
        H = H*(t-Z[i])+a[i]
    return H


def Collect_Data(PREFIX,METHOD,GRID,I):
    DATA = []
    PATH = "oppg3/"+PREFIX+"/"+METHOD+"_"+GRID+"/"
    for i in [2,4,6,8]:
        fname = PATH + PREFIX+"_"+METHOD+"_"+GRID+"_"+str(i)+".txt"

        with open(fname) as file:
            points = file.readline().split(";")
            DATA.append(points[:-1])# Remove final data point (it is empty)
    T = 8
    return DATA, T, PATH

def Plot_Error(PATH, DATA, f0, T, I):
    plt.figure()
    NAME = PATH.replace("oppg3/", "")
    NAME = NAME.replace("/", "_")
    NAME = NAME + "_error.png"
    plt.xlabel('Plynomial Degree')
    plt.ylabel('Error')

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
        xPlot.append(2*j + 2)  # Add degree of polynom (starts at 2 not 0, so add 2)
    plt.plot(xPlot, yPlot) 
    plt.savefig("img/"+NAME)
    return 0

def Plot_Polynomials(PATH,DATA,f0,I):
    plt.subplots()
    #plt.suptitle('Visualization for ' + PATH)
    NAME = PATH.replace("oppg3/", "")
    NAME = NAME.replace("/", "_")
    NAME = NAME + "_graph.png"
    xPlot = []
    yPlot = []
    zPlot = []

    for i in range(0, len(DATA)):
        # Loop through datasets (each different polynom)
        yPlot = []
        xPlot = []

        for j in range(0,len(DATA[i])):
            if i == 0:
                # Though first loop, add target value function and x values
                x = j * 0.01 + I[0]
                zPlot.append(f0(x))

            xPlot.append(j * 0.01 + I[0])
            ## Append value at point
            yPlot.append(float(DATA[i][j]))
        if i == 0:
            # Plot the target function
            plt.plot(xPlot, zPlot, '--', label='Target function', lw=2, color="black")
        # Plot the dataset
        plt.plot(xPlot, yPlot, label='Degree '+str(2*i + 2), lw=2)
    
    plt.legend()
    plt.savefig("img/"+NAME)
    return 0

for file in ["f1.xml", "f2.xml", "f3.xml", "f4.xml"]:
<<<<<<< HEAD
    for method in ["Hermite", "Lagrange"]:
=======
    for method in ["Lagrange", "Hermite"]:
>>>>>>> 363be2c569d94b26f1026d5dfeda239601d9be6c
        for partition in ["Chebyshev", "Uniform"]:

            for i in [2,4,6,8]:
                Interpolation_Program(file, "Evaluation", method, partition, i)
            #Interpolation_Program(file, "Error", method, partition, 4)
            Interpolation_Program(file, "Visualization", method, partition, 4)

#Interpolation_Program("f1.xml", "Evaluation", "Hermite", "Chebyshev", 4)