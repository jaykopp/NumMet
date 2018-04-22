from __future__ import division
################################################################################
import sys
import os
import numpy as np
import xml.etree.ElementTree as et
import matplotlib.pyplot as plt
from math import pi
################################################################################


def Interpolation_Program(XMLFILE, PROGRAM, METHOD, GRID, n):
    #if PROGRAM == "Evaluation":
    f0, f1, I = XML_Extraction(XMLFILE)
    if(PROGRAM == "Evaluation"):
        PREFIX = "f1"
        Compute_Points(PREFIX, METHOD, GRID, f0, f1, I, 0, n)
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
    #for i in range()
    #
    grid = Partition(GRID, I, 5)
    fname = "oppg3/"+PREFIX+"_"+METHOD+"_"+GRID+"_4.txt"
    file = open(fname, "w")
    x = I[0]
    while x < I[1]:
        x += 0.01
        file.write(str(f0(x)))
        file.write(";")
    file.close()
    return "H"


def Lagrange_Newton_Coefficients(f0, X, n):
    X.astype(float)
    f0.astype(float)
    a = []
    for i in range(n):
        a.append(f0[i])
    for j in range(1, n):
        for k in range(n - 1, j - 1, -1):
            a[k] = (a[k] - a[k - 1]) / (X[k] - X[k - j])
    return a


def Lagrange_Newton_Evaluation(X, a, n, t):
    ## TODO: Fix this
    X.astype(float)
    temp = a[n]
    for i in range(n - 1, 0, -1):
        temp = temp*(t - X[i]) + a[i]
    P = temp
    return P


def Hermite_Newton_Coefficients(f0,f1,X,n):
    #https://docs.scipy.org/doc/numpy-1.14.0/reference/generated/numpy.polynomial.hermite.hermfit.html
    return a,Z


def Hermite_Newton_Evaluation(Z,a,n,t):
    return H


def Collect_Data(PREFIX,METHOD,GRID,I):
    return DATA,T,PATH


def Plot_Error(PATH, DATA, f0, T):
    return 0

def Plot_Polynomials(PATH,DATA,f0,X):
    return 0

Interpolation_Program("f1.xml", "Evaluation", "Lagrange", "Uniform", 10)