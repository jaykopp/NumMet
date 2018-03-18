from __future__ import division
################################################################################
import sys
import os
import numpy as np
import xml.etree.ElementTree as et
import matplotlib.pyplot as plt
################################################################################


def Interpolation_Program(XMLFILE, PROGRAM, METHOD, GRID, n):
	#if PROGRAM == "Evaluation":
	f, I = XML_Extraction(XMLFILE)
	Partition("Uniform", I, 5)
	return 0

def XML_Extraction(XMLFILE):
    file = et.parse(XMLFILE)
    root = file.getroot()
    f = lambda x : eval(root[0][0].text)
    I = [int(root[1][0].text), int(root[1][1].text)]
    return f, I

def Partition(GRID, I, n):
		partition = []
		if GRID == "Uniform":
			print(range(I[0], I[1], 1))
			#for i in range(I[0], I[1], 0.01):
			#	partition.append
		#for i in range(0, n):
		#    GRID.append((1/2)*(I[0]+I[1]) + (1/2)*(I[1]-I[0]) * np.cos(((2*i+1) / (2*n+2)) * pi))
		return 0


def Compute_Points(PREFIX, METHOD, GRID, f0, f1, I, X, n):
	#for i in range()
	return Q


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

Interpolation_Program("f1.xml", "Evaluation", "", "", 10)