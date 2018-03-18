from __future__ import division
################################################################################
import sys
import numpy as np
import numpy.linalg as nl
import scipy.linalg as sl
import xml.etree.ElementTree as et
################################################################################


def Extremal_Points(XMLFILE,X0):
	n,f,Df,Hf = XML_Extraction(XMLFILE)
	#print(n, f, Df, Hf)
	#print(Hf([0.5, 0.5, 0.5]))
	X = Newton_Iteration(Df, Hf, [-2,-2,2])
	Classify_Point(f, Hf, X)
	#Classify_Point(f, Hf, [-2, 0])
	#Classify_Point(f, Hf, [0, 0])
def vectorAbs(v):
	s = 0
	for el in v:
		s += el*el
	return s

def generateMatrixString(d,u):
	s = "["
	n = len(d)
	j = 0
	for i in range(n):
		if i > 0:
			s += ","
		s += "["
		for k in range(0, i):
			s += "0,"
		s += d[i]
		for k in range(i, n -1):
			s += "," + u[j]
			j += 1
		s += "]"
	s += "]"
	return stringConvertVariablesToMultivariable(s)
def generateArrayString(gradient):
	s = "["
	for i in range(len(gradient)):
		if i != 0:
			s += ","
		s += gradient[i]
	s += "]"
	return stringConvertVariablesToMultivariable(s)

def stringConvertVariablesToMultivariable(s):
	s = s.replace("x", "X[0]")
	s = s.replace("y", "X[1]")
	s = s.replace("z", "X[2]")
	return s

def generateFunctionString(func):
	return stringConvertVariablesToMultivariable(func)
def generateAAndPrint(u, d):
	n = len(d)
	A = [[0]*n for i in range(n)]
	'''
	currU = len(u) - 1
	tempShift = 2
	for i in range(n-1, 0, -1):
		## Backwards
		t = u[currU]+"/"+d[i-1]
		for k in range(i, n):
			d[k] += " + "+t+"*"+u[currU + k - i]
		currU -= tempShift
		tempShift += 1'''
	## Adding diagonal
	k = 0
	## Total number: n + n(n+1)/2 -> O((n^2+n)/2)
	for i in range(n):
		A[i][i] = d[i]
		for j in range(i+1, n):
			A[i][j] = u[k]
			A[j][i] = u[k]
			k += 1
	print("Printing all the derivatives")
	for line in A:
		print(line)
	print("End printing all the derivatives")
def XML_Extraction(XMLFILE):
	tree = et.parse(XMLFILE)
	root = tree.getroot()
	n = int(root[0].text)
	d = []
	u = []
	grad = []
	for i in range(n):
		d.append(root[3][i].text)
		grad.append(root[2][i].text)
	for j in range( int((n*n - n)/2)):
		u.append(root[3][i+j+1].text)
	fs = generateFunctionString(root[1].text)
	Dfs = generateArrayString(grad)
	HS = generateMatrixString(d,u)
	Hf = lambda X: np.array(eval(HS))
	#Hf = np.array(HS)
	generateAAndPrint(u, d)
	print(HS)
	Df = lambda X : np.array(eval(Dfs))
	f = lambda X: np.array(eval(fs))
	#print(Hf([1,0]))
	return n,f,Df,Hf 

def backSubstitute(A, X):
	n = len(X)
	x = [0]*n
	n -= 1
	x[n] = X[n] / A[n][n]
	#x[n] = 5
	for i in range(n-1, -1, -1):
		s = X[i]
		for j in range(i+1, n+1):
			s -= A[i][j]*x[j]
		x[i] = s/A[i][i]
	return np.array(x)

def Newton_Iteration(Df,Hf,X0):
	fx = Df(X0)
	print(0, X0, fx)
	X = X0
	nmax = 500 # Set max iteration so it does not go on forever
	for n in range(10):
		Dx = Df(X)
		Xn = backSubstitute(Hf(X), Dx)
		X  = X - Xn
		print(n, X, Dx)
		if vectorAbs(Xn) < 10**(-14):
			print("Convergence")
			return X
	return X

def Classify_Point(f,Hf,X):
	##Firstly create eigenvalues of HF(X)
	A = Hf(X)
	evalues = nl.eigvals(A)
	#print("Evs",evalues, A)
	allPositives = True
	allNegatives = True
	noNulls = True
	for e in evalues:
		if e > 0:
			allNegatives = False
		elif e < 0:
			allPositives = False
		else:
			allPositives = False
			allNegatives = False
			noNulls = False
	if allPositives:
		print(X, "Is a minimum")
	elif allNegatives:
		print(X, "Is a maximum")
	elif noNulls:
		print(X, "is a Saddlepoint")
	else:
		print(X, "Is unclassified")


Extremal_Points("1f.xml", 0)
