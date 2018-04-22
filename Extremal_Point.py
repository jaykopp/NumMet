from __future__ import division
################################################################################
import sys
import numpy as np
import numpy.linalg as nl
import scipy.linalg as sl
import xml.etree.ElementTree as et
################################################################################
from utils.XMLHelper import generateMatrixString, generateFunctionString, generateArrayString, generateAAndPrint

def Extremal_Points(XMLFILE,X0):
	n,f,Df,Hf = XML_Extraction(XMLFILE)
	#print(n, f, Df, Hf)
	#print(Hf([0.5, 0.5, 0.5]))
	X = Newton_Iteration(Df, Hf, X0)
	Classify_Point(f, Hf, X)

def vectorAbs(v):
	s = 0
	for el in v:
		s += el*el
	return s

def XML_Extraction(XMLFILE):
	tree = et.parse(XMLFILE)
	root = tree.getroot()
	dimension = int(root[0].text)
	diagonalElements = []
	upperTriangularElements = []
	gradient = []

	## Extract diagonal elements and gradient from XML root
	for i in range(dimension):
		diagonalElements.append(root[3][i].text)
		gradient.append(root[2][i].text)

	## Extract upper triangular elements
	for j in range( int((dimension*dimension - dimension)/2)):
		upperTriangularElements.append(root[3][i+j+1].text)

	## create lambda ready function string from xml data 
	functionLambdaString = generateFunctionString(root[1].text)

	## Create lambda ready array function string from gradient
	gradientLambdaString = generateArrayString(gradient)

	## Create lambda ready Hessian matrix string from 
	## diagonal and upper triangular elements
	hessianLambdaString = generateMatrixString(diagonalElements,upperTriangularElements)
	
	## Simply a print statement to test
	generateAAndPrint(upperTriangularElements, diagonalElements)
	
	## Create lambda expressions from strings
	hessianLambda = lambda X: np.array(eval(hessianLambdaString))
	derivedLambda = lambda X : np.array(eval(gradientLambdaString))
	functionLambda = lambda X: np.array(eval(functionLambdaString))
	#print(Hf([1,0]))
	return dimension,functionLambda,derivedLambda,hessianLambda 

#
# Calculates A^-1 * X using back-substituttion
#	A must be upper triangular
def backSubstitute(A, X):
	n = len(X)
	x = [0]*n
	n -= 1
	x[n] = X[n] / A[n][n]
	#x[n] = 5
	for i in range(n-1, -1, -1):
		s = X[i]
		for j in range(i+1, n):
			s -= A[i][j]*x[j]
		x[i] = s/A[i][i]
	return np.array(x)

def Newton_Iteration(Df,Hf,X0):
	fx = Df(X0)
	print(0, X0, fx)
	X = X0
	nmax = 500 # Set max iteration so it does not go on forever
	## Newton iteration where X = X - Hf(X)^(-1)*Df(X)
	for n in range(nmax):
		Dx = Df(X)
		Xn = nl.tensorsolve(Hf(X), Dx)
		X  = X - Xn
		print(n, X, Dx)
		if vectorAbs(Xn) < 10**(-14):
			print("Convergence")
			return X
	return X

#
# Classify input points based on eigenvalues
# If all eigenvalues are positive - X is minimum
# If all eigenvalues are negative - X is maximum
# if there are no zero eigenvalues - X is saddlepoint
# Otherwise, unable to classify
def Classify_Point(f,Hf,X):
	A = Hf(X)
	evalues = nl.eigvals(A)
	print(A)
	print(evalues)
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

#Extremal_Points("1c.xml", [-.5, -.5])
#Extremal_Points("1c.xml", [-1.0, -0.5])
#Extremal_Points("1c.xml", [-1.95, -0.05])
