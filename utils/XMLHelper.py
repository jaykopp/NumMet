
def generateMatrixString(d,u):
	s = "["
	n = len(d)
	j = 0
	l = 0
	for i in range(n):
		if i > 0:
			s += ","
		s += "["
		for k in range(0, i):
			s += u[l] + ","
			l += 1
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
	k = 0
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
