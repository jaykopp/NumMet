

def subarray(A, k, n):
	x = []
	for i in range(k, n):
		x.append(A[i][k])
	return x
def Hessenberg(A, n):
	for k in range(0, n-2):
		z = subarray(A, k+1, n)
		e = [0]*(n-k)
		e[0] = 1
		u = z + sgn(z[])
def QR_Eigenvalue(A):
	n = len(A)
	L = [0]*n
	N = 0
	tol = tolerance...
	A = Hessenberg(A,n)