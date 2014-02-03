import sys
import math
import os

def coeff(st, c, k1, km1, k2, e0, s0, Km):
	ans = pow(st+1, c)*pow(e0+1, c)*1.0/math.factorial(c)
	ans = ans - pow(st+e0+2, c) + math.factorial(c)
	ans = ans / pow(Km, c)
	if (ans == 0):
		return 0.000000000001
	return ans

def coeff2 (s, c, k1, km1, k2, e0, s0, Km):
	p = 1
	if (c == 0):
		return 1
	for i in range (1, c+1):
		p *= (s-i+1)*(e0-i+1)*1.0/(i*(km1)/k1)
	if (p==0):
		return 0.000000000001
	return p	

def prop(s, k1, km1, k2, e0, s0, Km):
	if (s == s0+1):
		return 0
	r = (e0+Km+s)
	return 2*k2*e0*s/(r + math.sqrt(r*r-4*e0*s))
	
def main():

	fin = open (os.path.dirname(__file__) + "/parameters", 'r')

	k1 = float(fin.readline()) * 1.0
	km1 = float(fin.readline()) * 1.0
	k2 = float(fin.readline()) * 1.0

	e0 = int(fin.readline())
	s0 = int(fin.readline())

	Km = (km1+k2)/k1

	cwd = os.path.dirname( os.path.dirname(__file__) )
	f1_name = cwd + "/data/Cond_Exp_TQSSA"
	f2_name = cwd + "/data/Cond_Exp_Barik"
	f3_name = cwd + "/data/Cond_Exp_Substrate"

	f1 = open (f1_name, 'w')
	f2 = open (f2_name, 'w')
	f3 = open (f3_name, 'w')

	p = [0]*(s0+2)
	p[s0] = 1
	
	q = [0]*(s0+2)
	q[s0] = 1

	coeff_matrix = [[0 for _ in range(s0+1)] for _ in range(e0+1)]
	sum_coeff = [0 for _ in range(s0+1)]
	mean_coeff_matrix = [[0 for _ in range(s0+1)] for _ in range(e0+1)]
	expectation = [0 for _ in range(s0+1)]


	for i in range (0, s0+1):
		for j in range (0, e0+1):
			coeff_matrix[j][i] = coeff2(i,j, k1, km1, k2, e0, s0, Km)

	for i in range (0, s0+1):
		for j in range (0, min(i, e0)+1):
			sum_coeff[i] += coeff2(i, j, k1, km1, k2, e0, s0, Km)

	for i in range (0, s0+1):
		for j in range (0, e0+1):
			mean_coeff_matrix[j][i] = j*coeff2(i, j, k1, km1, k2, e0, s0, Km)

	for i in range (0, s0+1):
		for j in range (0, min(i, e0)+1):
			expectation[i] += mean_coeff_matrix[j][i]/sum_coeff[i]


	for i in range (0, s0+1):
		f1.write ( str(prop(i, k1, km1, k2, e0, s0, Km)/k2) + ' ')

	for i in range (0, s0+1):
		f2.write ( str(expectation[i]) + ' ')

	for i in range (0, s0+1):
		f3.write ( str(i) + ' ')

	f1.close()
	f2.close()
	f3.close()

if __name__ == '__main__':
	main()
