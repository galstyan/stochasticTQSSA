import sys
import math
import os
from copy import deepcopy

T_out = 25 # number of time points in output
N = 5000*(T_out-1) # number of steps in Euler method

def prop(s, k1, km1, k2, e0, s0):
	if (s == s0+1):
		return 0
	return k2*e0*s/(Km+s)
	
# mean number of molecules at time t over N sample path realizations
def mean (data, k1, km1, k2, e0, s0):
	mean = 0.0;
	for c in range (0, e0+1):
		for s in range (0,s0+1):
			mean = mean + s*data[s][c]
	return mean

# variance of molecules at time t over N sample path realizations	
def var (data, k1, km1, k2, e0, s0):
	var = 0.0
	m = mean(data, k1, km1, k2, e0, s0)
	for c in range (0, e0+1):
		for s in range (0, s0+1):
			var = var + (s-m)*(s-m)*data[s][c]
	
	return math.sqrt(var)

def expCgivenS (data, k1, km1, k2, e0, s0, s):
	sumS = 0.0
	for c in range (0, e0+1):
		sumS = sumS + data[s][c]
	cond = 0.0
	for c in range (0, e0+1):
		cond = cond + c * data[s][c]
	return cond/sumS

def ratioTest (data, k1, km1, k2, e0, s0, s):
	cond = expCgivenS (data, k1, km1, k2, e0, s0, s)
	sumS = 0.0
	for c in range (0, e0+1):
		sumS = sumS + data[s][c]
	sd = 0.0
	for c  in range (0, e0+1):
		sd = sd + (c-cond)*(c-cond)*data[s][c]
	sd = sd / sumS
	return sd/(e0*s)


def abs (s):
	if s < 0:
		return -s
	return s


def main():

	fin = open (os.path.dirname(__file__) + "/parameters", 'r')

	k1 = float(fin.readline())
	km1 = float(fin.readline())
	k2 = float(fin.readline())

	e0 = int(fin.readline())
	s0 = int(fin.readline())

	T = float(fin.readline())

	dt = T/N

	cwd = os.path.dirname( os.path.dirname(__file__) )
	f1_name = cwd + "/data/Exact_ODE_Substrate_Mean"
	f2_name = cwd + "/data/Exact_ODE_Substrate_SD"

	f3_1_name = cwd + "/data/Exact_ODE_Cond_1"
	f3_2_name = cwd + "/data/Exact_ODE_Cond_2"
	f3_3_name = cwd + "/data/Exact_ODE_Cond_3"

	f4_1_name = cwd + "/data/Exact_ODE_ratio_1"
	f4_2_name = cwd + "/data/Exact_ODE_ratio_2"
	f4_3_name = cwd + "/data/Exact_ODE_ratio_3"

	f5_name = cwd + "/data/Exact_ODE_ratio_tSubstrate"


	f1 = open (f1_name, 'w')
	f2 = open (f2_name, 'w')

	f3_1 = open (f3_1_name, 'w')
	f3_2 = open (f3_2_name, 'w')
	f3_3 = open (f3_3_name, 'w')

	f4_1 = open (f4_1_name, 'w')
	f4_2 = open (f4_2_name, 'w')
	f4_3 = open (f4_3_name, 'w')

	f5 = open (f5_name, 'w')


	# initializing the probability distribution
	p = [[0]*(e0+1) for i in range (s0+1)]
	p[s0][0] = 1
	
	q = [[0]*(e0+1) for i in range (s0+1)]
	q[s0][0] = 1
	
	s_mean_dynamics = []
	s_mean_dynamics.append(s0)
	
	s_var_dynamics = []	
	s_var_dynamics.append(0)
	
	for i in range (1,N+1):

		for c in range (1, e0):
			for s in range (c, s0):
				p[s][c] = q[s][c] + dt * ( 
					- (k1*(s-c)*(e0-c) + (km1+k2)*c)*q[s][c]
					+ k1 * (s-c+1)*(e0-c+1)*q[s][c-1]
					+ km1 * (c+1)*q[s][c+1]
					+ k2 * (c+1)*q[s+1][c+1]
					)

		c = 0
		s = 0
		p[s][c] = q[s][c] + dt * ( 
					- (k1*(s-c)*(e0-c) + (km1)*c)*q[s][c]
					+ km1 * (c+1)*q[s][c+1]
					+ k2 * (c+1)*q[s+1][c+1]
					)

		c = e0
		s = 0
		p[s][c] = q[s][c] + dt * ( 
					- ( (km1)*c)*q[s][c]
					)

		c = 0
		s = s0
		p[s][c] = q[s][c] + dt * ( 
					- (k1*(s-c)*(e0-c) + (km1+k2)*c)*q[s][c]
					+ km1 * (c+1)*q[s][c+1]
					)

		c = e0
		s = s0
		p[s][c] = q[s][c] + dt * ( 
					- ( (km1+k2)*c)*q[s][c]
					+ k1 * (s-c+1)*(e0-c+1)*q[s][c-1]
					)

		c = 0
		for s in range (1, s0):
			p[s][c] = q[s][c] + dt * ( 
					- (k1*(s-c)*(e0-c) + (km1+k2)*c)*q[s][c]
					+ km1 * (c+1)*q[s][c+1]
					+ k2 * (c+1)*q[s+1][c+1]
					)

		c = e0
		for s in range (e0, s0):
			p[s][c] = q[s][c] + dt * ( 
					- ( (km1+k2)*c)*q[s][c]
					+ k1 * (s-c+1)*(e0-c+1)*q[s][c-1]
					)

		s = s0
		for c in range (1, e0):
			p[s][c] = q[s][c] + dt * ( 
					- (k1*(s-c)*(e0-c) + (km1+k2)*c)*q[s][c]
					+ k1 * (s-c+1)*(e0-c+1)*q[s][c-1]
					+ km1 * (c+1)*q[s][c+1]
					)
		s = 0
		for c in range (1, e0):
			p[s][c] = q[s][c] + dt * ( 
					+ km1 * (c+1)*q[s][c+1]
					+ k2 * (c+1)*q[s+1][c+1]
					)

		s_mean_dynamics.append(mean(p, k1, km1, k2, e0, s0))
		s_var_dynamics.append(var(p, k1, km1, k2, e0, s0))

		q = deepcopy(p)
		
		sum = 0
		sum1 = 0
		for c in range (0, e0+1):
			for s in range (0, s0+1):

				sum = sum + p[s][c]
				sum1 = sum1 + abs(p[s][c])

	#	if (i % (N/10) == 0):
	#		print str(sum*0.5+sum1*0.5) + ":   " + str(i*100/N) + " percent completed"

		if (i % (N/10) == 0):
			print str(i*100/N) + " percent completed"


		if (i == 1*N/10):
			for l in range (0, s0+1):
				f3_1.write (str(expCgivenS (p, k1, km1, k2, e0, s0, l)) + ' ')

		if (i == 3*N/10):
			for l in range (0, s0+1):
				f3_2.write (str(expCgivenS (p, k1, km1, k2, e0, s0, l)) + ' ')

		if (i == 7*N/10):
			for l in range (0, s0+1):
				f3_3.write (str(expCgivenS (p, k1, km1, k2, e0, s0, l)) + ' ')



		if (i == 1*N/10):
			for l in range (1, s0+1):
				f4_1.write (str(ratioTest (p, k1, km1, k2, e0, s0, l)) + ' ')

		if (i == 3*N/10):
			for l in range (1, s0+1):
				f4_2.write (str(ratioTest (p, k1, km1, k2, e0, s0, l)) + ' ')

		if (i == 7*N/10):
			for l in range (1, s0+1):
				f4_3.write (str(ratioTest (p, k1, km1, k2, e0, s0, l)) + ' ')

	for l in range (1, s0+1):
			f5.write (str(l) + ' ')

	s_mean = []
	s_var = []

	for i in range (0,N+1):
		if (i % (N/(T_out-1)) == 0):
			s_mean.append (s_mean_dynamics[i])
			s_var.append (s_var_dynamics[i])

	interval = []
	# creating the time interval
	for i in range (0,T_out):
		interval.append(T*i/(T_out-1))
	
	str_s_mean = ''.join(str(e)+' ' for e in s_mean)
	str_s_var = ''.join(str(e)+' ' for e in s_var)

	f1.write(str_s_mean)
	f2.write(str_s_var)


if __name__ == '__main__':
	main()
