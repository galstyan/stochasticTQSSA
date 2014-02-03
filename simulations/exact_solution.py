import sys
import random
import os
import math

inf = 100000

T_out = 25 # number of time points in output
T = 10000 # number of Monte Carlo steps during each realization

# propensity functions
def prop1 (e,s, k1, km1, k2, e0, s0):
	return k1 * e * s
	
def propm1 (e,s, k1, km1, k2, e0, s0):
	return km1 * (e0-e)

def prop2 (e,s, k1, km1, k2, e0, s0):
	return k2 * (e0-e)
	
def sum_prop (e,s, k1, km1, k2, e0, s0):
	sum = prop1(e,s, k1, km1, k2, e0, s0) + propm1(e,s, k1, km1, k2, e0, s0) + prop2(e,s, k1, km1, k2, e0, s0)
	return sum

# rate functions normalized to 1
def rate1 (e,s, k1, km1, k2, e0, s0):
	rate = prop1(e,s, k1, km1, k2, e0, s0)/sum_prop(e,s, k1, km1, k2, e0, s0)
	return rate
	
def ratem1 (e,s, k1, km1, k2, e0, s0):
	rate = propm1(e,s, k1, km1, k2, e0, s0)/sum_prop(e,s, k1, km1, k2, e0, s0)
	return rate
	
def rate2 (e,s, k1, km1, k2, e0, s0):
	rate = prop2(e,s, k1, km1, k2, e0, s0)/sum_prop(e,s, k1, km1, k2, e0, s0)
	return rate
	
# change in number of molecules for the given reaction
def e_step (e,s,r, k1, km1, k2, e0, s0):
	if r < ratem1(e,s, k1, km1, k2, e0, s0):
		return 1
	if r < ratem1(e,s, k1, km1, k2, e0, s0) + rate1(e,s, k1, km1, k2, e0, s0):
		return -1
	return 1
	
def s_step (e,s,r, k1, km1, k2, e0, s0):
	if r < ratem1(e,s, k1, km1, k2, e0, s0):
		return 1
	if r < ratem1(e,s, k1, km1, k2, e0, s0) + rate1(e,s, k1, km1, k2, e0, s0):
		return -1
	return 0
	
# mean number of molecules at time t over N sample path realizations
def mean (matrix, t, N):
	mean = 0.0;
	for n in range (0,N):
		mean = mean + matrix[n][t]
	mean = mean/N;
	return mean

# variance of molecules at time t over N sample path realizations	
def var (matrix, t, N):
	var = 0.0
	for n in range (0,N):
		var = var + matrix[n][t] * matrix[n][t]
	var = var / N
	var = var - mean(matrix,t, N) * mean(matrix,t, N)
	
	return math.sqrt(var)
	
def main():

	fin = open (os.path.dirname(__file__) + "/parameters", 'r')

	k1 = float(fin.readline())
	km1 = float(fin.readline())
	k2 = float(fin.readline())

	e0 = int(fin.readline())
	s0 = int(fin.readline())

	time = float(fin.readline())

	N = int(fin.readline())

	cwd = os.path.dirname( os.path.dirname(__file__) )
	f1_name = cwd + "/data/Interval"
	f2_name = cwd + "/data/Exact_Substrate_Mean"
	f3_name = cwd + "/data/Exact_Substrate_SD"
	f4_name = cwd + "/data/Exact_Substrate_Mean_PSD"
	f5_name = cwd + "/data/Exact_Substrate_Mean_MSD"
	
	f1 = open (f1_name, 'w')
	f2 = open (f2_name, 'w')
	f3 = open (f3_name, 'w')
	f4 = open (f4_name, 'w')
	f5 = open (f5_name, 'w')
	
	time_all = [[] for i in range(N)]
	e_all = [[] for i in range(N)]
	s_all = [[] for i in range(N)]
	
	s_mean_dynamics = []
	s_var_dynamics = []	
	s_mean_psd_dynamics = []
	s_mean_msd_dynamics = []
	
	time_max = 0

	# N sample paths in Monte Carlo
	for n in range (0, N):
		
		# initialize molecule number variables
		e = e0
		s = s0
		
		time_all[n].append(0)
		e_all[n].append(e)
		s_all[n].append(s)		

		for t in range (0, T):
			if (sum_prop(e,s, k1, km1, k2, e0, s0) == 0): # reaction has ended
				
				time_all[n].extend ([inf]*(T-1-t))
				e_all[n].extend ([e0]*(T-1-t))
				s_all[n].extend ([0]*(T-1-t))
				
				break
			else:
				
				r = random.random()
				del_e = e_step(e,s,r, k1, km1, k2, e0, s0)
				del_s = s_step(e,s,r, k1, km1, k2, e0, s0)
				
				time_all[n].append(time_all[n][t]+1.0 / sum_prop(e,s, k1, km1, k2, e0, s0))
				
				if (time_all[n][t]+1.0 / sum_prop(e,s, k1, km1, k2, e0, s0) > time_max):
					time_max = time_all[n][t]+1.0 / sum_prop(e,s, k1, km1, k2, e0, s0)
				
				e = e + del_e
				s = s + del_s
				
				e_all[n].append(e)
				s_all[n].append(s)		
		
		# status report
		if (10*(n+1) % N == 0):
			print 'Monte Carlo Simulation ... ', 100.0*(n+1)/N, '%'
		
	s_data = [[] for i in range(N)]
	c_data = [[] for i in range(N)]
	
	interval = []
	
	# creating the time interval based on the maximum reaction duration
	for i in range (0,T_out):
		interval.append(time*i/(T_out-1))
	
	
	# fitting the data to the time scale
	# because time moments for different Monte Carlo realizations vary
	for n in range (0,N):
		j = 0
		i = 0
		while j < T_out:
			if time_all[n][i] <= interval[j]:
				i = i+1
			else:
				s_data[n].append( s_all[n][i-1])
				c_data[n].append (e0 - e_all[n][i-1])
				j = j+1
				
	# creating lists for the dynamics of the parameters of interest
	for t in range (0,T_out):
		s_mean_dynamics.append( mean(s_data, t, N) + mean(c_data,t,N))
		s_var_dynamics.append ( var(s_data,t, N))
		s_mean_psd_dynamics.append (mean(s_data, t, N) + var(s_data, t, N))
		s_mean_msd_dynamics.append (mean(s_data, t, N) - var(s_data, t, N))
	
	# converting into strings for output
	strr = ''.join(str(e)+' ' for e in interval)	
	str_s_mean = ''.join(str(e) + ' ' for e in s_mean_dynamics)
	str_s_var = ''.join(str(e) + ' ' for e in s_var_dynamics)
	str_s_mean_psd = ''.join(str(e) + ' ' for e in s_mean_psd_dynamics)
	str_s_mean_msd = ''.join(str(e) + ' ' for e in s_mean_msd_dynamics)

	f1.write(strr)
	f2.write(str_s_mean)
	f3.write(str_s_var)
	f4.write(str_s_mean_psd)
	f5.write(str_s_mean_msd)
	
	
if __name__ == '__main__':
	main()
