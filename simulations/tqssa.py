import sys
import math
import os

T_out = 25 # number of time points in output
N = 1000*(T_out-1) # number of steps in Euler method

def prop(s, k1, km1, k2, Km, e0, s0):
	if (s == s0+1):
		return 0
	r = (e0+Km+s)
	return 2*k2*e0*s/(r + math.sqrt(r*r-4*e0*s))
	
# mean number of molecules at time t over N sample path realizations

def mean_ss (data, k1, km1, k2, Km, e0, s0):
	mean = 0.0
	for s in range (0, s0+1):
		mean = mean + (s-complex(s, k1, km1, k2, Km, e0, s0))*data[s]
	return mean

def complex (s, k1, km1, k2, Km, e0, s0):
	r = e0 + Km + s
	r = r * 1.0
	return 2*e0*s/(r + math.sqrt(r*r-4*e0*s))

# variance of molecules at time t over N sample path realizations	

def var_ss (data, k1, km1, k2, Km, e0, s0):
	var = 0.0
	m = mean_ss(data, k1, km1, k2, Km, e0, s0)
	for s in range (0,s0+1):
		var = var + (s-m-complex(s, k1, km1, k2, Km, e0, s0))*(s-m-complex(s, k1, km1, k2, Km, e0, s0))*data[s]
	
	return math.sqrt(var)	

	
def main():

	fin = open (os.path.dirname(__file__) + "/parameters", 'r')

	k1 = float(fin.readline())
	km1 = float(fin.readline())
	k2 = float(fin.readline())

	e0 = int(fin.readline())
	s0 = int(fin.readline())

	T = float(fin.readline())

	Km = (km1+k2)/k1

	dt = T/N


	cwd = os.path.dirname( os.path.dirname(__file__) )
	f1_name = cwd + "/data/TQSSA_Substrate_Mean"
	f2_name = cwd + "/data/TQSSA_Substrate_SD"
	f3_name = cwd + "/data/TQSSA_Substrate_Mean_PSD"
	f4_name = cwd + "/data/TQSSA_Substrate_Mean_MSD"


	f1 = open (f1_name, 'w')
	f2 = open (f2_name, 'w')
	f3 = open (f3_name, 'w')
	f4 = open (f4_name, 'w')

	# initializing the probability distribution
	p = [0]*(s0+2)
	p[s0] = 1
	
	q = [0]*(s0+2)
	q[s0] = 1
	

	ss_mean_dynamics = []
	ss_mean_dynamics.append(s0)

	ss_var_dynamics = []	
	ss_var_dynamics.append(0)

	ss_mean_psd_dynamics = []
	ss_mean_psd_dynamics.append(s0)

	ss_mean_msd_dynamics = []
	ss_mean_msd_dynamics.append(s0)
	
	for i in range (1,N+1):

		for s in range (0,s0+1):
			p[s] = q[s] + dt * (prop(s+1, k1, km1, k2, Km, e0, s0)*q[s+1] - prop(s, k1, km1, k2, Km, e0, s0)*q[s])

		ss_mean_dynamics.append(mean_ss(p, k1, km1, k2, Km, e0, s0))
		ss_var_dynamics.append(var_ss(p, k1, km1, k2, Km, e0, s0))

		ss_mean_psd_dynamics.append(mean_ss(p, k1, km1, k2, Km, e0, s0) + var_ss(p, k1, km1, k2, Km, e0, s0))
		ss_mean_msd_dynamics.append(mean_ss(p, k1, km1, k2, Km, e0, s0) - var_ss(p, k1, km1, k2, Km, e0, s0))
		
		del q[:]
		q.extend(p)
		
		sum = 0
		for s in range (0,s0+1):
			sum = sum + p[s]
	


	ss_mean = []
	ss_var = []
	ss_mean_psd = []
	ss_mean_msd = []

	
	for i in range (0,N+1):
		if (i % (N/(T_out-1)) == 0):

			ss_mean.append (ss_mean_dynamics[i])
			ss_var.append (ss_var_dynamics[i])
			ss_mean_psd.append (ss_mean_psd_dynamics[i])
			ss_mean_msd.append (ss_mean_msd_dynamics[i])
			
			
			
	interval = []
	# creating the time interval
	for i in range (0,T_out):
		interval.append(T*i/(T_out-1))
	

	str_ss_mean = ''.join(str(e)+' ' for e in ss_mean)
	str_ss_var = ''.join(str(e)+' ' for e in ss_var)
	str_ss_mean_psd = ''.join(str(e)+' ' for e in ss_mean_psd)
	str_ss_mean_msd = ''.join(str(e)+' ' for e in ss_mean_msd)


	f1.write(str_ss_mean)
	f2.write(str_ss_var)
	f3.write(str_ss_mean_psd)
	f4.write(str_ss_mean_msd)

if __name__ == '__main__':
	main()
