import os
import sys

def ensure_dir (f):
	if not os.path.exists(f):
		os.makedirs(f)

def get_constants (f):

	print "Please enter the reaction constants:"

	sys.stdout.write("k1 = ")
	f.write(raw_input() + '\n')

	sys.stdout.write("km1 = ")
	f.write(raw_input() + '\n')

	sys.stdout.write("k2 = ")
	f.write(raw_input() + '\n')


def get_initials (f):

	print "\nPlease enter the initial number of enzyme and substrate molecules:"

	sys.stdout.write("E0 = ")
	f.write(raw_input() + '\n')

	sys.stdout.write("S0 = ")
	f.write(raw_input() + '\n')


def get_time (f):

	print "\nPleaes enter the total time of observation (sec):"

	sys.stdout.write("Total Time = ")
	f.write(raw_input() + '\n')


def get_steps (f):

	print "\nPlease enter the number of Monte Carlo simulations:"

	sys.stdout.write("N = ")
	f.write(raw_input() + '\n')

def main():

	param = os.getcwd() + '/simulations/parameters'
	f = open (param, 'w')

	sys.stdout.write ("Please enter the number (1 - 3) of the figure that you want to reproduce: ")
	figure = int(raw_input())
	sys.stdout.write ('\n')


	ensure_dir ("figures")
	ensure_dir ("data")

	pyDir = os.getcwd() + "/simulations"

	if (figure == 1):
		get_constants (f)
		get_initials (f)
		get_time (f)
		get_steps (f)
		sys.stdout.write ('\n')
		f.close()

		os.system ("python " + pyDir + "/exact_solution.py")
		print "\nGenerating approximate solutions. Please wait..."
		os.system ("python " + pyDir + "/tqssa.py")

	if (figure == 2):
		get_constants (f)
		get_initials (f)
		get_time (f)
		sys.stdout.write ('\n')
		f.close()

		os.system ("python " + pyDir + "/exact_solution_ode.py")

	if (figure == 3):
		get_constants (f)
		get_initials (f)
		get_time (f)
		sys.stdout.write ('\n')
		f.close()

		os.system ("python " + pyDir + "/exact_solution_ode.py")
		os.system ("python " + pyDir + "/tqssa_barik.py")		


if __name__ == '__main__':
	main()