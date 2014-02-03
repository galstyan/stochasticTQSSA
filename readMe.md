This package contains python and Matlab files that simulate the stochastic Michaelis-Menten kinetics and compare it with the results from simplified approximate models.

The Michaelis-Menten reaction is given by the following scheme:

E + S -> C (rate k1)
C -> E + S (rate km1)
C -> E + P (rate k2)

1. Run 'launch.py' file

2. Choose the figure you want to reproduce

3. Enter the reaction parameters

   (The launch program will run the simulator programs in '\simulations' folder and output the simualtion data into '\data' folder)

4. Go to '\plot_generators' folder and run the matlab program corresponding to the figure you chose initially

5. The figure in .EPS format will be created in '\figures' folder