#!/usr/bin/python3

## in order to use this code you have to have ms installed on your computer
## ms can be freely downloaded from:
## http://home.uchicago.edu/rhudson1/source/mksamples.html

import random
import os
import math
import shlex, subprocess
import numpy as np

##define a function to read ms' simulations and transform then into a NumPy array.    
def ms2nparray(xfile):
	g = list(xfile)
	k = [idx for idx,i in enumerate(g) if len(i) > 0 and i.startswith(b'//')]
	f = []
	for i in k:
		L = g[i+5:i+N_allpops+5]
		q = []
		for i in L:
			#originally: i = [int(j) for j in list(i)], need to decode as python3 loads as ASCII (0->48; 1->49).
			i = [int(j) for j in list(i.decode('utf-8'))]
			i = np.array(i, dtype=np.int8)
			q.append(i)
		q = np.array(q)
		q = q.astype("int8")
		f.append(np.array(q))
	return f

##define a function to read ms' coalescent trees simulations and add them to a list.
def get_newick(xfile):
	g = list(xfile)
	k = [idx for idx,i in enumerate(g) if len(i) > 0 and i.startswith(b'//')]
	t = []
	for i in k:
		n = g[i+1]
		t.append(n.decode('utf-8'))
	return t


### variable declarations

#define the number of simulations
Priorsize = 10000

## sample size of popA.
N_popA = 10
## sample size of popB.
N_popB = 10
## sample size of popC.
N_popC = 10
## sample size of popD.
N_popD = 10
## sample size of popE.
N_popE= 10
## sample size of popF.
N_popF = 10

## sample size for all pops combined.
N_allpops = N_popA + N_popB + N_popC + N_popD + N_popE + N_popF

## create files to store parameters, trees and the models
os.mkdir("trainingSims")
par_1sp = open("par_1sp.txt","w")
par_2sp = open("par_2sp.txt","w")
par_3sp = open("par_3sp.txt","w")
trees_1sp = []
trees_2sp = []
trees_3sp = []
Model_1sp = []
Model_2sp = []
Model_3sp = []

####Simulate the species delimitation scenarios####
### One species
for i in range(Priorsize):
	### Define parameters
	## Theta value of 0.005
	Theta = 0.005

	## Sample Pi values between 1/3 (minimum value -> 0.333) and 0.4 (same species).
	Pab_intra=random.uniform(1/3,0.4)
	## obtain divergence time (tau - τ) priors using the Pi values. Pi = 1−2/3*e^(−2τ/θ); τ = ln((Pi-1)*-3/2)*-θ/2
	tau5 = math.log((Pab_intra-1)*-3/2)*-Theta/2
	## Transform divergence times to 4Ne generations units (required by ms)
	T5 = 4*tau5/Theta
	## Sample the more recent splitting times with a smaller value than the more ancient one
	T4 = random.uniform(0,T5)
	T3 = random.uniform(0,T4)
	T2 = random.uniform(0,T3)
	T1 = random.uniform(0,T2)

	## Set the Pi values for interpecific nodes to a default value (0)
	Pab_inter1 = 0
	Pab_inter2 = 0

	## ms command
	com=subprocess.Popen("./ms %d 1000 -s 1 -t %f -I 6 %d %d %d %d %d %d -ej %f 2 1 -ej %f 4 3 -ej %f 6 5 -ej %f 5 3 -ej %f 3 1 -T" % (N_allpops, Theta, N_popA, N_popB, N_popC, N_popD, N_popE, N_popF, T1, T2, T3, T4, T5), shell=True, stdout=subprocess.PIPE).stdout
	output = com.read().splitlines()
	Model_1sp.append(np.array(ms2nparray(output)).swapaxes(0,1).reshape(N_allpops,-1).T)

	## save parameter values and trees
	par_1sp.write("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n" % (Pab_intra, Pab_inter1, Pab_inter2,T1,T2,T3,T4,T5))
	#Randomly save a number of tree equivalent to the number of traits
	trees_1sp.append(random.sample(get_newick(output),100))
	print("Completed %d %% of Model 1 simulations" % (float(i)/Priorsize*100))

#Save the simulated SNP data
Model_1sp=np.array(Model_1sp)
np.savez_compressed('trainingSims/Model_1sp.npz', Model_1sp=Model_1sp)
del(Model_1sp)

### Two species, AB and CDEF
for i in range(Priorsize):

	## Theta value of 0.005
	Theta = 0.005

	## Sample Pi values between 0.5 and 1 (different species).
	Pab_inter1=random.uniform(0.5,1)
	## obtain divergence time (tau - τ) priors using the Pi values. Pi = 1−2/3*e^(−2τ/θ); τ = ln((Pi-1)*-3/2)*-θ/2
	tau5 = math.log((Pab_inter1-1)*-3/2)*-Theta/2
	## Transform divergence times to 4Ne generations units (required by ms)
	T5 = 4*tau5/Theta

	## Sample Pi values between 1/3 (minimum value -> 0.333) and 0.4 (same species).
	Pab_intra=random.uniform(1/3,0.4)
	## obtain divergence time (tau - τ) priors using the Pi values. Pi = 1−2/3*e^(−2τ/θ); τ = ln((Pi-1)*-3/2)*-θ/2
	tau4 = math.log((Pab_intra-1)*-3/2)*-Theta/2
	## Transform divergence times to 4Ne generations units (required by ms)
	T4 = 4*tau4/Theta
	## Sample the more recent splitting times with a smaller value than the more ancient one
	T3 = random.uniform(0,T4)
	T2 = random.uniform(0,T3)
	T1 = random.uniform(0,T2)

	## Set the Pi values for the second interspecific node to a default value (0)
	Pab_inter2 = 0

	## ms command
	com=subprocess.Popen("./ms %d 1000 -s 1 -t %f -I 6 %d %d %d %d %d %d -ej %f 2 1 -ej %f 4 3 -ej %f 6 5 -ej %f 5 3 -ej %f 3 1 -T" % (N_allpops, Theta, N_popA, N_popB, N_popC, N_popD, N_popE, N_popF, T1, T2, T3, T4, T5), shell=True, stdout=subprocess.PIPE).stdout
	output = com.read().splitlines()
	Model_2sp.append(np.array(ms2nparray(output)).swapaxes(0,1).reshape(N_allpops,-1).T)

	## save parameter values and models
	par_2sp.write("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n" % (Pab_intra, Pab_inter1, Pab_inter2,T1,T2,T3,T4,T5))
	#Randomly save a number of tree equivalent to the number of traits
	trees_2sp.append(random.sample(get_newick(output),100))
	print("Completed %d %% of Model 2 simulations" % (float(i)/Priorsize*100))

#Save the simulated SNP data
Model_2sp=np.array(Model_2sp)
np.savez_compressed('trainingSims/Model_2sp.npz', Model_2sp=Model_2sp)
majority_trees = []
del(Model_2sp)

### Three species
for i in range(Priorsize):

	## Theta value of 0.005
	Theta = 0.005

	## Sample Pi values between 0.5 and 1 (same species).
	Pab_inter1=random.uniform(0.5,1)
	Pab_inter2=random.uniform(0.5,Pab_inter1)
	## obtain divergence time (tau - τ) priors using the Pi values. Pi = 1−2/3*e^(−2τ/θ); τ = ln((Pi-1)*-3/2)*-θ/2
	tau5 = math.log((Pab_inter1-1)*-3/2)*-Theta/2
	tau4 = math.log((Pab_inter2-1)*-3/2)*-Theta/2
	## Transform divergence times to 4Ne generations units (required by ms)
	T5 = 4*tau5/Theta
	T4 = 4*tau4/Theta

	## Sample Pi values between 1/3 (minimum value -> 0.333) and 0.4 (same species).
	Pab_intra=random.uniform(1/3,0.4)
	## obtain divergence time (tau - τ) priors using the Pi values. Pi = 1−2/3*e^(−2τ/θ); τ = ln((Pi-1)*-3/2)*-θ/2
	tau3 = math.log((Pab_intra-1)*-3/2)*-Theta/2
	## Transform divergence times to 4Ne generations units (required by ms)
	T3 = 4*tau3/Theta
	## Sample the more recent splitting times with a smaller value than the more ancient one
	T2 = random.uniform(0,T3)
	T1 = random.uniform(0,T2)

	## ms command
	com=subprocess.Popen("./ms %d 1000 -s 1 -t %f -I 6 %d %d %d %d %d %d -ej %f 2 1 -ej %f 4 3 -ej %f 6 5 -ej %f 5 3 -ej %f 3 1 -T" % (N_allpops, Theta, N_popA, N_popB, N_popC, N_popD, N_popE, N_popF, T1, T2, T3, T4, T5), shell=True, stdout=subprocess.PIPE).stdout
	output = com.read().splitlines()
	Model_3sp.append(np.array(ms2nparray(output)).swapaxes(0,1).reshape(N_allpops,-1).T)
	## save parameter values
	par_3sp.write("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n" % (Pab_intra, Pab_inter1, Pab_inter2,T1,T2,T3,T4,T5))
	#Randomly save a number of tree equivalent to the number of traits
	trees_3sp.append(random.sample(get_newick(output),100))
	print("Completed %d %% of Model 3 simulations" % (float(i)/Priorsize*100))

#Save the simulated SNP data
Model_3sp=np.array(Model_3sp)
np.savez_compressed('trainingSims/Model_3sp.npz', Model_3sp=Model_3sp)
#Save trees from all scenarios
trees=np.concatenate((trees_1sp,trees_2sp,trees_3sp),axis=0)
np.savez_compressed('trees.npz', trees=trees)
