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

## sample size of Laqu.
N_Laqu = 78*2
## sample size of Lmeg.
N_Lmeg = 47*2
## sample size of Loua.
N_Loua = 10*2
## sample size of Lozk.
N_Lozk = 24*2
## sample size of Lpel.
N_Lpel = 29*2
## sample size of Lsol.
N_Lsol = 41*2

## sample size for all pops combined.
N_allpops = N_Laqu + N_Lmeg + N_Loua + N_Lozk + N_Lpel + N_Lsol

## create files to store parameters, trees and the models
os.mkdir("trainingSims")
par_2sp = open("par_2sp.txt","w")
par_6sp = open("par_6sp.txt","w")
par_6spMig = open("par_6spMig.txt","w")
trees_2sp = []
trees_6sp = []
trees_6spMig = []
Model_2sp = []
Model_6sp = []
Model_6spMig = []


####Simulate the species delimitation scenarios####
### Two species, Lpel and all other (note this hypothesis violates the topology recovered in Kim et al. 2022)
for i in range(Priorsize):

	### Define parameters
	Theta = random.uniform(0.001,0.02)

	## Sample Pi values between 0.5 and 1 (different species).
	Pab_inter=random.uniform(0.5,1)
	## obtain divergence time (tau - τ) priors using the Pi values. Pi = 1−2/3*e^(−2τ/θ); τ = ln((Pi-1)*-3/2)*-θ/2
	tau5 = math.log((Pab_inter-1)*-3/2)*-Theta/2
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

	## ms command
	com=subprocess.Popen("./ms %d 1000 -s 1 -t %f -I 6 %d %d %d %d %d %d -ej %f 4 2 -ej %f 2 3 -ej %f 3 6 -ej %f 6 1 -ej %f 5 1 -T" % (N_allpops, Theta, N_Laqu, N_Lmeg, N_Loua, N_Lozk, N_Lpel, N_Lsol, T1, T2, T3, T4, T5), shell=True, stdout=subprocess.PIPE).stdout
	output = com.read().splitlines()
	Model_2sp.append(np.array(ms2nparray(output)).swapaxes(0,1).reshape(N_allpops,-1).T)

	## save parameter values and trees
	par_2sp.write("%f\t%f\t%f\t%f\t%f\t%f\t%f\n" % (Pab_intra, Pab_inter,T1,T2,T3,T4,T5))
	#Randomly save a number of tree equivalent to the number of traits
	trees_2sp.append(random.sample(get_newick(output),28))
	print("Completed %d %% of Model 1 simulations" % (float(i)/Priorsize*100))

#Save the simulated SNP data
Model_2sp=np.array(Model_2sp)
np.savez_compressed('trainingSims/Model_2sp.npz', Model_2sp=Model_2sp)
del(Model_2sp)


### Six species
for i in range(Priorsize):

	### Define parameters
	Theta = random.uniform(0.001,0.02)

	## Sample Pi values between 0.5 and 1 (different species).
	Pab_inter5=random.uniform(0.5,1)
	Pab_inter4=random.uniform(0.5,Pab_inter5)
	Pab_inter3=random.uniform(0.5,Pab_inter4)
	Pab_inter2=random.uniform(0.5,Pab_inter3)
	Pab_inter=random.uniform(0.5,Pab_inter2)

	## obtain divergence time (tau - τ) priors using the Pi values. Pi = 1−2/3*e^(−2τ/θ); τ = ln((Pi-1)*-3/2)*-θ/2
	tau5 = math.log((Pab_inter5-1)*-3/2)*-Theta/2
	tau4 = math.log((Pab_inter4-1)*-3/2)*-Theta/2
	tau3 = math.log((Pab_inter3-1)*-3/2)*-Theta/2
	tau2 = math.log((Pab_inter2-1)*-3/2)*-Theta/2
	tau1 = math.log((Pab_inter-1)*-3/2)*-Theta/2
	## Transform divergence times to 4Ne generations units (required by ms) 
	T5 = 4*tau5/Theta
	T4 = 4*tau4/Theta
	T3 = 4*tau3/Theta
	T2 = 4*tau2/Theta
	T1 = 4*tau1/Theta


	## Set the Pi values for intraspecific nodes to a default value (0)
	Pab_intra=0

	## ms command
	com=subprocess.Popen("./ms %d 1000 -s 1 -t %f -I 6 %d %d %d %d %d %d -ej %f 5 4 -ej %f 4 2 -ej %f 2 3 -ej %f 3 6 -ej %f 6 1 -T" % (N_allpops, Theta, N_Laqu, N_Lmeg, N_Loua, N_Lozk, N_Lpel, N_Lsol, T1, T2, T3, T4, T5), shell=True, stdout=subprocess.PIPE).stdout
	output = com.read().splitlines()
	Model_6sp.append(np.array(ms2nparray(output)).swapaxes(0,1).reshape(N_allpops,-1).T)

	## save parameter values and trees
	par_6sp.write("%f\t%f\t%f\t%f\t%f\t%f\t%f\n" % (Pab_intra, Pab_inter,T1,T2,T3,T4,T5))
	#Randomly save a number of tree equivalent to the number of traits
	trees_6sp.append(random.sample(get_newick(output),28))
	print("Completed %d %% of Model 2 simulations" % (float(i)/Priorsize*100))

#Save the simulated SNP data
Model_6sp=np.array(Model_6sp)
np.savez_compressed('trainingSims/Model_6sp.npz', Model_6sp=Model_6sp)
del(Model_6sp)

### Six species + M
for i in range(Priorsize):

	### Define parameters
	Theta = random.uniform(0.001,0.02)

	## Sample Pi values between 0.5 and 1 (different species).
	Pab_inter5=random.uniform(0.5,1)
	Pab_inter4=random.uniform(0.5,Pab_inter5)
	Pab_inter3=random.uniform(0.5,Pab_inter4)
	Pab_inter2=random.uniform(0.5,Pab_inter3)
	Pab_inter=random.uniform(0.5,Pab_inter2)

	## obtain divergence time (tau - τ) priors using the Pi values. Pi = 1−2/3*e^(−2τ/θ); τ = ln((Pi-1)*-3/2)*-θ/2
	tau5 = math.log((Pab_inter5-1)*-3/2)*-Theta/2
	tau4 = math.log((Pab_inter4-1)*-3/2)*-Theta/2
	tau3 = math.log((Pab_inter3-1)*-3/2)*-Theta/2
	tau2 = math.log((Pab_inter2-1)*-3/2)*-Theta/2
	tau1 = math.log((Pab_inter-1)*-3/2)*-Theta/2
	## Transform divergence times to 4Ne generations units (required by ms) 
	T5 = 4*tau5/Theta
	T4 = 4*tau4/Theta
	T3 = 4*tau3/Theta
	T2 = 4*tau2/Theta
	T1 = 4*tau1/Theta
	
	## Sample the migration rate
	M=random.uniform(0,1)


	## Set the Pi values for intraspecific nodes to a default value (0)
	Pab_intra=0

	## ms command
	com=subprocess.Popen("./ms %d 1000 -s 1 -t %f -I 6 %d %d %d %d %d %d -m 1 3 %f -m 1 5 %f -m 2 6 %f -m 4 1 %f -m 4 2 %f -m 4 2 %f -m 4 2 %f -m 6 2 %f -ej %f 5 4 -ej %f 4 2 -ej %f 2 3 -ej %f 3 6 -ej %f 6 1 -T" % (N_allpops, Theta, N_Laqu, N_Lmeg, N_Loua, N_Lozk, N_Lpel, N_Lsol, M, M, M, M, M, M, M, M, T1, T2, T3, T4, T5), shell=True, stdout=subprocess.PIPE).stdout
	output = com.read().splitlines()
	Model_6spMig.append(np.array(ms2nparray(output)).swapaxes(0,1).reshape(N_allpops,-1).T)

	## save parameter values and trees
	par_6spMig.write("%f\t%f\t%f\t%f\t%f\t%f\t%f\n" % (Pab_intra, Pab_inter,T1,T2,T3,T4,T5))
	#Randomly save a number of tree equivalent to the number of traits
	trees_6spMig.append(random.sample(get_newick(output),28))
	print("Completed %d %% of Model 3 simulations" % (float(i)/Priorsize*100))

#Save the simulated SNP data
Model_6spMig=np.array(Model_6spMig)
np.savez_compressed('trainingSims/Model_6spMig.npz', Model_6spMig=Model_6spMig)
del(Model_6spMig)

#Save trees from all scenarios
trees=np.concatenate((trees_2sp,trees_6sp,trees_6spMig),axis=0)
np.savez_compressed('trees.npz', trees=trees)