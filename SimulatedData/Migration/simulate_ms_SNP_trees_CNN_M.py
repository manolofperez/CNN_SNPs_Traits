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
par_2sp_01 = open("par_2sp_01.txt","w")
par_3sp_01 = open("par_3sp_01.txt","w")
par_2sp_05 = open("par_2sp_05.txt","w")
par_3sp_05 = open("par_3sp_05.txt","w")
par_2sp_1 = open("par_2sp_1.txt","w")
par_3sp_1 = open("par_3sp_1.txt","w")
trees_2sp_01 = []
trees_3sp_01 = []
trees_2sp_05 = []
trees_3sp_05 = []
trees_2sp_1 = []
trees_3sp_1 = []
Model_2sp_01 = []
Model_3sp_01 = []
Model_2sp_05 = []
Model_3sp_05 = []
Model_2sp_1 = []
Model_3sp_1 = []

###Migration 0.1
## Two species, AB and CDEF
for i in range(Priorsize):

	## Theta value of 0.005
	Theta = 0.005

	## Sample Pi values for the deepest node (Pi1) between 0.5 and 1 (different species).
	Pi1=random.uniform(0.5,1)
	## obtain divergence time (tau - τ) priors using the Pi values. Pi = 1−2/3*e^(−2τ/θ); τ = -(θ*ln(-(3*Pi-3)/2)/2)
	tau1 = -(Theta*math.log(-(3*Pi1-3)/2)/2)  
	## Transform divergence times to 4Ne generations units (required by ms)
	T1 = 4*tau5/Theta

	## Sample Pi values for the second deepest node (Pi2) between 1/3 (minimum value -> 0.333) and 0.4 (same species).
	Pi2=random.uniform(1/3,0.4)
	## obtain divergence time (tau - τ) priors using the Pi values. Pi = 1−2/3*e^(−2τ/θ); τ = -(θ*ln(-(3*Pi-3)/2)/2)
	tau2 = -(Theta*math.log(-(3*Pi2-3)/2)/2)  
	## Transform divergence times to 4Ne generations units (required by ms)
	T2 = 4*tau2/Theta

	## Sample the more recent splitting times with a smaller value than the more ancient ones and bigger than 0 (present)
	T3 = random.uniform(0,T2)
	T4 = random.uniform(0,T3)
	T5 = random.uniform(0,T4)

	## ms command
	com=subprocess.Popen("./ms %d 1000 -s 1 -t %f -I 6 %d %d %d %d %d %d -ej %f 2 1 -ej %f 4 3 -ej %f 6 5 -ej %f 5 3 -eM %f 0.1 -ej %f 3 1 -T" % (N_allpops, Theta, N_popA, N_popB, N_popC, N_popD, N_popE, N_popF, T5, T4, T3, T2, T2, T1), shell=True, stdout=subprocess.PIPE).stdout
	output = com.read().splitlines()
	Model_2sp_01.append(np.array(ms2nparray(output)).swapaxes(0,1).reshape(N_allpops,-1).T)

	## save parameter values and trees
	par_2sp_01.write("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n" % (Pi1, Pi2,T1,T2,T3,T4,T5))
	#Randomly save a number of tree equivalent to the number of traits
	trees_2sp_01.append(random.sample(get_newick(output),100))
	print("Completed %d %% of Model 2sp_01 simulations" % (float(i)/Priorsize*100))

#Save the simulated SNP data
Model_2sp_01=np.array(Model_2sp_01)
np.savez_compressed('trainingSims/Model_2sp_01.npz', Model_2sp_01=Model_2sp_01)
del(Model_2sp_01)


## Three species
for i in range(Priorsize):

	## Theta value of 0.005
	Theta = 0.005

	## Sample Pi values for the deepest node (Pi1) between 0.5 and 1 (different species).
	Pi1 = random.uniform(0.5,1)
	## obtain divergence time (tau - τ) priors using the Pi values. Pi = 1−2/3*e^(−2τ/θ); τ = -(θ*ln(-(3*Pi-3)/2)/2)
	tau1 = -(Theta*math.log(-(3*Pi1-3)/2)/2)  
	## Transform divergence times to 4Ne generations units (required by ms)
	T1 = 4*tau1/Theta

	## Sample Pi values for the second deepest node (Pi2) that are higer than 0.5 (different species) but lower than the deepest node (Pi1).
	Pi2 = random.uniform(0.5,Pi1)
	## obtain divergence time (tau - τ) priors using the Pi values. Pi = 1−2/3*e^(−2τ/θ); τ = -(θ*ln(-(3*Pi-3)/2)/2)
	tau2 = -(Theta*math.log(-(3*Pi2-3)/2)/2)  
	## Transform divergence times to 4Ne generations units (required by ms)
	T2 = 4*tau2/Theta

	## Sample Pi values for the oldest intraspecific node (Pi3) between 1/3 (minimum value -> 0.333) and 0.4 (same species).
	Pi3=random.uniform(1/3,0.4)
	## obtain divergence time (tau - τ) priors using the Pi values. Pi = 1−2/3*e^(−2τ/θ); τ = -(θ*ln(-(3*Pi-3)/2)/2)
	tau3 = -(Theta*math.log(-(3*Pi3-3)/2)/2)  
	## Transform divergence times to 4Ne generations units (required by ms)
	T3 = 4*tau3/Theta

	## Sample the more recent splitting times with a smaller value than the more ancient ones and bigger than 0 (present)
	T4 = random.uniform(0,T4)
	T5 = random.uniform(0,T5)

	## ms command
	com=subprocess.Popen("./ms %d 1000 -s 1 -t %f -I 6 %d %d %d %d %d %d -ej %f 2 1 -ej %f 4 3 -ej %f 6 5 -eM %f 0.1 -ej %f 5 3 -eM %f 0.1 -ej %f 3 1 -T" % (N_allpops, Theta, N_popA, N_popB, N_popC, N_popD, N_popE, N_popF, T5, T4, T3, T3, T2, T2, T1), shell=True, stdout=subprocess.PIPE).stdout
	output = com.read().splitlines()
	Model_3sp_01.append(np.array(ms2nparray(output)).swapaxes(0,1).reshape(N_allpops,-1).T)
	
	## save parameter values and trees
	par_3sp_01.write("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n" % (Pi1, Pi2,T1,T2,T3,T4,T5))
	#Randomly save a number of tree equivalent to the number of traits
	trees_3sp_01.append(random.sample(get_newick(output),100))
	print("Completed %d %% of Model 3sp_01 simulations" % (float(i)/Priorsize*100))


#Save the simulated SNP data
Model_3sp_01=np.array(Model_3sp_01)
np.savez_compressed('trainingSims/Model_3sp_01.npz', Model_3sp_01=Model_3sp_01)
del(Model_3sp_01)

###Migration 0.5
## Two species, AB and CDEF
for i in range(Priorsize):

	## Theta value of 0.005
	Theta = 0.005

	## Sample Pi values for the deepest node (Pi1) between 0.5 and 1 (different species).
	Pi1 = random.uniform(0.5,1)
	## obtain divergence time (tau - τ) priors using the Pi values. Pi = 1−2/3*e^(−2τ/θ); τ = -(θ*ln(-(3*Pi-3)/2)/2)
	tau1 = -(Theta*math.log(-(3*Pi1-3)/2)/2)  
	## Transform divergence times to 4Ne generations units (required by ms)
	T1 = 4*tau1/Theta

	## Sample Pi values for the second deepest node (Pi2) between 1/3 (minimum value -> 0.333) and 0.4 (same species).
	Pi2 = random.uniform(1/3,0.4)
	## obtain divergence time (tau - τ) priors using the Pi values. Pi = 1−2/3*e^(−2τ/θ); τ = -(θ*ln(-(3*Pi-3)/2)/2)
	tau2 = -(Theta*math.log(-(3*Pi2-3)/2)/2)  
	## Transform divergence times to 4Ne generations units (required by ms)
	T2 = 4*tau2/Theta

	## Sample the more recent splitting times with a smaller value than the more ancient ones and bigger than 0 (present)
	T3 = random.uniform(0,T2)
	T4 = random.uniform(0,T3)
	T5 = random.uniform(0,T4)

	## ms command
	com=subprocess.Popen("./ms %d 1000 -s 1 -t %f -I 6 %d %d %d %d %d %d -ej %f 2 1 -ej %f 4 3 -ej %f 6 5 -ej %f 5 3 -eM %f 0.5 -ej %f 3 1 -T" % (N_allpops, Theta, N_popA, N_popB, N_popC, N_popD, N_popE, N_popF, T5, T4, T3, T2, T2, T1), shell=True, stdout=subprocess.PIPE).stdout
	output = com.read().splitlines()
	Model_2sp_05.append(np.array(ms2nparray(output)).swapaxes(0,1).reshape(N_allpops,-1).T)

	## save parameter values and trees
	par_2sp_05.write("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n" % (Pi1, Pi2,T1,T2,T3,T4,T5))
	#Randomly save a number of tree equivalent to the number of traits
	trees_2sp_05.append(random.sample(get_newick(output),100))
	print("Completed %d %% of Model 2sp_05 simulations" % (float(i)/Priorsize*100))

#Save the simulated SNP data
Model_2sp_05=np.array(Model_2sp_05)
np.savez_compressed('trainingSims/Model_2sp_05.npz', Model_2sp_05=Model_2sp_05)
del(Model_2sp_05)


## Three species
for i in range(Priorsize):

	## Theta value of 0.005
	Theta = 0.005

	## Sample Pi values for the deepest node (Pi1) between 0.5 and 1 (different species).
	Pi1 = random.uniform(0.5,1)
	## obtain divergence time (tau - τ) priors using the Pi values. Pi = 1−2/3*e^(−2τ/θ); τ = -(θ*ln(-(3*Pi-3)/2)/2)
	tau1 = -(Theta*math.log(-(3*Pi1-3)/2)/2)  
	## Transform divergence times to 4Ne generations units (required by ms)
	T1 = 4*tau1/Theta

	## Sample Pi values for the second deepest node (Pi2) that are higer than 0.5 (different species) but lower than the deepest node (Pi1).
	Pi2 = random.uniform(0.5,Pi1)
	## obtain divergence time (tau - τ) priors using the Pi values. Pi = 1−2/3*e^(−2τ/θ); τ = -(θ*ln(-(3*Pi-3)/2)/2)
	tau2 = -(Theta*math.log(-(3*Pi2-3)/2)/2)  
	## Transform divergence times to 4Ne generations units (required by ms)
	T2 = 4*tau2/Theta

	## Sample Pi values for the oldest intraspecific node (Pi3) between 1/3 (minimum value -> 0.333) and 0.4 (same species).
	Pi3=random.uniform(1/3,0.4)
	## obtain divergence time (tau - τ) priors using the Pi values. Pi = 1−2/3*e^(−2τ/θ); τ = -(θ*ln(-(3*Pi-3)/2)/2)
	tau3 = -(Theta*math.log(-(3*Pi3-3)/2)/2)  
	## Transform divergence times to 4Ne generations units (required by ms)
	T3 = 4*tau3/Theta

	## Sample the more recent splitting times with a smaller value than the more ancient ones and bigger than 0 (present)
	T4 = random.uniform(0,T4)
	T5 = random.uniform(0,T5)

	## ms command
	com=subprocess.Popen("./ms %d 1000 -s 1 -t %f -I 6 %d %d %d %d %d %d -ej %f 2 1 -ej %f 4 3 -ej %f 6 5 -eM %f 0.5 -ej %f 5 3 -eM %f 0.5 -ej %f 3 1 -T" % (N_allpops, Theta, N_popA, N_popB, N_popC, N_popD, N_popE, N_popF, T5, T4, T3, T3, T2, T2, T1), shell=True, stdout=subprocess.PIPE).stdout
	output = com.read().splitlines()
	Model_3sp_05.append(np.array(ms2nparray(output)).swapaxes(0,1).reshape(N_allpops,-1).T)
	
	## save parameter values and trees
	par_3sp_05.write("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n" % (Pi1, Pi2,T1,T2,T3,T4,T5))
	#Randomly save a number of tree equivalent to the number of traits
	trees_3sp_05.append(random.sample(get_newick(output),100))
	print("Completed %d %% of Model 3sp_05 simulations" % (float(i)/Priorsize*100))

#Save the simulated SNP data
Model_3sp_05=np.array(Model_3sp_05)
np.savez_compressed('trainingSims/Model_3sp_05.npz', Model_3sp_05=Model_3sp_05)
del(Model_3sp_05)

###Migration 1
## Two species, AB and CDEF
for i in range(Priorsize):

	## Theta value of 0.005
	Theta = 0.005

	## Sample Pi values for the deepest node (Pi1) between 0.5 and 1 (different species).
	Pi1 = random.uniform(0.5,1)
	## obtain divergence time (tau - τ) priors using the Pi values. Pi = 1−2/3*e^(−2τ/θ); τ = -(θ*ln(-(3*Pi-3)/2)/2)
	tau1 = -(Theta*math.log(-(3*Pi1-3)/2)/2)  
	## Transform divergence times to 4Ne generations units (required by ms)
	T1 = 4*tau1/Theta

	## Sample Pi values for the second deepest node (Pi2) between 1/3 (minimum value -> 0.333) and 0.4 (same species).
	Pi2 = random.uniform(1/3,0.4)
	## obtain divergence time (tau - τ) priors using the Pi values. Pi = 1−2/3*e^(−2τ/θ); τ = -(θ*ln(-(3*Pi-3)/2)/2)
	tau2 = -(Theta*math.log(-(3*Pi2-3)/2)/2)  
	## Transform divergence times to 4Ne generations units (required by ms)
	T2 = 4*tau2/Theta

	## Sample the more recent splitting times with a smaller value than the more ancient ones and bigger than 0 (present)
	T3 = random.uniform(0,T2)
	T4 = random.uniform(0,T3)
	T5 = random.uniform(0,T4)

	## ms command
	com=subprocess.Popen("./ms %d 1000 -s 1 -t %f -I 6 %d %d %d %d %d %d -ej %f 2 1 -ej %f 4 3 -ej %f 6 5 -ej %f 5 3 -eM %f 1 -ej %f 3 1 -T" % (N_allpops, Theta, N_popA, N_popB, N_popC, N_popD, N_popE, N_popF, T5, T4, T3, T2, T2, T1), shell=True, stdout=subprocess.PIPE).stdout
	output = com.read().splitlines()
	Model_2sp_1.append(np.array(ms2nparray(output)).swapaxes(0,1).reshape(N_allpops,-1).T)

	## save parameter values and trees
	par_2sp_1.write("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n" % (Pi1, Pi2,T1,T2,T3,T4,T5))
	#Randomly save a number of tree equivalent to the number of traits
	trees_2sp_1.append(random.sample(get_newick(output),100))
	print("Completed %d %% of Model 2sp_1 simulations" % (float(i)/Priorsize*100))

#Save the simulated SNP data
Model_2sp_1=np.array(Model_2sp_1)
np.savez_compressed('trainingSims/Model_2sp_1.npz', Model_2sp_1=Model_2sp_1)
del(Model_2sp_1)


## Three species
for i in range(Priorsize):

	## Theta value of 0.005
	Theta = 0.005

	## Sample Pi values for the deepest node (Pi1) between 0.5 and 1 (different species).
	Pi1 = random.uniform(0.5,1)
	## obtain divergence time (tau - τ) priors using the Pi values. Pi = 1−2/3*e^(−2τ/θ); τ = -(θ*ln(-(3*Pi-3)/2)/2)
	tau1 = -(Theta*math.log(-(3*Pi1-3)/2)/2)  
	## Transform divergence times to 4Ne generations units (required by ms)
	T1 = 4*tau1/Theta

	## Sample Pi values for the second deepest node (Pi2) that are higer than 0.5 (different species) but lower than the deepest node (Pi1).
	Pi2 = random.uniform(0.5,Pi1)
	## obtain divergence time (tau - τ) priors using the Pi values. Pi = 1−2/3*e^(−2τ/θ); τ = -(θ*ln(-(3*Pi-3)/2)/2)
	tau2 = -(Theta*math.log(-(3*Pi2-3)/2)/2)  
	## Transform divergence times to 4Ne generations units (required by ms)
	T2 = 4*tau2/Theta

	## Sample Pi values for the oldest intraspecific node (Pi3) between 1/3 (minimum value -> 0.333) and 0.4 (same species).
	Pi3=random.uniform(1/3,0.4)
	## obtain divergence time (tau - τ) priors using the Pi values. Pi = 1−2/3*e^(−2τ/θ); τ = -(θ*ln(-(3*Pi-3)/2)/2)
	tau3 = -(Theta*math.log(-(3*Pi3-3)/2)/2)  
	## Transform divergence times to 4Ne generations units (required by ms)
	T3 = 4*tau3/Theta

	## Sample the more recent splitting times with a smaller value than the more ancient ones and bigger than 0 (present)
	T4 = random.uniform(0,T4)
	T5 = random.uniform(0,T5)

	## ms command
	com=subprocess.Popen("./ms %d 1000 -s 1 -t %f -I 6 %d %d %d %d %d %d -ej %f 2 1 -ej %f 4 3 -ej %f 6 5 -eM %f 1 -ej %f 5 3 -eM %f 1 -ej %f 3 1 -T" % (N_allpops, Theta, N_popA, N_popB, N_popC, N_popD, N_popE, N_popF, T5, T4, T3, T3, T2, T2, T1), shell=True, stdout=subprocess.PIPE).stdout
	output = com.read().splitlines()
	Model_3sp_1.append(np.array(ms2nparray(output)).swapaxes(0,1).reshape(N_allpops,-1).T)
	
	## save parameter values and trees
	par_3sp_1.write("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n" % (Pi1, Pi2,T1,T2,T3,T4,T5))
	#Randomly save a number of tree equivalent to the number of traits
	trees_3sp_1.append(random.sample(get_newick(output),100))
	print("Completed %d %% of Model 3sp_1 simulations" % (float(i)/Priorsize*100))

#Save the simulated SNP data
Model_3sp_1=np.array(Model_3sp_1)
np.savez_compressed('trainingSims/Model_3sp_1.npz', Model_3sp_1=Model_3sp_1)
del(Model_3sp_1)

#Save trees from all scenarios
trees=np.concatenate((trees_2sp_01,trees_3sp_01,trees_2sp_05,trees_3sp_05,trees_2sp_1,trees_3sp_1),axis=0)
np.savez_compressed('trees.npz', trees=trees)
