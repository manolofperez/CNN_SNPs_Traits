#!/usr/bin/python3

## in order to use this code you have to have ms installed on your computer
## ms can be freely downloaded from:
## http://home.uchicago.edu/rhudson1/source/mksamples.html

#Import required libraries
import random
import os
import math
import shlex, subprocess
import numpy as np

##define a function to read ms' simulations and transform them into a NumPy array.    
def ms2nparray(xfile):
	g = list(xfile)
	k = [idx for idx,i in enumerate(g) if len(i) > 0 and i.startswith(b'//')]
	f = []
	for i in k:
		L = g[i+5:i+N_allpops+5]
		q = []
		for i in L:
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

## sample size of Euphorbia balsamifera subsp. adenensis.
N_adenensis = 19
## sample size of Euphorbia balsamifera subsp. balsamifera.
N_balsamifera = 80
## sample size of Euphorbia balsamifera subsp. sepium.
N_sepium = 10

## sample size for all pops combined.
N_allpops = N_adenensis + N_balsamifera + N_sepium

## create files to store parameters, trees and the models
os.mkdir("trainingSims")
par_1sp = open("par_1sp.txt","w")
par_2spMorph = open("par_2spMorph.txt","w")
par_2spPhylo = open("par_2spPhylo.txt","w")
par_3sp = open("par_3sp.txt","w")
## create lists to store trees from each scenario
trees_1sp = []
trees_2spMorph = []
trees_2spPhylo = []
trees_3sp = []
## create lists to store simulations from each scenario
Model_1sp = []
Model_2spMorph = []
Model_2spPhylo = []
Model_3sp = []

####Simulate the species delimitation scenarios####

### Scenario 1: All populations belong to one species
for i in range(Priorsize):
	### Define parameters
	#Theta per site per generation (based on the estimates of Rincón-Barrado 2022). We divided the values by the average length of the sequences.
	Theta = random.uniform(1,5)/400
	
	## Sample Pi values for the deepest node (Pi1) between 1/3 (minimum value -> 0.333) and 0.4 (same species).
	Pi1 = random.uniform(1/3,0.4)
	
	## obtain divergence time (tau - τ) priors using the Pi values. Pi = 1−2/3*e^(−2τ/θ); τ = -(θ*ln(-(3*Pi-3)/2)/2)
	tau1 = -(Theta*math.log(-(3*Pi1-3)/2)/2) 
	## Transform divergence times to 4Ne generations units (required by ms)
	T1 = 4*tau1/Theta

	## Sample Pi values for the second deepest node (Pi2) with a smaller value than deepest one (Pi1).
	Pi2 = random.uniform(1/3,Pi1)
	## obtain divergence time (tau - τ) priors using the Pi values. Pi = 1−2/3*e^(−2τ/θ); τ = -(θ*ln(-(3*Pi-3)/2)/2)
	tau2 = -(Theta*math.log(-(3*Pi2-3)/2)/2)  
	## Transform divergence times to 4Ne generations units (required by ms)
	T2 = 4*tau2/Theta

	## ms command
	com=subprocess.Popen("./ms %d 429 -s 1 -t %f -I 3 %d %d %d -ej %f 1 2 -ej %f 2 3 -T" % (N_allpops, Theta, N_adenensis, N_balsamifera, N_sepium, T2, T1), shell=True, stdout=subprocess.PIPE).stdout
	# read ms output
	output = com.read().splitlines()
	# save the SNPs output as a NumPy array
	Model_1sp.append(np.array(ms2nparray(output)).swapaxes(0,1).reshape(N_allpops,-1).T)

	## save parameter values
	par_1sp.write("%f\t%f\t%f\t%f\n" % (Pi1, Pi2,T1,T2))
	#Randomly save a number of trees equivalent to the number of traits
	trees_1sp.append(random.sample(get_newick(output),4))
	print("Completed %d %% of Model 1 simulations" % (float(i)/Priorsize*100))

#Compress and save the simulated SNP data
Model_1sp=np.array(Model_1sp)
np.savez_compressed('trainingSims/Model_1sp.npz', Model_1sp=Model_1sp)
del(Model_1sp)

### Scenario 2: Two species, according to morphology
for i in range(Priorsize):

	### Define parameters
	#Theta per site per generation (based on the estimates of Rincón-Barrado 2022). We divided the values by the average length of the sequences.
	Theta = random.uniform(1,5)/400
	
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

	## ms command
	com=subprocess.Popen("./ms %d 429 -s 1 -t %f -I 3 %d %d %d -ej %f 2 3 -ej %f 1 3 -T" % (N_allpops, Theta, N_adenensis, N_balsamifera, N_sepium, T2, T1), shell=True, stdout=subprocess.PIPE).stdout
	# read ms output
	output = com.read().splitlines()
	# save the SNPs output as a NumPy array
	Model_2spMorph.append(np.array(ms2nparray(output)).swapaxes(0,1).reshape(N_allpops,-1).T)

	## save parameter values
	par_2spMorph.write("%f\t%f\t%f\t%f\n" % (Pi1, Pi2,T1,T2))
	#Randomly save a number of tree equivalent to the number of traits
	trees_2spMorph.append(random.sample(get_newick(output),4))
	print("Completed %d %% of Model 2 simulations" % (float(i)/Priorsize*100))

#Compress and save the simulated SNP data
Model_2spMorph=np.array(Model_2spMorph)
np.savez_compressed('trainingSims/Model_2spMorph.npz', Model_2spMorph=Model_2spMorph)
del(Model_2spMorph)

### Scenario 3: Two species, Phylogenomics
for i in range(Priorsize):

	### Define parameters
	#Theta per site per generation (based on the estimates of Rincón-Barrado 2022). We divided the values by the average length of the sequences.
	Theta = random.uniform(1,5)/400
	
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

	## ms command
	com=subprocess.Popen("./ms %d 429 -s 1 -t %f -I 3 %d %d %d -ej %f 1 2 -ej %f 2 3 -T" % (N_allpops, Theta, N_adenensis, N_balsamifera, N_sepium, T2, T1), shell=True, stdout=subprocess.PIPE).stdout
	# read ms output
	output = com.read().splitlines()
	# save the SNPs output as a NumPy array
	Model_2spPhylo.append(np.array(ms2nparray(output)).swapaxes(0,1).reshape(N_allpops,-1).T)
	
	## save parameter values and trees
	par_2spPhylo.write("%f\t%f\t%f\t%f\n" % (Pi1, Pi2,T1,T2))
	#Randomly save a number of tree equivalent to the number of traits
	trees_2spPhylo.append(random.sample(get_newick(output),4))
	print("Completed %d %% of Model 3 simulations" % (float(i)/Priorsize*100))

#Compress and save the simulated SNP data
Model_2spPhylo=np.array(Model_2spPhylo)
np.savez_compressed('trainingSims/Model_2spPhylo.npz', Model_2spPhylo=Model_2spPhylo)
del(Model_2spPhylo)


### Scenario 4: Three species
for i in range(Priorsize):

	### Define parameters
	#Theta per site per generation (based on the estimates of Rincón-Barrado 2022). We divided the values by the average length of the sequences.
	Theta = random.uniform(1,5)/400

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

	## ms command
	com=subprocess.Popen("./ms %d 429 -s 1 -t %f -I 3 %d %d %d -ej %f 1 2 -ej %f 2 3 -T" % (N_allpops, Theta, N_adenensis, N_balsamifera, N_sepium, T2, T1), shell=True, stdout=subprocess.PIPE).stdout
	# read ms output
	output = com.read().splitlines()
	# save the SNPs output as a NumPy array
	Model_3sp.append(np.array(ms2nparray(output)).swapaxes(0,1).reshape(N_allpops,-1).T)
	
	## save parameter values and trees
	par_3sp.write("%f\t%f\t%f\t%f\n" % (Pi1, Pi2,T1,T2))
	#Randomly save a number of tree equivalent to the number of traits
	trees_3sp.append(random.sample(get_newick(output),4))
	print("Completed %d %% of Model 4 simulations" % (float(i)/Priorsize*100))

#Compress and save the simulated SNP data
Model_3sp=np.array(Model_3sp)
np.savez_compressed('trainingSims/Model_3sp.npz', Model_3sp=Model_3sp)
# Compress and save trees from all scenarios
trees=np.concatenate((trees_1sp,trees_2spMorph,trees_2spPhylo,trees_3sp),axis=0)
np.savez_compressed('trees.npz', trees=trees)
