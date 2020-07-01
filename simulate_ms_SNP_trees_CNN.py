#!/usr/bin/python

## in order to use this code you have to have ms installed on your computer
## ms can be freely downloaded from:
## http://home.uchicago.edu/rhudson1/source/mksamples.html

import random
import os
import math
import shlex, subprocess
import numpy as np

def ms2nparray(xfile):
	g = list(xfile)
	k = [idx for idx,i in enumerate(g) if len(i) > 0 and i.startswith(b'//')]
	f = []
	for i in k:
		L = g[i+5:i+N_allpops+5]
		q = []
		for i in L:
			i = [int(j-48) for j in list(i)]
			i = np.array(i, dtype=np.int8)
			q.append(i)
		q = np.array(q)
		q = q.astype("int8")
		f.append(np.array(q))   
	return f

def get_newick(xfile):
	g = list(xfile)
	k = [idx for idx,i in enumerate(g) if len(i) > 0 and i.startswith(b'//')]
	for i in k:
		n = g[i+1]
	return n


### variable declarations

#define the number of simulations
Priorsize = 1000

## sample size of pop1.
N_popA = 10
## nDNA sample size of ITA.
N_popB = 10
## nDNA sample size of PMN.
N_popC = 10

## sample size for pops 1 and 2 combined.
N_popAB = N_popA + N_popB

## sample size for all pops combined.
N_allpops = N_popA + N_popB + N_popC

simModel1 = []
simModel2 = []
simModel3 = []
simModel4 = []

## create a file to store parameters and one to store the models
parameters = open("parameters.txt","w")
models = open("models.txt","w")
trees = open("trees.txt","w")
os.mkdir("trainingSims")

#Define default values for priors absent in some models.
T2=0
T1=0
M=0

### One species
for i in range(Priorsize):
	### Define parameters
	## Theta values from 1 to 15
	Theta = random.uniform(1,10)

	## ms command
	com=subprocess.Popen("./ms %d 1000 -s 1 -t %f -T" % (N_allpops, Theta), shell=True, stdout=subprocess.PIPE).stdout
	output = com.read().splitlines()
	simModel1.append(np.array(ms2nparray(output)).swapaxes(0,1).reshape(N_allpops,-1).T)
	## save parameter values and models
	parameters.write("%f\t%f\t%f\t%f\n" % (Theta, T2, T1, M))
	models.write("1\n")
	trees.write("%s\n" % (get_newick(output).decode('utf-8')))
	print("Completed %d %% of Model 1 simulations" % (float(i)/Priorsize*100))


simModel1=np.array(simModel1)
np.savez_compressed('trainingSims/simModel1.npz', simModel1=simModel1)
del(simModel1)

### Two species, A and B collapsed
for i in range(Priorsize):

	## Theta values from 1 to 15
	Theta = random.uniform(1,10)
	## divergence time prior following an uniform distribution from 0.5 to 5.
	T2 = random.uniform(2,5)

	## ms command
	com=subprocess.Popen("./ms %d 1000 -s 1 -t %f -I 2 %d %d -ej %f 1 2 -T" % (N_allpops, Theta, N_popAB, N_popC ,T2), shell=True, stdout=subprocess.PIPE).stdout
	output = com.read().splitlines()
	simModel2.append(np.array(ms2nparray(output)).swapaxes(0,1).reshape(N_allpops,-1).T)

	## save parameter values and models
	parameters.write("%f\t%f\t%f\t%f\n" % (Theta, T2, T1, M))
	models.write("2\n")
	trees.write("%s\n" % (get_newick(output).decode('utf-8')))
	print("Completed %d %% of Model 2 simulations" % (float(i)/Priorsize*100))



simModel2=np.array(simModel2)
np.savez_compressed('trainingSims/simModel2.npz', simModel2=simModel2)
del(simModel2)


### Three species, isolation only
for i in range(Priorsize):

	## Theta values from 1 to 15
	Theta = random.uniform(1,10)
	## divergence time prior following an uniform distribution from 0.5 to 5.
	T2 = random.uniform(2,5)
	T1 = random.uniform(0.5,T2)

	## nDNA ms's command
	com=subprocess.Popen("./ms %d 1000 -s 1 -t %f -I 3 %d %d %d -ej %f 1 2 -ej %f 3 2 -T" % (N_allpops, Theta, N_popA, N_popB ,N_popC, T1, T2), shell=True, stdout=subprocess.PIPE).stdout
	output = com.read().splitlines()
	simModel3.append(np.array(ms2nparray(output)).swapaxes(0,1).reshape(N_allpops,-1).T)
	## save parameter values
	parameters.write("%f\t%f\t%f\t%f\n" % (Theta, T2, T1, M))
	models.write("3\n")
	trees.write("%s\n" % (get_newick(output).decode('utf-8')))
	print("Completed %d %% of Model 3 simulations" % (float(i)/Priorsize*100))

simModel3=np.array(simModel3)
np.savez_compressed('trainingSims/simModel3.npz', simModel3=simModel3)

### Three species, isolation with migration
for i in range(Priorsize):

	## Theta values from 1 to 15
	Theta = random.uniform(1,10)
	## divergence time prior following an uniform distribution from 0.5 to 5.
	T2 = random.uniform(2,5)
	T1 = random.uniform(0.5,T2)
	M = random.uniform(1,5)

	## nDNA ms's command
	com=subprocess.Popen("./ms %d 1000 -s 1 -t %f -I 3 %d %d %d %f -ej %f 1 2 -ej %f 3 2 -T" % (N_allpops, Theta, N_popA, N_popB ,N_popC, M, T1, T2), shell=True, stdout=subprocess.PIPE).stdout
	output = com.read().splitlines()
	simModel4.append(np.array(ms2nparray(output)).swapaxes(0,1).reshape(N_allpops,-1).T)
	## save parameter values
	parameters.write("%f\t%f\t%f\t%f\n" % (Theta, T2, T1, M))
	models.write("4\n")
	trees.write("%s\n" % (get_newick(output).decode('utf-8')))
	print("Completed %d %% of Model 4 simulations" % (float(i)/Priorsize*100))

simModel3=np.array(simModel3)
np.savez_compressed('trainingSims/simModel4.npz', simModel4=simModel4)

