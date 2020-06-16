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
		L = g[i+5:i+nDNANsam+5]
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

## nDNA sample size of JFE.
nDNAJFE = 8
## nDNA sample size of ITA.
nDNAITA = 8
## nDNA sample size of PMN.
nDNAPMN = 8
## nDNA sample size of EDB.
nDNAEDB = 8
## nDNA sample size of BOV.
nDNABOV = 8
## nDNA sample size of COC.
nDNACOC = 8
## nDNA sample size of INA.
nDNAINA = 8
## nDNA sample size of ODA.
nDNAODA = 8
## nDNA sample size of MEN.
nDNAMEN = 8

## nDNA sample size of ITAPMN.
nDNAITAPMN = nDNAITA + nDNAPMN
## nDNA sample size of BOVEDB.
nDNABOVEDB = nDNAEDB + nDNABOV
## nDNA sample size of Central.
nDNACentral = nDNACOC + nDNAINA + nDNAODA + nDNAMEN

## nDNA sample sizes (number of alleles).
nDNANsam = nDNAITA + nDNAPMN + nDNAEDB + nDNABOV + nDNAJFE + nDNACOC + nDNAINA + nDNAODA + nDNAMEN
## number of years per generation
genlen = 15

simModel1 = []
simModel2 = []
simModel3 = []
## create a file to store parameters and one to store the models
parameters = open("parameters.txt","w")
models = open("models.txt","w")
trees = open("trees.txt","w")
os.mkdir("trainingSims")

#Define default values for priors absent in some models.
T5=0
T6=0

### Clade ES Geneland Model
for i in range(Priorsize):
	### Define parameters
	## Theta values from 1 to 15
	Theta = random.uniform(1,15)
	## divergence time prior set to 0 in this model.
	coalRootDivTime = 0
	## nDNA ms's command
	com=subprocess.Popen("./ms %d 500 -s 1 -t %f -I 4 %d %d %d %d -ej 0 2 4 -ej 0 3 4 -ej 0 1 4 -T" % (nDNANsam, Theta, nDNAJFE, nDNAITAPMN ,nDNABOVEDB , nDNACentral), shell=True, stdout=subprocess.PIPE).stdout
	output = com.read().splitlines()
	simModel1.append(np.array(ms2nparray(output)).swapaxes(0,1).reshape(nDNANsam,-1).T)
	## save parameter values and models
	parameters.write("%f\t%f\n" % (Theta, coalRootDivTime))
	models.write("1\n")
	trees.write("%s\n" % (get_newick(output).decode('utf-8')))
	print("Completed %d %% of Model 1 simulations" % (float(i)/Priorsize*100))


simModel1=np.array(simModel1)
np.savez_compressed('trainingSims/simModel1.npz', simModel1=simModel1)
del(simModel1)

### Clade ES Splitter Model
for i in range(Priorsize):

	## Theta values from 1 to 15
	Theta = random.uniform(1,15)
	## divergence time prior following an uniform distribution from 0.5 to 5.
	coalRootDivTime = random.uniform(0.5,5)
	coalT1=random.uniform(0,coalRootDivTime)
	coalT2=random.uniform(0,coalT1)
	## nDNA ms's command
	com=subprocess.Popen("./ms %d 500 -s 1 -t %f -I 4 %d %d %d %d -ej %f 2 4 -ej %f 3 4 -ej %f 1 4 -T" % (nDNANsam, Theta, nDNAJFE, nDNAITAPMN ,nDNABOVEDB , nDNACentral, coalT2, coalT1, coalRootDivTime), shell=True, stdout=subprocess.PIPE).stdout
	output = com.read().splitlines()
	simModel2.append(np.array(ms2nparray(output)).swapaxes(0,1).reshape(nDNANsam,-1).T)

	## save parameter values and models
	parameters.write("%f\t%f\n" % (Theta, coalRootDivTime))
	models.write("2\n")
	trees.write("%s\n" % (get_newick(output).decode('utf-8')))
	print("Completed %d %% of Model 2 simulations" % (float(i)/Priorsize*100))



simModel2=np.array(simModel2)
np.savez_compressed('trainingSims/simModel2.npz', simModel2=simModel2)
del(simModel2)


### Clade ES Lumper Model
for i in range(Priorsize):

	## Theta values from 1 to 15
	Theta = random.uniform(1,15)
	## divergence time prior following an uniform distribution from 0.5 to 5.
	coalRootDivTime = random.uniform(0.5,5)

	## nDNA ms's command
	com=subprocess.Popen("./ms %d 500 -s 1 -t %f -I 4 %d %d %d %d -ej 0 2 4 -ej 0 3 4 -ej %f 1 4 -T" % (nDNANsam, Theta, nDNAJFE, nDNAITAPMN ,nDNABOVEDB , nDNACentral, coalRootDivTime), shell=True, stdout=subprocess.PIPE).stdout
	output = com.read().splitlines()
	simModel3.append(np.array(ms2nparray(output)).swapaxes(0,1).reshape(nDNANsam,-1).T)
	## save parameter values
	parameters.write("%f\t%f\n" % (Theta, coalRootDivTime))
	models.write("3\n")
	trees.write("%s\n" % (get_newick(output).decode('utf-8')))
	print("Completed %d %% of Model 3 simulations" % (float(i)/Priorsize*100))

simModel3=np.array(simModel3)
np.savez_compressed('trainingSims/simModel3.npz', simModel3=simModel3)
