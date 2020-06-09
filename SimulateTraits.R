#!/usr/bin/env Rscript

# Script By: Manolo Perez
# A script to simulate continuous traits on several newick trees saved in a file.

library(ape)
library(geiger)
library(caper)
library(phytools)
library(diversitree)

nsim <- 10 # Number of traits to simulate per tree per treeset
dir.create("traits")
dir.create("traits/BM")
dir.create("traits/discrete")
dir.create("traits/OU")

### Simulate traits ###

# BM

trees <- read.tree("trees.txt")
for(j in 1:length(trees)) {
  tree <- trees[[j]]
  
  # Trait simulations 
  sims <- fastBM(tree, a = 0, sig2 = 0.06, nsim = nsim) 		
  ordered_sims <- sims[order(as.numeric(row.names(sims))), ]
  
  write.table(as.matrix(ordered_sims), paste0("./traits/BM/", j, ".txt"),row.names=F,col.names=F)
}	


# Discrete trait - multivariate

trees <- read.tree("trees.txt")
for(j in 1:length(trees)) {
		tree <- trees[[j]]
				
		# Trait simulations 		
		statecols <- vector()
		for(k in 1:nsim) {
			qq <- list(rbind(c(-.5, .5), c(.5, -.5)))
			msims <- sim.char(tree, qq, model="discrete")
			statecols <- c(statecols, msims)	
		}
		
		sims <- matrix(unlist(statecols), ncol=nsim)		
		rownames(sims) <- tree$tip.label	
		ordered_sims <- sims[order(as.numeric(row.names(sims))), ]
								
		write.table(as.matrix(ordered_sims), paste0("./traits/discrete/", j, ".txt"),row.names=F,col.names=F)
}	


# OU process

trees <- read.tree("trees.txt")
for(j in 1:length(trees)) {
		tree <- trees[[j]]
		
		# Trait simulations 
		sims <- fastBM(tree, a = 0, sig2 = 0.06, alpha = 0.2, theta = 0, nsim = nsim) 		
		ordered_sims <- sims[order(as.numeric(row.names(sims))), ]
		
		write.table(as.matrix(ordered_sims), paste0("./traits/OU/", j, ".txt"),row.names=F,col.names=F)
	}
