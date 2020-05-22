setwd("/Users/manolo/Desktop/SimulateTraits/test")
getwd()

# Script By: Manolo Perez - Modified from Michael Harvey
# A script to simulate continuous traits on several newick trees saved in a file.

library(ape)
library(caper)
library(phytools)
library(diversitree)

nsim <- 10 # Number of traits to simulate per tree per treeset
dir.create("traits")
dir.create("traits/BM")
dir.create("traits/discrete")
dir.create("traits/OUstrong")
dir.create("traits/OUweak")

### Simulate traits ###

# BM (uncorrelated)

trees <- read.tree("trees.txt")
for(j in 1:length(trees)) {
  tree <- trees[[j]]
  
  # Trait simulations 
  vv <- vcv.phylo(tree)
  sims <- t(rmvnorm(nsim, sigma=0.06*vv))
  write.table(as.matrix(sims), paste0("./traits/BM/", j, ".txt"),row.names=F,col.names=F)
}	


# Discrete trait (2 normal distributions) #not discrete

trees <- read.tree("trees.txt")
for(j in 1:length(trees)) {
		tree <- trees[[j]]
				
		# Trait simulations 		
		statecols <- vector()
		for(k in 1:nsim) {
			states <- sim.character(tree, c(0.1,0.1), x0=0, model="mk2")
			states[states == 0] <- rnorm(length(states[states == 0]), mean = 0, sd = 2) # Normal dist. for state 0
			states[states == 1] <- rnorm(length(states[states == 1]), mean = 5, sd = 4) # Normal dist. for state 1
			statecols <- c(statecols, states[tree$tip.label])	
		}
		
		sims <- matrix(unlist(statecols), ncol=nsim)		
		rownames(sims) <- tree$tip.label	
		new_sims <- sims[ order(as.numeric(row.names(sims))), ]
								
		write.table(as.matrix(new_sims), paste0("./traits/discrete/", j, ".txt"),row.names=F,col.names=F)
}	


# OU process (strong)

trees <- read.tree("trees.txt")
for(j in 1:length(trees)) {
		tree <- trees[[j]]
		
		# Trait simulations 
		sims <- fastBM(tree, a = 0, sig2 = 0.06, alpha = 0.2, theta = 0, nsim = nsim) 		
		new_sims <- sims[ order(as.numeric(row.names(sims))), ]
		
		write.table(as.matrix(new_sims), paste0("./traits/OUstrong/", j, ".txt"),row.names=F,col.names=F)
	}

# OU process (weak)

trees <- read.tree("trees.txt")
for(j in 1:length(trees)) {
		tree <- trees[[j]]
		
		# Trait simulations 
		sims <- fastBM(tree, a = 0, sig2 = 0.06, alpha = 0.002, theta = 0, nsim = nsim) 		
		new_sims <- sims[ order(as.numeric(row.names(sims))), ]
		
		write.table(as.matrix(new_sims), paste0("./traits/OUweak/", j, ".txt"),row.names=F,col.names=F)
	}


##################################

## BM (multirate)
#
#treefiles <- list.files("./trees/", pattern="*.tre", full.names=T)
#for(i in 1:length(treefiles)){
#  trees <- read.tree(treefiles[i])
#  treeset <- strsplit(strsplit(strsplit(treefiles[i], "//")[[1]][2], ".tre")[[1]][1], "_")[[1]][2]
#  if(length(trees) < 50) {
#    trees <- list(trees)
#  }
#  for(j in 1:length(trees)) {
#    tree <- trees[[j]]
#    
#    # Isolate clade of > 10% of tips
#    subtrees <- subtrees(tree)
#    large.subtrees <- vector()
#    for (k in 1:length(subtrees)) {
#      if((subtrees[[k]]$Nnode+1) > ((0.1*tree$Nnode)+1)) {
#        large.subtrees <- c(large.subtrees, subtrees[k])				
#      }
#    }
#    subclade <- large.subtrees[sample(1:length(large.subtrees), 1)][[1]]
#    
#    # Increase branch lengths for that clade
#    branches <- which.edge(tree, subclade$tip.label)
#    tree$edge.length[branches] <- tree$edge.length[branches]*5
#    
#    # Trait simulations 
#    vv <- vcv.phylo(tree)
#    sims <- t(rmvnorm(nsim, sigma=0.06*vv))
#    rownames(sims) <- rownames(vv)	
#    
#    write.table(as.matrix(sims), paste0("./traits/traits_", treeset, "_BMmultirate_", j, ".txt"))
#  }	
#}
#
## BM - jump in mean
#
#treefiles <- list.files("./trees/", pattern="*.tre", full.names=T)
#for(i in 1:length(treefiles)){
#  trees <- read.tree(treefiles[i])
#  treeset <- strsplit(strsplit(strsplit(treefiles[i], "//")[[1]][2], ".tre")[[1]][1], "_")[[1]][2]
#  if(length(trees) < 50) {
#    trees <- list(trees)
#  }
#  for(j in 1:length(trees)) {
#    tree <- trees[[j]]
#    
#    # Trait simulations 
#    vv <- vcv.phylo(tree)
#    sims <- t(rmvnorm(nsim, sigma=0.06*vv))
#    rownames(sims) <- rownames(vv)	
#    
#    # Isolate clade of > 10% of tips
#    subtrees <- subtrees(tree)
#    large.subtrees <- vector()
#    for (k in 1:length(subtrees)) {
#      if((subtrees[[k]]$Nnode+1) > ((0.1*tree$Nnode)+1)) {
#        large.subtrees <- c(large.subtrees, subtrees[k])				
#      }
#    }
#    subclade <- large.subtrees[sample(1:length(large.subtrees), 1)][[1]]
#    
#    # Jump in mean for that clade
#    sims[rownames(sims) %in% subclade$tip.label,] <- sims[rownames(sims) %in% subclade$tip.label,]+0.3
#    
#    write.table(as.matrix(sims), paste0("./traits/traits_", treeset, "_BMjump_", j, ".txt"))
#  }	
#}
#
## No phylogenetic signal
#
#treefiles <- list.files("./trees/", pattern="*.tre", full.names=T)
#for(i in 1:length(treefiles)){
#  trees <- read.tree(treefiles[i])
#  treeset <- strsplit(strsplit(strsplit(treefiles[i], "//")[[1]][2], ".tre")[[1]][1], "_")[[1]][2]
#  if(length(trees) < 50) {
#    trees <- list(trees)
#  }
#  for(j in 1:length(trees)) {
#    tree <- trees[[j]]
#    
#    # Trait simulations 
#    sim.vals <- rnorm(length(tree$tip.label)*50, mean = 0, sd = max(branching.times(tree)))
#    sims <- matrix(unlist(sim.vals), ncol=50)
#    rownames(sims) <- tree$tip.label	
#    
#    write.table(as.matrix(sims), paste0("./traits/traits_", treeset, "_nosignal_", j, ".txt"))
#  }	
#}
#
## 1 clade with no phylogenetic signal
#
#treefiles <- list.files("./trees/", pattern="*.tre", full.names=T)
#for(i in 1:length(treefiles)){
#  trees <- read.tree(treefiles[i])
#  treeset <- strsplit(strsplit(strsplit(treefiles[i], "//")[[1]][2], ".tre")[[1]][1], "_")[[1]][2]
#  if(length(trees) < 50) {
#    trees <- list(trees)
#  }
#  for(j in 1:length(trees)) {
#    tree <- trees[[j]]
#    
#    # Trait simulations 
#    vv <- vcv.phylo(tree)
#    sims <- t(rmvnorm(nsim, sigma=0.06*vv))
#    rownames(sims) <- rownames(vv)		
#    
#    # Isolate clade of > 10% of tips
#    subtrees <- subtrees(tree)
#    large.subtrees <- vector()
#    for (k in 1:length(subtrees)) {
#      if((subtrees[[k]]$Nnode+1) > ((0.1*tree$Nnode)+1)) {
#        large.subtrees <- c(large.subtrees, subtrees[k])				
#      }
#    }
#    subclade <- large.subtrees[sample(1:length(large.subtrees), 1)][[1]]
#    
#    # No signal in that one clade
#    norm.vals <- rnorm(length(subclade$tip.label)*50, mean = 0, sd = max(branching.times(subclade)))
#    sims[rownames(sims) %in% subclade$tip.label,] <- norm.vals
#    
#    write.table(as.matrix(sims), paste0("./traits/traits_", treeset, "_1cladenosignal_", j, ".txt"))
#  }	
#}
#
## 1 clade fixed
#
#treefiles <- list.files("./trees/", pattern="*.tre", full.names=T)
#for(i in 1:length(treefiles)){
#  trees <- read.tree(treefiles[i])
#  treeset <- strsplit(strsplit(strsplit(treefiles[i], "//")[[1]][2], ".tre")[[1]][1], "_")[[1]][2]
#  if(length(trees) < 50) {
#    trees <- list(trees)
#  }
#  for(j in 1:length(trees)) {
#    tree <- trees[[j]]
#    
#    # Trait simulations 
#    vv <- vcv.phylo(tree)
#    sims <- t(rmvnorm(nsim, sigma=0.06*vv))
#    rownames(sims) <- rownames(vv)		
#    
#    # Isolate clade of > 10% of tips
#    subtrees <- subtrees(tree)
#    large.subtrees <- vector()
#    for (k in 1:length(subtrees)) {
#      if((subtrees[[k]]$Nnode+1) > ((0.1*tree$Nnode)+1)) {
#        large.subtrees <- c(large.subtrees, subtrees[k])				
#      }
#    }
#    subclade <- large.subtrees[sample(1:length(large.subtrees), 1)][[1]]
#    
#    # Fix that clade to the value of one of its tips
#    val <- sims[sample(1:length(subclade$tip.label), 1),]
#    rep.val <- rep(val, each=nrow(sims[rownames(sims) %in% subclade$tip.label,]))
#    sims[rownames(sims) %in% subclade$tip.label,] <- rep.val
#    
#    write.table(as.matrix(sims), paste0("./traits/traits_", treeset, "_1cladefixed_", j, ".txt"))
#  }	
#}
#
#
#
#