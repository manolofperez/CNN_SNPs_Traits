#!/usr/bin/env Rscript

# Script By: Manolo Perez
# A script to simulate continuous traits on several newick trees saved in a file.

#devtools::install_github("KlausVigo/phangorn")
#install.packages("ape",repos="https://cloud.r-project.org")
#install.packages("diversitree",repos="https://cloud.r-project.org")
#install.packages("geiger",repos="https://cloud.r-project.org")
#install.packages("phytools",repos="https://cloud.r-project.org")
#install.packages("reticulate")
#install.packages(c( "foreach", "doParallel") )

library(ape)
library(geiger)
#library(caper)
library(phytools)
library(reticulate)
#library(diversitree)
library(foreach)
library(doParallel)

#setup parallel backend to use many processors
cores=detectCores()
#cl <- makeCluster(cores[1]-20) #not to overload your computer
cl <- makeCluster(20) #not to overload your computer
registerDoParallel(cl)

nsims <- 100 # Number of traits to simulate per tree per treeset
dir.create("traits")
#dir.create("traits/BM")
#dir.create("traits/discrete")
#dir.create("traits/OU")
np <- import("numpy",convert=F)
ntrees <- np$load("trees.npz")

trees<-matrix(ntrees$f[["trees"]], ncol=100,byrow=TRUE)

### Simulate traits ###
# BM
traits_BM={}
traits_BM<-foreach(i=1:nrow(trees),.combine=rbind,.packages = c("phytools","ape","geiger")) %:%
  foreach(j=1:nsims,.combine=cbind,.packages = c("phytools","ape","geiger")) %dopar% {
    #tree <- read.tree(text=trees[i,sample(seq_len(ncol(trees)), size=1)])
    tree <- read.tree(text=trees[i,j])
    tree <-drop.tip(tree, as.character(seq(0,length(tree$tip.label),by=2)))
    # Trait simulations 
    sims <- fastBM(tree, a = 0, sig2 = 0.06, nsim = 1) 		
    ordered_sims <- sims[order(as.numeric(names(sims)))]
    return(ordered_sims)
    #write.table(as.matrix(ordered_sims), paste0("./traits/BM/", formatC(j, width = 6, format = "d", flag = "0"), ".txt"),row.names=F,col.names=F)
}	

write.table(format(as.matrix(traits_BM), width = 6),"./traits/traits_BM.txt",row.names=F,col.names=F,quote = F)

traits_OU={}
traits_OU<-foreach(i=1:nrow(trees),.combine=rbind,.packages = c("phytools","ape","geiger")) %:%
  foreach(j=1:nsims,.combine=cbind,.packages = c("phytools","ape","geiger")) %dopar% {
    #tree <- read.tree(text=trees[i,sample(seq_len(ncol(trees)), size=1)])
    tree <- read.tree(text=trees[i,j])
    tree <-drop.tip(tree, as.character(seq(0,length(tree$tip.label),by=2)))
    # Trait simulations 
	sims <- fastBM(tree, a = 0, sig2 = 0.06, alpha = 0.2, theta = 0, nsim = 1) 		
	ordered_sims <- sims[order(as.numeric(names(sims)))]
	return(ordered_sims)		
	
	#write.table(as.matrix(ordered_sims), paste0("./traits/OU/", formatC(j, width = 6, format = "d", flag = "0"), ".txt"),row.names=F,col.names=F)
}

write.table(format(as.matrix(traits_OU), width = 6, flag = "0"),"./traits/traits_OU.txt",row.names=F,col.names=F,quote = F)

# Discrete trait - multivariate (3 states)

traits_disc={}
traits_disc<-foreach(i=1:nrow(trees),.combine=rbind,.packages = c("phytools","ape","geiger")) %:%
  foreach(j=1:nsims,.combine=cbind,.packages = c("phytools","ape","geiger")) %dopar% {
    #tree <- read.tree(text=trees[i,sample(seq_len(ncol(trees)), size=1)])
    tree <- read.tree(text=trees[i,j])
    tree <-drop.tip(tree, as.character(seq(0,length(tree$tip.label),by=2)))
    # Trait simulations 		
		statecols <- vector()
			qq <- matrix(c(-1, .5, .5, .5, -1, .5, .5, .5, -1),3,3)
			msims <- sim.history(tree, qq, message=FALSE)
			statecols <- c(statecols, msims$states)	
		
		sims <- matrix(unlist(statecols), ncol=1)		
		rownames(sims) <- tree$tip.label	
		ordered_sims <- sims[order(as.numeric(row.names(sims))), ]
								
		return(ordered_sims)
		
		#write.table(as.matrix(ordered_sims), paste0("./traits/discrete/", formatC(j, width = 6, format = "d", flag = "0"), ".txt"),row.names=F,col.names=F,quote=F)
}	

write.table(as.matrix(traits_disc),"./traits/traits_disc.txt",row.names=F,col.names=F,quote=F)
stopImplicitCluster()
## Discrete trait - multivariate (2 states)
#
#trees <- read.tree("trees.txt")
#for(j in 1:length(trees)) {
#		tree <- trees[[j]]
#				
#		# Trait simulations 		
#		statecols <- vector()
#		for(k in 1:nsim) {
#			qq <- list(rbind(c(-.5, .5), c(.5, -.5)))
#			msims <- sim.char(tree, qq, model="discrete")
#			statecols <- c(statecols, msims)	
#		}
#		
#		sims <- matrix(unlist(statecols), ncol=nsim)		
#		rownames(sims) <- tree$tip.label	
#		ordered_sims <- sims[order(as.numeric(row.names(sims))), ]
#								
#		write.table(as.matrix(ordered_sims), paste0("./traits/discrete/", formatC(j, width = 6, format = "d", flag = "0"), ".txt"),row.names=F,col.names=F)
#}	

