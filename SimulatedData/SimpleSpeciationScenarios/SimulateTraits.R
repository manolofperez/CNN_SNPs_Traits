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
library(phytools)
library(reticulate)
library(foreach)
library(doParallel)

#setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(20) #not to overload your computer
registerDoParallel(cl)

nsims <- 100 # Number of traits to simulate per tree per treeset
dir.create("traits")
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
    #drop one tip for each diploid individual
    tree <-drop.tip(tree, as.character(seq(0,length(tree$tip.label),by=2)))
    # Trait simulations 
    sims <- fastBM(tree, a = 0, sig2 = 0.06, nsim = 1) 		
    #order the traits according to the SNPs order
    ordered_sims <- sims[order(as.numeric(names(sims)))]
    return(ordered_sims)
}	

write.table(format(as.matrix(traits_BM), width = 6),"./traits/traits_BM.txt",row.names=F,col.names=F,quote = F)

traits_OU={}
traits_OU<-foreach(i=1:nrow(trees),.combine=rbind,.packages = c("phytools","ape","geiger")) %:%
  foreach(j=1:nsims,.combine=cbind,.packages = c("phytools","ape","geiger")) %dopar% {
    tree <- read.tree(text=trees[i,j])
    #drop one tip for each diploid individual
    tree <-drop.tip(tree, as.character(seq(0,length(tree$tip.label),by=2)))
    # Trait simulations 
    sims <- fastBM(tree, a = 0, sig2 = 0.06, alpha = 0.2, theta = 0, nsim = 1) 		
    #order the traits according to the SNPs order
    ordered_sims <- sims[order(as.numeric(names(sims)))]
    return(ordered_sims)		
}

write.table(format(as.matrix(traits_OU), width = 6, flag = "0"),"./traits/traits_OU.txt",row.names=F,col.names=F,quote = F)

# Discrete trait - multivariate (3 states)

traits_disc={}
traits_disc<-foreach(i=1:nrow(trees),.combine=rbind,.packages = c("phytools","ape","geiger")) %:%
  foreach(j=1:nsims,.combine=cbind,.packages = c("phytools","ape","geiger")) %dopar% {
    #tree <- read.tree(text=trees[i,sample(seq_len(ncol(trees)), size=1)])
    tree <- read.tree(text=trees[i,j])
    #drop one tip for each diploid individual
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
}	

write.table(as.matrix(traits_disc),"./traits/traits_disc.txt",row.names=F,col.names=F,quote=F)
stopImplicitCluster()	

