library(ape)
library(phangorn)
library(geiger)
library(tidyverse)
library(castor)

calls = "food"
q = 1 + 64 + 32
q = 128

args = commandArgs(TRUE)

if(length(args) > 0)
{	
	calls = commandArgs(TRUE)[1]
	
	if(length(args) > 1)
	{	
		q = as.numeric(commandArgs(TRUE)[2])
	}
}

##############################
### Definition of the data ###
##############################

sem = read.tree(paste(paste("semantic_", calls, sep = ""), ".tre", sep = ""))

# Construct the lexicons (as lines in the Lexicons data frame)
n_meanings = length(c(sem$tip.label,sem$node.label)) 
Lexicons = expand.grid(rep(list(c(0,1)),n_meanings))

Q_ind = matrix(data = 0, nrow = 2^n_meanings, ncol = 2^n_meanings)


for (n in seq(n_meanings)) # Loop over all the meanings
{	 				
	
	for (lex in seq(nrow(Lexicons))) # Loop over all the lexicons
	{		
		
		if (Lexicons[lex,n]) {
		
			##########################
			#### Death of the call ###
			##########################
			
			Q_ind[lex, lex - 2^(n-1)] = Q_ind[lex, lex - 2^(n-1)] + 13
			
			
		} else
		{
		
			#########################
			### Birth of the call ###
			#########################
		
			Q_ind[lex, lex + 2^(n-1)] = Q_ind[lex, lex + 2^(n-1)] + 1 + sum(Lexicons[lex,])
		
		}		
	}
}

u = unique(c(Q_ind))
u = u[-1]
u = sort(u)
u = u[-length(u)]

for(i in 1:nrow(Q_ind))
{
	for(j in 1:ncol(Q_ind))
	{
		if(Q_ind[i,j] != 0)
		{
			Q_ind[i,j] = paste("mu", Q_ind[i,j], sep = "")
		}
	}
}


write.table(Q_ind, "Q_tmp", sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)
sink("n_param_tmp")
cat(length(u)-1)
sink()	
sink("param_tmp")
cat(u)
sink()	
