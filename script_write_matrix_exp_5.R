library(ape)
library(phangorn)
library(geiger)
library(tidyverse)
library(castor)

calls = "alarm"
q = 97

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
			
			Q_ind[lex, lex - 2^(n-1)] = Q_ind[lex, lex - 2^(n-1)] + 1
			
			#########################################
			#### Death depending on the hierarchy ###
			#########################################	
			
			AC = Ancestors(sem, n, type = "all")
			if (length(AC) == 0)
			{
				Q_ind[lex, lex - 2^(n-1)] = Q_ind[lex, lex - 2^(n-1)] + 2
			}		

			AC = Ancestors(sem, n, type = "all")
			if (length(AC) > 0)
			{
				Q_ind[lex, lex - 2^(n-1)] = Q_ind[lex, lex - 2^(n-1)] + 4
			}
			
			AC = Ancestors(sem, n, type = "all")
			if (length(AC) > 1)
			{
				Q_ind[lex, lex - 2^(n-1)] = Q_ind[lex, lex - 2^(n-1)] + 8
			}	
					
		} else
		{
		
			#########################
			### Birth of the call ###
			#########################
		
			Q_ind[lex, lex + 2^(n-1)] = Q_ind[lex, lex + 2^(n-1)] + 16
		
			########################################
			### Birth depending on the hierarchy ###
			########################################
		
			AC = Ancestors(sem, n, type = "all")
			if (length(AC) == 0)
			{
				Q_ind[lex, lex + 2^(n-1)] = Q_ind[lex, lex + 2^(n-1)] + 32
			}		

			AC = Ancestors(sem, n, type = "all")
			if (length(AC) > 0)
			{
				Q_ind[lex, lex + 2^(n-1)] = Q_ind[lex, lex + 2^(n-1)] + 64
			}
			
			AC = Ancestors(sem, n, type = "all")
			if (length(AC) > 1)
			{
				Q_ind[lex, lex + 2^(n-1)] = Q_ind[lex, lex + 2^(n-1)] + 128
			}	
		}		
	}
}


Q_ind[Q_ind == 1 + 2] = 1
Q_ind[Q_ind == 1 + 4] = 2
Q_ind[Q_ind == 1 + 4 + 8] = 3
Q_ind[Q_ind == 16 + 32] = 4
Q_ind[Q_ind == 16 + 64] = 5
Q_ind[Q_ind == 16 + 64 + 128] = 6

u = unique(c(Q_ind))
u = u[-1]
u = sort(u)

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
cat(length(u))
sink()	
sink("param_tmp")
cat(u)
sink()	
