library(ape)
library(phangorn)
library(geiger)
library(tidyverse)
library(castor)

calls = "alarm"
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

param = c(1,2,4,8,16,32,64)
n_param = length(param)

### Code for q : ###
# q = 1
# lambda_death 		=	+ 64 	= + 2^(n_param - 1)
# lambda_birth 		= 	+ 32 	= + 2^(n_param - 2)
# lambda_root_death = 	+ 16	= + 2^(n_param - 3)
# lambda_root_birth =	+ 8		etc.
# lambda_birth_dom  = 	+ 4
# lambda_split 		= 	+ 2
# lambda_specialize = 	+ 1
# etc.

lambda_death 		= 	1  *    (q-1) %/% (2^(n_param-1) )
lambda_birth 		= 	2  * ( ((q-1) %% (2^(n_param-1))) %/% (2^(n_param-2) ) )
lambda_root_death 	= 	4  * ( ((q-1) %% (2^(n_param-2))) %/% (2^(n_param-3) ) )
lambda_root_birth 	= 	8  * ( ((q-1) %% (2^(n_param-3))) %/% (2^(n_param-4) ) )
lambda_birth_dom 	= 	16 * ( ((q-1) %% (2^(n_param-4))) %/% (2^(n_param-5) ) )
lambda_split 		= 	32 * ( ((q-1) %% (2^(n_param-5))) %/% (2^(n_param-6) ) )
lambda_specialize 	= 	64 * ( ((q-1) %% (2^(n_param-6))) %/% (2^(n_param-7) ) )

{
	print(paste("lambda_death :", 		lambda_death))
	print(paste("lambda_birth :", 		lambda_birth))
	print(paste("lambda_root_death :", 	lambda_root_death))
	print(paste("lambda_root_birth :", 	lambda_root_birth))
	print(paste("lambda_birth_dom :", 	lambda_birth_dom))
	print(paste("lambda_split :", 		lambda_split))
	print(paste("lambda_specialize :", 	lambda_specialize))
}

Q_ind = matrix(data = 0, nrow = 2^n_meanings, ncol = 2^n_meanings)


for (n in seq(n_meanings)) # Loop over all the meanings
{	 				
	
	for (lex in seq(nrow(Lexicons))) # Loop over all the lexicons
	{		
		
		if (Lexicons[lex,n]) {
		
			##########################
			#### Death of the call ###
			##########################
			
			Q_ind[lex, lex - 2^(n-1)] = Q_ind[lex, lex - 2^(n-1)] + lambda_death	

			##########################
			#### Death of the root ###
			##########################
			
			AC = Ancestors(sem, n, type = "all")
			if (length(AC) == 0)
			{
				Q_ind[lex, lex - 2^(n-1)] = Q_ind[lex, lex - 2^(n-1)] + lambda_root_death
			}		
			
			#########################
			### Split of the call ###
			#########################
			
			CH = Children(sem, n)
			
			if (length(CH) > 0)
			{
				CH_0 = c()
				for (i in 1:length(CH))
				{
					if(Lexicons[lex,CH[i]] == 0)
					{
						CH_0 = c(CH_0, CH[i])
					}
				}
				
				if(length(CH_0) > 1)
				{
					for(i in 1:length(CH_0))
					{
						
						for (j in 1:length(CH_0))
						{
							index = 0
							index = index + 2^(CH_0[i]-1)
							
							if(i != j)
							{
								index = index + 2^(CH_0[j]-1)
							}
							
							if (Q_ind[lex, lex - 2^(n-1) + index] == 0)
							{
								Q_ind[lex, lex - 2^(n-1) + index] = Q_ind[lex, lex - 2^(n-1) + index] + lambda_split		
							}
						}	
					}
				}
			}
			
			
			##################################
			### Specialisation of the call ###
			##################################
		
			CH = Children(sem, n)
			
			if (length(CH) > 0)
			{	
				CH_0 = c()
				for (i in 1:length(CH))
				{
					if(Lexicons[lex,CH[i]] == 0)
					{
						Q_ind[lex, lex - 2^(n-1) + 2^(CH[i]-1)] = Q_ind[lex, lex - 2^(n-1) + 2^(CH[i]-1)] + lambda_specialize
					}
				}								
			}
		
		} else
		{
		
			#########################
			### Birth of the call ###
			#########################
		
			Q_ind[lex, lex + 2^(n-1)] = Q_ind[lex, lex + 2^(n-1)] + lambda_birth
		
			#########################
			### Birth of the root ###
			#########################
		
			AC = Ancestors(sem, n, type = "all")
			if (length(AC) == 0)
			{
				Q_ind[lex, lex + 2^(n-1)] = Q_ind[lex, lex + 2^(n-1)] + lambda_root_birth
			}		

			#################################
			### Birth of a dominated call ###
			#################################

			AC = Ancestors(sem, n, type = "all")
			DOMINATED = FALSE
			#DOMINATED = (length(AC)<1) # Set to FALSE, unless there is no ancestor

			for (ac in AC) 
			{ 
				if (Lexicons[lex, ac]) 
				{ 
					DOMINATED = TRUE 
				}
			}
		
			if (DOMINATED)
			{
				Q_ind[lex, lex + 2^(n-1)] = Q_ind[lex, lex + 2^(n-1)] + lambda_birth_dom
			}	
		}		
	}
}

Q_ind[Q_ind == 5] = 4
Q_ind[Q_ind == 10] = 8
Q_ind[Q_ind == 18] = 16

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
