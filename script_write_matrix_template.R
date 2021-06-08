library(ape)
library(phangorn)
library(geiger)
library(tidyverse)
library(castor)

calls = "alarm"
q = 1 + 64 + 32

args = commandArgs(TRUE)

if(length(args) > 0)
{	
	q = as.numeric(commandArgs(TRUE)[1])
	
	if(length(args) > 1)
	{	
		calls = commandArgs(TRUE)[2]
	}
}

##############################
### Definition of the data ###
##############################

sem = read.tree(paste(paste("semantic_", calls, sep = ""), ".tre", sep = ""))

# Construct the lexicons (as lines in the Lexicons data frame)
n_meanings = length(c(sem$tip.label,sem$node.label)) 
Lexicons = expand.grid(rep(list(c(0,1)),n_meanings))

####################################
### Quadrillage de la likelihood ###
####################################

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

Q_ind = matrix(data = 1, nrow = 2^n_meanings, ncol = 2^n_meanings)

u = unique(c(Q_ind))

for(i in 1:nrow(Q_ind))
{
	for(j in 1:ncol(Q_ind))
	{
		if(Q_ind[i,j] != 0)
		{
			Q_ind[i,j] = paste("mu", Q_ind[i,j], sep = "")
			if (i == j)
			{
			  Q_ind[i,j] = 0
			}
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

n = n_meanings
n_param_norm_tmp =  0.1*sum(Q_ind != "0") * (n/2) * (2^n) / (2^n -1)
sink("n_param_norm_tmp")
cat(n_param_norm_tmp)
sink()
