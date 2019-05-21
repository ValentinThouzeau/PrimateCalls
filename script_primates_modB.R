library(ape)
library(mgcv)
library(MCMCpack)
library(MASS)


#######################
### DATA PARAMETERS ###
#######################

# Phylogenetic tree of the species
tree = read.tree("primates.tre")

# Semantic tree of the repertoire
sem = read.tree("semantic.tre")

#####################################
### REPETITION OF THE SIMULATIONS ###
#####################################

n_rep = 40000

for (rep in 1:n_rep)
{
	
	########################
	### PRIOR PARAMETERS ###
	########################
	
	# Birth rate
	lambda = rgamma(1, 0.2)
	
	# Death rate
	delta = lambda * rlnorm(1, 0)


	###################################################
	### INITIALISATION OF THE SIMULATION PARAMETERS ###
	###################################################
	
	# Initialisation
	nb.call.root = 1
	traits = c()
	traits[[1]] = matrix(0,nb.call.root,length(sem$edge[,2]) + 1)
	
	tips = c()
	tips2 = c()
	nodes = c()
	

	# Vector of the nodes growing during the loop
	node.stock = 0

	# Index of the current node
	i = 1
	
	# Index of the children nodes
	l = 1
	
	########################################
	### LOOP OF EVOLUTION ALONG THE TREE ###
	########################################
	
	while (length(node.stock) > 0)
	{		
		################################
		### EVOLUTION ALONG A BRANCH ###
		################################
 		
		if (i > 1)
		{ 
			# Length of the branch 
			wait.time = tree$edge.length[node.stock[1]]*1000
		}
		else 
		{ 
			# Time for an ancestral equilibrium
			wait.time = 1000
		}
		
		time = 0
		
		while (time < wait.time)
		{							
			
			#########################################
			### CONSTRUCTION OF THE SEMANTIC TREE ###	
			#########################################
							
			mat.new.words = c()
			vect.wait.time = c()

			tree.language = sem
			tree.language$node.names = c(tree.language$tip.label, tree.language$node.label)
	
			# Root of the semantic tree
			root.label.language = length(tree.language$tip.label) + 1

			# Start is the vector of the children of the Root
			start.language <- which(tree.language$edge[, 1] == root.label.language)

			# Start.stock is the vector of the nodes, growing in the loop
			start.language.stock = c()
			start.language.stock = c(start.language.stock, start.language)

			index.language.stock = c()


			################################
			### APPARITION OF THE ROOT ? ###
			################################
			
			new.word = traits[[i]][1,]
			new.word[length(tree.language$tip.label) + 1] = 1
			mat.new.words = rbind(mat.new.words, new.word)
			vect.wait.time = c(vect.wait.time, rexp(1,lambda))
			
			vect.index.word = c()
			vect.index.word = c(vect.index.word, 0)
			c = length(start.language)
			p = 1
			while (c > 0)
			{
				c = c - 1
				vect.index.word = c(vect.index.word, p)
			}	
		
			##############################
			### APPARITION OF A WORD ? ###
			##############################
			
			while (length(start.language.stock) > 0)
			{	
				p = p + 1

				# Index is the next Node used
				index.language = start.language.stock[1]
				index.language.stock = c(index.language.stock, index.language)

				# Start is the vector of the children of the Node
				start.language = which(tree.language$edge[, 1] == tree.language$edge[index.language, 2])

				c = length(start.language)
				while (c > 0)
				{
					c = c - 1
					vect.index.word = c(vect.index.word, p)
				}
			
				if (length(which(apply(traits[[i]], 1, function(x) identical(x, mat.new.words[vect.index.word[p],])))) > 0)
				{
					new.word = mat.new.words[vect.index.word[p],]
					new.word[tree.language$edge[index.language, 2]] = 1
					mat.new.words = rbind(mat.new.words, new.word)
					vect.wait.time = c(vect.wait.time, rexp(1,lambda))
				
					start.language.stock = c(start.language.stock, start.language)	
				}
							
				start.language.stock = start.language.stock[-1]
				
			}

			nb.words = dim(traits[[i]])[1]


			################################
			### DISAPPARTION OF A WORD ? ###
			################################
			
			if (nb.words > 1)
			{
				for (j in 1:(nb.words-1))
				{
					mat.new.words = rbind(mat.new.words, matrix(-1,1,length(sem$node.label) + length(sem$tip.label)))
					vect.wait.time = c(vect.wait.time, rexp( 1 , delta) )
				}
			}

			# Event
			index.event = which(vect.wait.time == min(vect.wait.time))
		
			if ((time + min(vect.wait.time)) < wait.time)
			{
				if (mat.new.words[index.event,1] == -1)
				{
					traits[[i]] = traits[[i]][-(index.event - (dim(mat.new.words)[1] - nb.words)),, drop = FALSE]
				}
				else
				{
					traits[[i]] = rbind(traits[[i]], mat.new.words[index.event,])
					traits[[i]] = uniquecombs(traits[[i]])
				}
			}
			nb.words = dim(traits[[i]])[1]
	
			time = time + min(vect.wait.time)
			
		}
		
		
		####################################
		### RECUPERATION OF THE CHILDREN ###
		####################################
		
		if (i > 1)
		{			
			# Node is the vector of the children of the node
			node = which(tree$edge[, 1] == tree$edge[node.stock[1], 2])
						
			if (length(node) < 1 )
			{
				tips = c(tips, tree$edge[node.stock[1],2])
				tips2 = c(tips2, i)
				nodes = c(nodes, tree$tip.label[tree$edge[node.stock[1],2]])
			}
			
			node.stock = c(node.stock, node)
		}
		else
		{
			# Root
			root.label = length(tree$tip.label) + 1

			# Vector of the children of the root
			node <- which(tree$edge[, 1] == root.label)
			node.stock = c(node.stock, node)
		}

		
		#############################
		### BRANCHING OF THE TREE ###
		#############################
	
		while(length(node) > 0)
		{
			l = l + 1
			
			# Transmission from the parent to the child node
			traits[[l]] = traits[[i]]
			
			node = node[-1]			
		}
		
	
		node.stock = node.stock[-1]
		
		i = i + 1
		
	}
	
	#############################
	### WRITING OF THE OUTPUT ###
	#############################
	
	traits_out = list()
	for (i in 1:length(tips2))
	{
		traits_out[[i]] = traits[[ tips2[ tips == i ] ]]
	} 
	#Writing of the output
	dput(traits_out, paste(paste('primates_output_', rep , sep = ""),'.csv', sep = "") ) 

}

