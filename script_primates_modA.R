library(ape)
library(mgcv)

tree = read.tree("primates.tre")
sem = read.tree("semantic.tre")

args <- commandArgs(TRUE)
index_set = as.numeric(args[1])
nb.call.root = 1

n_rep = 40000

for (rep in 1:n_rep)
{
	# Esp 1/10E6
	# exp -lamba*10E6 = 1/2
	lambda = rgamma(1, 0.2)
	# delta => mu / 1/5 lambda (si 5 traits) 
	delta = lambda * rlnorm(1, 0)

	# Initialisation of the languages
	language = matrix(0,nb.call.root,length(sem$edge[,2]) + 1)

	leaves = c()

	tree$node.names = c(tree$tip.label, tree$node.label)
	# Root
	root.label = length(tree$tip.label) + 1

	# Start is the vector of the children of the Root
	start <- which(tree$edge[, 1] == root.label)

	# Start.stock is the vector of the nodes, growing in the loop
	start.stock = c()
	start.stock = c(start.stock, start)

	index.stock = c()

	# traits is the list of the evolving characters
	traits = list()

	l = 0

	i = -1

	# Loop walking along the tree of the species/languages
	while (length(start.stock) > 0)
	{	
		i = i + 1

		if (i > 0)
		{
			# Index is the next Node used
			index = start.stock[1]
			index.stock = c(index.stock, index)

			# Start is the vector of the children of the Node
			start = which(tree$edge[, 1] == tree$edge[index, 2])
			start.stock = c(start.stock, start)

			if(length(start)  == 0)
			{
				leaves = c(leaves, i)
			}
		}


		while(length(start > 0))
		{
			l = l + 1

			# Transmission of the trait from the parent node to the children
			if (i == 0)
			{
				traits[[l]] = language
				wait.time = 1
			}
			else
			{
				traits[[l]] = traits[[i]]
				wait.time = tree$edge.length[index.stock[i]]*1000
			}

			time = 0
			
			while (time < wait.time)
			{
				
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

				# The root may appear
		
				new.word = traits[[l]][1,]
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
			
				# Loop walking along the semantic tree
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
				
					new.word = mat.new.words[vect.index.word[p],]
					new.word[tree.language$edge[index.language, 2]] = 1
					mat.new.words = rbind(mat.new.words, new.word)
					vect.wait.time = c(vect.wait.time, rexp(1,lambda))
			
					start.language.stock = c(start.language.stock, start.language)				
					start.language.stock = start.language.stock[-1]
				
				}

				nb.words = dim(traits[[l]])[1]

				# A word may disappear ?
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
						traits[[l]] = traits[[l]][-(index.event - (dim(mat.new.words)[1] - nb.words)),, drop = FALSE]
					}
					else
					{
						traits[[l]] = rbind(traits[[l]], mat.new.words[index.event,])
						traits[[l]] = uniquecombs(traits[[l]])
					}
				}
				nb.words = dim(traits[[l]])[1]
		
				time = time + min(vect.wait.time)
		
			}
	

			start = start[-1]
		}
	
		if (i > 0)
		{
			start.stock	= start.stock[-1]
		}
	}

	traits_out = list()
	for (i in 1:length(leaves))
	{
		traits_out[[i]] = traits[[leaves[i]]]
	} 
	#Writing of the output
	dput(traits_out, paste(paste('primates_output_', (index_set - 1)*n_rep + rep , sep = ""),'.csv', sep = "") ) 

}


