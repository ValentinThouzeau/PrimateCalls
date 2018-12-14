library(ape)
library(mgcv)

wd = getwd()
folderA = paste(wd, "/primatesA", sep = "")
folderB = paste(wd, "/primatesB", sep = "")
t = read.table("Data_alarm.csv", header = TRUE, sep = ",")

setwd(folderA)
traits_zero = dget("primates_output_0_alarm.csv")
sem = read.tree("semantic.tre")

nb_simulations = 40000


# WRITING OF THE CANONICAL MATRIX

tree.language = sem
tree.language$node.names = c(tree.language$tip.label, tree.language$node.label)

dim.tree.language = length(tree.language$node.names)

# Root of the semantic tree
root.label.language = length(tree.language$tip.label) + 1

# Start is the vector of the children of the Root
start.language <- which(tree.language$edge[, 1] == root.label.language)

# Start.stock is the vector of the nodes, growing in the loop
start.language.stock = c()
start.language.stock = c(start.language.stock, start.language)

index.language.stock = c()

vect.index.word = c()
vect.index.word = c(vect.index.word, 0)
c = length(start.language)
p = 1
while (c > 0)
{
	c = c - 1
	vect.index.word = c(vect.index.word, p)
}

mat.canonique = list() 
new.word = matrix(0,1,dim.tree.language)
new.word[tree.language$edge[which(tree.language$edge[,1] == root.label.language)][1]] = 1
mat.new.words = new.word
mat.canonique[[p]] = rbind(matrix(0,1,dim.tree.language), new.word) 

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
	mat.new.words = rbind(mat.new.words,t(as.matrix( new.word )))

	mat.canonique[[p]] = rbind(mat.canonique[[ vect.index.word[p] ]] , new.word)
	rownames(mat.canonique[[p]]) = NULL

	tree.language$edge[index.language, 2]

	# Next nodes are stocked in start.language.stock
	start.language.stock = c(start.language.stock, start.language)				
	start.language.stock = start.language.stock[-1]

}


mat.comp = list()

for (i in 1:length(mat.canonique))
{
	mat.comp[[i]] = combn(mat.canonique, i, simplify = FALSE)
	
}

mat.comp.bis = list()

c = 1
for (i in 1:length(mat.comp))
{
	for(j in 1:length(mat.comp[[i]]))
	{
		c = c + 1
		mat.temp = matrix(0,0,dim.tree.language)

		for (k in 1:length(mat.comp[[i]][[j]]))
		{
			mat.temp = rbind(mat.temp, mat.comp[[i]][[j]][[k]])	
		}
		mat.comp.bis[[c]] = unique(mat.temp)
	}
}

mat.comp.bis[[1]] = matrix(0,1,dim.tree.language)
mat.comp.bis = unique(mat.comp.bis)





# COMPUTING OF THE SUMSTAT

setwd(folderA)
sumstats = c()

for (numero_simulation in 1:nb_simulations)
{
	traits = dget(paste(paste("primates_output_",numero_simulation,sep =""),".csv",sep=""))

	if(identical(traits, traits_zero))
	{

	}
	else
	{
		vec.dist = c()

		for (i in 1:length(traits))
		{
			vec.dist = c(vec.dist, 0)
			c = 0
			for (j in 1:length(mat.comp.bis))
			{
				bol = c()
				if(dim(traits[[i]])[1] == dim(mat.comp.bis[[j]])[1])
				{					
					for(k in 1:dim(traits[[i]])[1])
					{

						bol = c(bol, any(apply(mat.comp.bis[[j]], 1, function(x, want) isTRUE(all.equal(x, want)), traits[[i]][k,])) )
					}
				}
				else
				{
					bol = c(bol, FALSE)
				}

				if(all(bol))
				{
					c = c + 1
				}
			}
			vec.dist[i] = c	
		}

		dist.sum = sum(vec.dist)
		dist.var = var(vec.dist)

		vec.nb.words.species = c()

		for (i in 1:length(traits))
		{
			vec.nb.words.one.species = dim(traits[[i]])[1]
		
			vec.nb.words.species = c(vec.nb.words.species, vec.nb.words.one.species)	
		}
	
		nb.mean = mean(vec.nb.words.species)
		nb.var = var(vec.nb.words.species)
		nb.min = min(vec.nb.words.species)
		nb.max = max(vec.nb.words.species)
	
		sumstats = rbind(sumstats, c(vec.dist, dist.sum, dist.var, vec.nb.words.species, nb.mean, nb.var, nb.min, nb.max))
	}

}

sA = sumstats




setwd(folderB)
sumstats = c()

for (numero_simulation in 1:nb_simulations)
{
	traits = dget(paste(paste("primates_output_",numero_simulation,sep =""),".csv",sep=""))

	if(identical(traits, traits_zero))
	{

	}
	else
	{
		vec.dist = c()

		for (i in 1:length(traits))
		{
			vec.dist = c(vec.dist, 0)
			c = 0
			for (j in 1:length(mat.comp.bis))
			{
				bol = c()
				if(dim(traits[[i]])[1] == dim(mat.comp.bis[[j]])[1])
				{					
					for(k in 1:dim(traits[[i]])[1])
					{

						bol = c(bol, any(apply(mat.comp.bis[[j]], 1, function(x, want) isTRUE(all.equal(x, want)), traits[[i]][k,])) )
					}
				}
				else
				{
					bol = c(bol, FALSE)
				}

				if(all(bol))
				{
					c = c + 1
				}
			}
			vec.dist[i] = c	
		}

		dist.sum = sum(vec.dist)
		dist.var = var(vec.dist)

		vec.nb.words.species = c()

		for (i in 1:length(traits))
		{
			vec.nb.words.one.species = dim(traits[[i]])[1]
		
			vec.nb.words.species = c(vec.nb.words.species, vec.nb.words.one.species)	
		}
	
		nb.mean = mean(vec.nb.words.species)
		nb.var = var(vec.nb.words.species)
		nb.min = min(vec.nb.words.species)
		nb.max = max(vec.nb.words.species)
	
		sumstats = rbind(sumstats, c(vec.dist, dist.sum, dist.var, vec.nb.words.species, nb.mean, nb.var, nb.min, nb.max))
	}

}

sB = sumstats






setwd(wd)

t[is.na(t)] = 0

real = list()
real[[1]] = matrix(0,1,dim(t)[2]-1)

t = unique(t)


c = 1

for (i in 1:(dim(t)[1]))
{
	row1 = t[i,1]
	
	if(i != dim(t)[1])
	{
		row2 = t[i+1,1]
	}
	else
	{
		row2 = "void"
	}


	d = c()
	for(j in 2:dim(t)[2])
	{
		d = c(d, as.numeric(t[i, j]))
	}

	real[[c]] = unique(rbind(real[[c]], d))


	if (row1 != row2)
	{
		rownames(real[[c]]) = NULL
	
		c = c + 1
	
		if(row2 != "void")
		{
			real[[c]] = matrix(0,1,dim(t)[2]-1)
		}
	}	
}

traits = real

vec.dist = c()

for (i in 1:length(traits))
{
	vec.dist = c(vec.dist, 0)
	c = 0
	for (j in 1:length(mat.comp.bis))
	{
		bol = c()
		if(dim(traits[[i]])[1] == dim(mat.comp.bis[[j]])[1])
		{					
			for(k in 1:dim(traits[[i]])[1])
			{

				bol = c(bol, any(apply(mat.comp.bis[[j]], 1, function(x, want) isTRUE(all.equal(x, want)), traits[[i]][k,])) )
			}
		}
		else
		{
			bol = c(bol, FALSE)
		}

		if(all(bol))
		{
			c = c + 1
		}
	}
	vec.dist[i] = c	
}

dist.sum = sum(vec.dist)
dist.var = var(vec.dist)

vec.nb.words.species = c()

for (i in 1:length(traits))
{
	vec.nb.words.one.species = dim(traits[[i]])[1]
	
	vec.nb.words.species = c(vec.nb.words.species, vec.nb.words.one.species)	
}

nb.mean = mean(vec.nb.words.species)
nb.var = var(vec.nb.words.species)
nb.min = min(vec.nb.words.species)
nb.max = max(vec.nb.words.species)

sR = c(vec.dist, dist.sum, dist.var, vec.nb.words.species, nb.mean, nb.var, nb.min, nb.max)





###


s = rbind(sA,sB)

indA = c()
for (i in 1:dim(sA)[1])
{
	indA = c(indA, "A")
}

indB = c()
for (i in 1:dim(sB)[1])
{
	indB = c(indB, "B")
}

ind = c(indA, indB)

install.packages("nnet")
library(abc)

cv = cv4postpr(ind, s, nval=100, tols=c(0.5, 0.1, 0.01), method = "neuralnet")

summary(cv)

ppr = postpr(sR, ind, s, tol=0.5, method = "neuralnet")

summary(ppr)





