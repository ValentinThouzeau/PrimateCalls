# Checking for phylogenetic signal in the data

# Load necessary libraries
library(ape)          # For working with phylogenetic trees
library(phytools)     # For calculating phylogenetic signal
library(geiger)       # For additional phylogenetic tools
library(FactoMineR)   # For performing PCA

set.seed(1989)

# Load phylogenetic tree
tree <- read.tree("primatesx100.tre")  # Replace with your tree file path

# Load data
dat1 = read.nexus.data("lexicons_primates_food_binary.nex")
dat1 = t(data.frame(dat1))
mode(dat1)="integer"
data1 = as.data.frame(dat1)

dat2 = read.nexus.data("lexicons_primates_alarm_binary.nex")
dat2 = t(data.frame(dat2))
mode(dat2)="integer"
data2 = as.data.frame(dat2)

dat3 = read.nexus.data("lexicons_primates_non_disturb_binary.nex")
dat3 = t(data.frame(dat3))
mode(dat3)="integer"
data3 = as.data.frame(dat3)

# check that rows are in the same order
all(rownames(data1) == rownames(data2))
all(rownames(data1) == rownames(data3))

data = cbind(data1, data2, data3)


# Ensure the species in the data match the tree's tip labels
data <- data[match(tree$tip.label, rownames(data)), ]

# Perform PCA on binary data
# Using FactoMineR, scale.unit = FALSE prevents scaling of binary variables
pca_result <- PCA(data, graph = FALSE, scale.unit = FALSE)

# Extract principal components (e.g., first 2 PCs)
pc1 <- pca_result$ind$coord[,1]

# Quantify phylogenetic signal for PC1 using Blomberg's K
phylosig(tree, pc1, method = "K", test = TRUE)


