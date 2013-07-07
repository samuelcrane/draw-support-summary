# labelNodeSupport: map node support values from parsimony, likelihood, and Bayesian phylogenies onto 
# a single target tree as filled circles. This script uses some bits of code from various people, 
# including Rich Glor, Luke Harmon, and Emmanuel Paradis. 
# For more info on the filled circle plotting, see https://github.com/samuelcrane/draw-node-support-symbol 
# For more info on the subclade comparison and node support mapping, see https://github.com/samuelcrane/label-node-support

# Load necessary libraries
library(ape)
library(geiger)

# The subfunction that atomizes a tree into each individual subclade and was provided compliments of Luke Harmon.
getAllSubtrees <- function(phy, minSize=2) 
{
    res <- list() 
    count = 1 
    ntip <- length(phy$tip.label) 
    for(i in 1:phy$Nnode) 
    { 
        l <- tips(phy, ntip+i) 
        bt <- match(phy$tip.label, l) 
        if(sum(is.na(bt)) == 0) 
        {
            st <- phy 
        } 
        else st <- drop.tip(phy, phy$tip.label[is.na(bt)]) 
        if(length(st$tip.label)>=minSize) 
        { 
            res[[count]] <- st 
            count <- count+1 
        }
    } 
    res
}

# The plotSupportSummary function will add a grey-filled circle on the maximum likelihood tree if at 
# least one of the ML, parsimony, or Bayesian trees supports that node and will draw a black-filled
# circle if all three do. Parts of this come from label-node-support and parts from draw-node-support-symbol 
plotSupportSummary <- function(targetTree, suppTree1, suppTree2) 
{
    targetTree <- ladderize(root(targetTree, "Outgroup_species"))
    getAllSubtrees(targetTree) -> targetSub
    getAllSubtrees(suppTree1) -> suppSub1
    getAllSubtrees(suppTree2) -> suppSub2
    suppList1 <- matrix("", Nnode(targetTree), 1)
    suppList2 <- matrix("", Nnode(targetTree), 1)
    
    #The commands below compare all the subclades in the Bayes tree to all the subclades in the bootstrap tree, and vice versa, and identifies all those clades that are identical.
    for(i in 1:Nnode(targetTree)) 
    {
        for(j in 1:Nnode(suppTree1)) 
        {
            match(targetSub[[i]]$tip.label[order(targetSub[[i]]$tip.label)], suppSub1[[j]]$tip.label[order(suppSub1[[j]]$tip.label)]) -> shared
            match(suppSub1[[j]]$tip.label[order(suppSub1[[j]]$tip.label)], targetSub[[i]]$tip.label[order(targetSub[[i]]$tip.label)]) -> shared2
            if(sum(is.na(c(shared,shared2)))==0) 
            {
                round(as.numeric(suppTree1$node.label[j]), digits=2) -> suppList1[i]
            }
        }
        for(j in 1:Nnode(suppTree2)) 
        {
            match(targetSub[[i]]$tip.label[order(targetSub[[i]]$tip.label)], suppSub2[[j]]$tip.label[order(suppSub2[[j]]$tip.label)]) -> shared
            match(suppSub2[[j]]$tip.label[order(suppSub2[[j]]$tip.label)], targetSub[[i]]$tip.label[order(targetSub[[i]]$tip.label)]) -> shared2
            if(sum(is.na(c(shared,shared2)))==0) 
            {
                suppTree2$node.label[j] -> suppList2[i]
            }
        }
    }
    plot(targetTree, cex=0.5, lwd=0.5, direction='r', use.edge.length=FALSE, label.offset=1, no.margin=TRUE, x.lim=c(0.0000, 300), y.lim=c(0, 110)) 

    # Color palette for filled circles
    co <- c("black", "grey", "white")

    # Initilize character matrix for drawing node circles
    p <- character(length(targetTree$node.label))
    p[] <- co[3] # Set all nodes to white
    p[(as.numeric(targetTree$node.label) >= 70) | (suppList1 >= 0.90) | (suppList2 >= 70)] <- co[2] # If at least one tree provides good support, set node to grey
    p[(as.numeric(targetTree$node.label) >= 70) & (suppList1 >= 0.90) & (suppList2 >= 70)] <- co[1] # If all three trees provide good support, set node to black
    
    # This for loop draws filled circles on only those nodes with good support from at least one method. 
    for(j in 1:Nnode(targetTree)) 
        {
            if(targetTree$node.label[[j]] != "" & p[j] != "white") 
            {
                nodelabels("", j+length(targetTree$tip.label), cex=0.75, bg=p[j], pch=21, frame="n")
            }
        }
}

# Read in the trees
bayesTree <- read.nexus("MrBayes.tre")
mlTree <- read.tree("RAxML.tre")
mpTree <- read.nexus("TNT.tre")

# Call the ploting function
plotSupportSummary(mlTree, bayesTree, mpTree)
