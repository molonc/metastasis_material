library(ape)
library(phytools)
text.string<-
  "(((((((cow, pig),whale),(bat,(lemur,human))),(robin,iguana)),coelacanth),gold_fish),shark);"
tree<-read.tree(text=text.string)
plotTree(tree,offset=1)
tiplabels() # add number to tip labels
nodelabels() # add internal nodes labels

str(tree) # tree structure
tree$edge
# E(tree)

tree$tip.label
Ntip(tree)

tree$Nnode # number of interior nodes in the tree


plotTree(tree,ftype="i",fsize=0.6,lwd=1)
plotTree(anolis.tree,type="fan",fsize=0.7,lwd=1,
         ftype="i")
add.arrow(anolis.tree,tip=ii,arrl=1)


anolis.noPR<-drop.tip(anolis.tree,pr.species)
plotTree(anolis.noPR,type="fan",fsize=0.7,lwd=1,
         ftype="i")


anolis.tree<-read.tree(file="/home/htran/Projects/hakwoo_project/corrupt_tree/src/corrupt_grow/anole.tre")
anolis.tree
pr.species<-c("cooki","poncensis","gundlachi","pulchellus","stratulus",
              "krugi","evermanni","occultus","cuvieri","cristatellus")
ii<-sapply(pr.species,grep,anolis.tree$tip.label)
ii
pr.species<-anolis.tree$tip.label[ii]
pr.species
pr.tree<-drop.tip(anolis.tree, setdiff(anolis.tree$tip.label,pr.species))
plotTree(pr.tree,ftype="i")

anolis.noPR<-drop.tip(anolis.tree,pr.species)
plotTree(anolis.noPR,type="fan",fsize=0.7,lwd=1,
         ftype="i")

node<-fastMRCA(anolis.tree,"Anolis_evermanni",
               "Anolis_cristatellus")
pr.clade<-extract.clade(anolis.tree,node)
plotTree(pr.clade,ftype="i")
?fastMRCA

anolis.pruned<-collapseTree(anolis.tree)
plotTree(anolis.pruned,type="fan",fsize=0.7,lwd=1,
         ftype="i")
write.tree(tree,file="example.tre")
## this is what our file looks like (you can open it to check)
cat(readLines("example.tre"),sep="\n")


## ok, now let's re-root the tree at node #67
rr.67<-root(tree,node=67)
plotTree(rr.67)

## check if tree & rt.all are equal
all.equal(tree,rt.all)



set.seed(123)## generate tree
tree<-pbtree(b=0.1, n=10)
## plot original tree
plot(tree)
axisPhylo()
## add an extant tip ("t_extant") sister to taxon't5'## with divergence time of 4.5 Ma
node <- which(tree$tip.label=="t5")
tiplabels() # add number to tip labels
nodelabels() # add internal nodes labels

tree <- bind.tip(tree, tip.label="t_extant",where=node, position=1) #, position=4.5
for(i in 1:10){
  # print(i)
  tree <- bind.tip(tree, tip.label=paste0("t_extant_",i),where=node) #, position=4.5
}
tree <- bind.tip(tree, tip.label=paste0("t_extant_",11),where=node) #, position=4.5

tree <- drop.tip(tree, "t5")
?drop.tip

# plot to see the result
plot(tree)

axisPhylo()
## add an extinct tip ("t_extinct") sister to't2'with## divergence time of 7.8 Ma and duration (edge length) of## 3.3 Ma
node <- which(tree$tip.label=="t2")
tree <- bind.tip(tree, tip.label="t_extinct", where=node,position=7.8, edge.length=3.3)
## plot to see the result
plot(tree)
axisPhylo()


# https://biology.stackexchange.com/questions/37057/construct-phylogeny-form-edgelist
phylo_from_el <- function(el){
  i <- 0
  for(n in el[2:length(el)]){
    i <- i+1
    parents[i] <- n 
    children[i] <- i+1
  }
  nnode <- sum(children %in% parents) + 1
  nleaf <- length(el)-nnode
  # Translate preorder sequence into ape sequence
  ttable <- matrix(c(1:length(el), rep(0, length(el))), ncol=2)
  ttable[unique(parents),2] <- 1:nnode + nleaf
  ttable[-unique(parents),2] <- 1:nleaf
  
  parents_t <- sapply(parents, FUN = function(x)ttable[x,2])
  children_t <- sapply(children, FUN = function(x)ttable[x,2])
  
  edge <- matrix(c(parents_t, children_t), ncol = 2)
  tr <- list(edge = edge, tip.label = 1:nleaf, Nnode = nnode)
  class(tr) <- "phylo"
  return(tr)
}
children <- c(1, 2, 3, 4)
parents <- c(1, 1, 2, 2)
el <- c(1,1,1,3,3)

tr <- phylo_from_el(el)

plot(tr)

