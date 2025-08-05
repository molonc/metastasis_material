library(igraph)
# Random bipartite graph
inc <- matrix(sample(0:1, 50, replace = TRUE, prob=c(2,1)), 10, 5)
g <- graph_from_incidence_matrix(inc)
plot(g, layout = layout_as_bipartite,
     vertex.color=c("green","cyan")[V(g)$type+1])

# Two columns
g %>%
  add_layout_(as_bipartite()) %>%
  plot()



library(igraph)
library(dplyr)
# Random bipartite graph
inc <- matrix(sample(0:1, 50, replace = TRUE, prob=c(2,1)), 10, 5)
g <- graph_from_incidence_matrix(inc)
plot(g, layout = layout_as_bipartite,
     vertex.color=c("green","cyan")[V(g)$type+1])

# Two columns
g %>%
  add_layout_(as_bipartite()) %>%
  plot()
LO = layout_as_bipartite(g)
LO = LO[,c(2,1)]
plot(g, layout = LO, vertex.color=c("green","cyan")[V(g)$type+1])


# https://stackoverflow.com/questions/65989711/arrange-nodes-at-specific-location
# https://stackoverflow.com/questions/31366066/how-to-plot-a-bipartite-graph-in-r
rm(g)
V(g)$x <- c(1, 1, 1, 2, 2, 2, 2)
V(g)$y <- c(3, 2, 1, 3.5, 2.5, 1.5, 0.5)
V(g)$shape <- shape[as.numeric(V(g)$type) + 1]
V(g)$color <- c('red', 'blue', 'green', 'steelblue', 'steelblue', 'steelblue', 'steelblue')
E(g)$color <- 'gray'
# E(g)$color[E(g)['A' %--% 'C']] <- 'yellow'
# E(g)$color[E(g)['B' %--% 'C']] <- 'purple'
E(g)$color[E(g)['A' %--% V(g)]] <- 'red'
E(g)$color[E(g)['B' %--% V(g)]] <- 'blue'
E(g)$color[E(g)['C' %--% V(g)]] <- 'green'
plot(g)


# https://rpubs.com/lgadar/load-bipartite-graph
A <- matrix(c(2,4,0,0,1,2,1,0,0,0,4,5,0,0,1,2,0,0,6,2), byrow = T, nrow = 5, ncol = 4)
rownames(A) <- letters[c(1:nrow(A))]
colnames(A) <- LETTERS[c(1:ncol(A))]
print(A)

g <- graph.incidence(A, weighted = T)
g
V(g)$type
V(g)$name

# g <- graph.empty(directed = F)
# node.out <- unique(edges$from) #stringsAsFactor = F in data frame
# node.in <- unique(edges$to) #stringsAsFactor = F in data frame
# g <- add.vertices(g,nv=length(node.out),attr=list(name=node.out),type=rep(FALSE,length(node.out)))
# g <- add.vertices(g,nv=length(node.in),attr=list(name=node.in),type=rep(TRUE,length(node.in)))
# edge.list.vec <- as.vector(t(as.matrix(data.frame(edges))))
# g <- add.edges(g,edge.list.vec)
# g
# g <- graph.data.frame(edges, directed = F)
# V(g)$type <- V(g)$name %in% edges[,2] #the second column of edges is TRUE type
# g


A <- matrix(c(2,4,0,0,1,2,1,0,0,0,4,5,0,0,1,2,0,0,6,2), byrow = T, nrow = 5, ncol = 4)
rownames(A) <- letters[c(1:nrow(A))]
colnames(A) <- LETTERS[c(1:ncol(A))]
print(A)

g <- graph.incidence(A, weighted = T)

get.incidence(g)
V(g)$color <- V(g)$type
V(g)$color=gsub("FALSE","red",V(g)$color)
V(g)$color=gsub("TRUE","blue",V(g)$color)
# plot(g, edge.color="gray30",edge.width=E(g)$weight, layout=layout_as_bipartite)
E(g)$color[E(g)['D' %--% 'e']] <- 'yellow'
E(g)$color[E(g)['B' %--% 'a']] <- 'purple'

LO = layout_as_bipartite(g)
LO = LO[,c(2,1)]

plot(g, layout = LO,edge.width=E(g)$weight)
