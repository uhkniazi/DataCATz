# Name: graphWeights.R
# Auth: umar niazi
# Date: 10/09/2019
# Desc: weighting graph observed to expected probabilities

library(igraph)
setwd('countDataPoissonOthers/')
######### simulate data to create a graph
#########
p.old = par()
set.seed(123)
# universe of possible type 2 vertices
type.2.universe = LETTERS[1:26]
type.1.universe = 1:7
# graph data frame
dfGraph = NULL
# randomly assign labels to type 1 vertices
for (i in seq_along(type.1.universe)){
  s = sample(type.2.universe, runif(1, 1, length(type.2.universe)), replace = F)
  df = data.frame(i, s)
  dfGraph = rbind(dfGraph, df)
}
## assign some labels non randomly
i = 8:10
for (x in sample(1:26, 5, replace = F)){
  s = LETTERS[x]
  df = data.frame(i, s)
  dfGraph = rbind(dfGraph, df)
}
########
####### create bipartite graph
oIGbp = graph.data.frame(dfGraph, directed = F)
# set the vertex type variable to make graph bipartite
f = rep(c(T, F), times = c(length(unique(dfGraph[,1])),length(unique(dfGraph[,2]))))
V(oIGbp)$type = f
# sanity check - is graph bipartite
if (!is.bipartite(oIGbp)) {
  stop(paste('Graph is not bipartite'))
}

# make the type 2 vertices square
fType = V(oIGbp)$type
V(oIGbp)[fType]$shape = 'circle'
V(oIGbp)[!fType]$shape = 'square'

plot(oIGbp, layout=layout_as_bipartite, vertex.size=10)
######## end simulate data

####### calculate probabilities and project graph
#######
# vertex of the first kind will be assigned probabilities
# based on their relations with the vertices of the second kind
# flag to identify vertex types
f = V(oIGbp)$type
d = degree(oIGbp)
d = d[f]
# r is the total numbers of vertices of the second kind
r = sum(!f)
p = d/r
V(oIGbp)[f]$prob_marginal = p
bVertexType = f

# project the graph in one dimension and
# assign weights based on observed to expected ratios
g.p = bipartite.projection(oIGbp, which = 'TRUE')
# get the matrix with rows representing each edge
m = get.edgelist(g.p)
w = E(g.p)$weight
# calculate observed ratio
# weight / r
ob = w / r
# calculate expected 
mExp = cbind(V(g.p)[m[,1]]$prob_marginal, V(g.p)[m[,2]]$prob_marginal)
ex = mExp[,1] * mExp[,2]
# E(g.p)$observed = ob
# E(g.p)$expected = ex
# E(g.p)$ob_to_ex = ob / ex
# obj@ig.p = g.p
# return(obj)
# #######
####### end projection

dfData = data.frame(obs=w, ex=r*ex)

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

stanDso = rstan::stan_model(file='negBinomialmultiPhi.stan')

lStanData = list(Ntotal=nrow(dfData),
                 #mu=dfData$ex, 
                 y=dfData$ob)

fit.stan = sampling(stanDso, data=lStanData, iter=1000, chains=2,
                    cores=2)
print(fit.stan, digits=3)
