#new analysis
install.packages("C:/Users/Ambrose/Desktop/subgraphMining_1.0.tar.gz", repos =NULL, type = "source")
install.packages("C:/Users/Ambrose/Desktop/igraph0_0.5.7.tar.gz", repos =NULL, type = "source")
install.packages("C:/Users/Ambrose/Desktop/Courses/Bana 7047/LinkAnalysis.zip", repos = NULL)
setwd("C:/Users/Ambrose/Desktop/")


library(igraph)
g <- random.graph.game(20, 5/20, directed=TRUE)
page.rank(g)$vector

plot(g)


g2 <- graph.star(10)
page.rank(g2)$vector
plot(g2)

library(subgraphMining)
library(igraph0)
# The cslogs data set
data(cslogs)
# The matabolicInteractions data
graph = data(metabolicInteractions)
summary(graph)
plot(lc, type = "summary")

data = data(metabolicInteractions)
gspan(data, support="80%")
str(data)


library(linkcomm)

demo(topic = "linkcomm", package = "linkcomm")
A set of 56 yeast proteins involved in 651 interactions related to transcription.


forest.fire.game (nodes, fw.prob, bw.factor = 1, ambs = 1, directed = TRUE)

g <- forest.fire.game(10000, fw.prob=0.37, bw.factor=0.32/0.37)
dd1 <- degree.distribution(g, mode="in")
dd2 <- degree.distribution(g, mode="out")
if (interactive()) {
  plot(seq(along=dd1)-1, dd1, log="xy")
  points(seq(along=dd2)-1, dd2, col=2, pch=2)
}

#Yeast Analysis
yeast_pp <- pp_rnapol
yeast_pp <- human_pp
head(yeast_pp)

lc <- getLinkCommunities(yeast_pp, hcmethod = "single")
class(lc)
print(lc)

plot(lc, type = "graph", layout = layout.fruchterman.reingold)
plot(lc, type = "graph", layout = "spencer.circle")
plot(lc, type = "graph", shownodesin = 2, node.pies = TRUE)
plot(lc, type = "members")

plot(lc, type = "summary")
plot(lc, type = "dend")
getAllNestedComm(lc)

getNestedHierarchies(lc, clusid = 9)
#or
plot(lc, type = "graph", clusterids = c(9,11))
cr <- getClusterRelatedness(lc, hcmethod = "ward.D")
	


Ambrose <- read.graph("C:/Users/Ambrose/Desktop/Courses/Bana 7047/Project/data/Facebook.gml",format=c("gml"))
E(Ambrose)
V(Ambrose)$label
V(Ambrose)[ sex == "male" ]
plot(degree.distribution(Ambrose, mode="in"), log="xy")
plot(degree.distribution(Ambrose, mode="out"), log="xy")

library(intergraph)
Ambrose.g = as.matrix(Ambrose, matrix.type = c("adjacency", "edgelist"))
class(g_ambrose)
class(Ambrose)

Ambrose.el <- get.edgelist(Ambrose)

Ambrose.lc <- getLinkCommunities(Ambrose.el,directed = FALSE,  hcmethod = "single")

plot(Ambrose.lc, type = "members")

plot(Ambrose.lc, type = "dend")
plot(Ambrose.lc, type = "summary")


#display nodes that belong to 3 or more communities
plot(Ambrose.lc, type = "graph", layout = "spencer.circle", shownodesin = 5)
plot(Ambrose.lc, type = "graph", shownodesin = 2, node.pies = TRUE )

#community analysis
getAllNestedComm(Ambrose.lc)
plot(Ambrose.lc, type = "graph", clusterids = c(25,24)	)

cr <- getClusterRelatedness(Ambrose.lc, hcmethod = "ward.D")
cr <- getClusterRelatedness(Ambrose.lc, hcmethod = "ward.D2")


cm <- getCommunityConnectedness(Ambrose.lc, conn = "modularity")
plotLinkCommSummComm(Ambrose.lc, type = "commsumm", summary = "modularity")

#modularity
fgreedy<-fastgreedy.community(Ambrose,merges=TRUE, modularity=TRUE)
memberships <-community.to.membership(Ambrose, fgreedy$merges, steps=which.max(fgreedy$modularity)-1)
print(paste('Number of detected communities=',length(memberships$csize)))
# Community sizes:
print(memberships$csize)
# modularity:
max(fgreedy$modularity)

plot(fgreedy, Ambrose, vertex.size=0.1, vertex.label.font=0.01)
plot(Ambrose, vertex.color=membership(fgreedy),vertex.size=2, vertex.label.font=0.2, vertex.label.color=membership(fgreedy))

dendPlot(fgreedy)


wc = edge.betweenness.community(Ambrose)
modularity(wc)
plot(Ambrose, vertex.color=membership(wc),vertex.size=0.8, vertex.label.font=0.2, vertex.label.color=membership(wc) )

V(Ambrose)$color=clusters(g2)$membership

library(intergraph)
#more conversion
net <- asNetwork(Ambrose)
net2 <- asDF(net)
summary(net2)
cbind(net2$vertexes["id"], net2$vertexes["label"])

net2$edges
plot(net)

plot.igraph(Ambrose,vertex.frame.color='blue', vertex.label.font=0.5, vertex.label=V(Ambrose)$label	)

class(lc)


lc <- getLinkCommunities(net,directed = FALSE,  hcmethod = "single")

Yeast = read.table("C:/Users/Ambrose/Desktop/Courses/Bana 7047/Project/data/yeast.txt")
lc_yeast <- getLinkCommunities(Yeast,directed = FALSE,  hcmethod = "single")
plot(lc_yeast, type = "graph", layout = "spencer.circle", shownodesin = 5)
plot(lc_yeast, type = "members")
getAllNestedComm(lc_yeast)
plot(degree.distribution(Yeast, mode="in"), log="xy")






