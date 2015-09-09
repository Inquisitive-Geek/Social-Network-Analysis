install.packages("statnet")
install.packages("igraph")
install.packages("intergraph")
library(statnet)
library(igraph)
library(intergraph)

#Reading the graph
fbAarti <- read.graph(file = "Aart_FB.gml", format = "gml")

#Converting the file to a network file
fb <- asNetwork(fbAarti)
plot(fb)
summary(fb)

#Mixing Matrices based on attribute names
mixingmatrix(fb, "locale")
mixingmatrix(fb, "agerank")
mixingmatrix(fb, "label")
mixingmatrix(fb, "sex")

#Network components
table(component.dist(fb)$csize)
network.size(fb)

#Erdos-Renyi model
fb.model.er <- ergm(fb~edges)
summary(fb.model.er)
fb.model.er$mle.lik

#Assortative mixing edge model
fb.model.assor <- ergm(fb~edges + nodematch("locale") + nodematch("agerank") + nodematch("sex"))
summary(fb.model.assor)
fb.model.assor$mle.lik

fb.model.assor.2 <- ergm(fb~edges + nodematch("locale") + nodematch("sex"))
summary(fb.model.assor.2)
fb.model.assor.2$mle.lik

#Simulate the network
sim.assor <- simulate(fb.model.assor.2, burnin = 1e+6, verbose = TRUE, seed = 9)
mixingmatrix(sim.assor, "locale")
mixingmatrix(sim.assor, "sex")

#Network components
table(component.dist(fb)$csize)
network.size(fb)

#Erdos-Renyi model
fb.model.er <- ergm(fb~edges)
summary(fb.model.er)
fb.model.er$mle.lik

#Assortative mixing edge model
fb.model.assor <- ergm(fb~edges + nodematch("locale") + nodematch("agerank") + nodematch("sex"))
summary(fb.model.assor)
fb.model.assor$mle.lik

fb.model.assor.2 <- ergm(fb~edges + nodematch("locale") + nodematch("sex"))
summary(fb.model.assor.2)
fb.model.assor.2$mle.lik

#Simulate the network
sim.assor <- simulate(fb.model.assor.2, burnin = 1e+6, verbose = TRUE, seed = 9)
mixingmatrix(sim.assor, "locale")
mixingmatrix(sim.assor, "sex")

#Degree distribution
summary(fb ~ degree(0:120))
summary(sim.assor ~ degree(0:120))


c(fb = summary(fb ~ triangle), sim.assor = summary(sim.assor ~ triangle))

fb.model.ed.tr <- ergm(fb ~ edges + triangle, maxit = 4, control.ergm =  control.ergm(seed = 9), verbose = TRUE)

summary(fb.model.ed.tr)
fb.model.ed.tr$mle.lik

mcmc.diagnostics(fb.model.ed.tr)
gof.fb.model.ed.tr <- gof(fb.model.ed.tr)
plot(gof.fb.model.ed.tr)


fb.model.ed.tr.1 <- ergm(fb ~ edges + triangle + nodematch("sex") + nodematch("locale") + degree(1),
                       maxit = 4, MCMCsamplesize = 20000, control.ergm =  control.ergm(seed = 9, MCMC.interval=500), 
                       calc.mcmc.se = FALSE, parallel = 4, verbose = TRUE)




gof.deg.fb.model.ed.tr <- gof(fb.model.ed.tr~degree)
plot(gof.deg.fb.model.ed.tr)
#Simulate the network
sim.ed.tr <- simulate(fb.model.ed.tr, burnin = 1e+6, verbose = TRUE, seed = 9, nsim = 5)

mixingmatrix(sim.ed.tr[[2]], "sex")
mixingmatrix(sim.ed.tr[[2]], "locale")

# GWESP model with assortative fields
fb.gwesp <- ergm(fb ~ edges + nodematch("sex") + nodematch("locale") + gwesp(0.25, fixed = TRUE), 
               MCMCsamplesize = 10000, maxit = 4, verbose = TRUE, 
               control = control.ergm(steplength = 0.25), seed = 9, parallel = 4)

#Comparing degree distrbution of Observed and Simulated network
plot(summary(sim.assor ~ degree(0:120)), type = "l", lty = 1, lwd = 2, xlab = "Degree", ylab = "Count", col = "black")
lines(summary(fb ~ degree(0:120)), lty = 1, lwd = 3, col = "blue")
lines(summary(sim.ed.tr[[1]] ~ degree(0:120)), lty = 1, lwd = 2, col = "red")
legend("topright", legend = c("Observed", "Simulated Assortative Edge", "Simulated Edge Triangle"), lwd = 3, col = c("black", "blue", "red"))

