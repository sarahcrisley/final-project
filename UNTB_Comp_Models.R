#Ecology & Evolutionary Theory Final Project
#Name: Phoebe Jekielek and Sarah Risley


##Installing roleLite
#devtools::install_github('role-model/roleLite')

library(roleLite)


########################
###NEUTRAL SIMULATION###
########################

##Run simulation
neut <- untbSim(J = 1000, # number of individuals in local community
                Sm = 50, # number of species in meta community
                Jm = 10000, # number of individuals in meta community
                nu = 0, # speciation probability
                m = 0.1, # immigration probability 
                niter = 10000) # number of iteration to run

##Passing to getAbund
neutAbund <- getAbund(neut)
#log-series in meta community, initialize with random pull form meta comm for start of local comm

# returns a data.frame, see `?getAbund` for more detials
head(neutAbund)



############################
###COMPETITIVE SIMULATION###
############################

S <- 50 # number of species
#best practice to establish variable outside of fxn call 

# make a competition matrix
a <- matrix(runif(S^2, 0, 5), nrow = S) 
diag(a) <- runif(S, 1, 10) #increasing avg niche differentiation 

# run sim
comp <- compSim(J = 1000, # number of individuals in local community
                Sm = S, # number of species in meta community
                Jm = 10000, # number of individuals in meta community
                nu = 0, # speciation probability
                m = 0.1, # immigration probability
                alpha = a, # competition matrix
                niter = 10000) # number of iteration to run

compAbund <- getAbund(comp)
head(compAbund)



#######################
###SAD VISUALIZATION###
#######################

##Maxes
max(neutAbund$abund) #277
max(compAbund$abund) #262


##Looking at the final SAD for neut
finalSADneut <- neutAbund$abund[neutAbund$tstep == max(neutAbund$tstep)]

plot(sort(finalSADneut, decreasing = TRUE), log = 'y', 
     xlab = 'Rank', ylab = 'Abundance',
     main = 'Neutral SAD')



##Looking at the final SAD for comp
finalSADComp <- compAbund$abund[compAbund$tstep == max(compAbund$tstep)]

plot(sort(finalSADComp, decreasing = TRUE), log = 'y', 
     xlab = 'Rank', ylab = 'Abundance',
     main = 'Competitive SAD')


##Both
#par(mfrow=c(2,1))
plot(sort(finalSADneut, decreasing = TRUE), log = 'y', 
     xlab = 'Rank', ylab = 'Abundance',
     main = 'Neutral SAD', col = 'blue',
     xlim = c(0,30))
plot(sort(finalSADComp, decreasing = TRUE), log = 'y', 
     xlab = 'Rank', ylab = 'Abundance',
     main = 'Competitive SAD', col = 'red',
     xlim = c(0,30))


##################
###PURTUBATIONS###
##################

##Neutral: rare species becomes abundant
# parameter set-up
S <- 50 
J <- 1000
nu <- 0
m <- 0.1
niter <- 10000


# make meta SAD outside of `untbSim` so we can manipulate it to reflect 
# invasion...hang tight, you'll see
metaSAD <- pika::rlseries(S, 0.001) 

# run sim pre-invasion
untbPre <- untbSim(J = J, nu = nu, m = m, niter = niter, initMeta = metaSAD)

#don't need to specify Sm and Jm 

##Modifying metacommunity so least common species is 2x's as common
newMetaSAD <- metaSAD
irare <- which.min(newMetaSAD)[1] # index of a rare species
newMetaSAD[irare] <- 2 * max(newMetaSAD)

# now we can pass this new meta SAD to the sim function, and pass the final 
# state of the pre-invasion simulation as the initial condition for this sim:
untbPost <- untbSim(J = J, nu = nu, m = m, niter = niter, 
                    initLocal = untbPre[nrow(untbPre), ], 
                    initMeta = newMetaSAD)



##Combining the two simulations
abundPre <- getAbund(untbPre)
abundPost <- getAbund(untbPost)

# note: first timestep of post is last of pre, so we can discard one
abundPost <- abundPost[abundPost$tstep != 0, ]

# update timestep values of post to come *after* pre
abundPost$tstep <- abundPost$tstep + max(abundPre$tstep)

# combine
abundAll <- rbind(abundPre, abundPost)

# visualize abundance timeseries for invader
abundInvader <- abundAll[abundAll$spID == irare, ]

plot(abundInvader$tstep, abundInvader$abund, type = 'b',
     xlab = 'Timestep', ylab = 'Abundance')




###################
###INVASION SADS###
###################

##Maxes
max(abundPre$abund[abundPre$tstep == max(abundPre$tstep)]) #431
max(abundPost$abund[abundPost$tstep == 13000]) #409
max(abundPost$abund[abundPost$tstep == max(abundPost$tstep)]) #354


#Right before invasion
plot(sort(abundPre$abund[abundPre$tstep == max(abundPre$tstep)], 
          decreasing = TRUE), 
     log = 'y', 
     xlab = 'Rank', ylab = 'Abundance', main = 'Right before invasion')


#Shortly after invasion
plot(sort(abundPost$abund[abundPost$tstep == 13000], 
          decreasing = TRUE), 
     log = 'y', 
     xlab = 'Rank', ylab = 'Abundance', main = 'Shortly after invasion')


#Long after invasion
plot(sort(abundPost$abund[abundPost$tstep == max(abundPost$tstep)], 
          decreasing = TRUE), 
     log = 'y', 
     xlab = 'Rank', ylab = 'Abundance', main = 'Long after invasion')


##Getting final SAD's for time steps 
preSADinvas <- abundPre$abund[abundPre$tstep == max(abundPre$tstep)]
postSADinvas <- abundPost$abund[abundPost$tstep == 13000]
longSADinvas <- abundPost$abund[abundPost$tstep == max(abundPost$tstep)]


##################
###HILL NUMBERS###
##################

#hill package 
library(hillR)

# Evenness 
q <- 2 #Simpsons evenness 

neut_2 <- hillR::hill_taxa(finalSADneut, q)
comp_2 <- hillR::hill_taxa(finalSADComp, q)
pre_2 <- hillR::hill_taxa(preSADinvas, q)
post_2 <- hillR::hill_taxa(postSADinvas, q)
long_2 <- hillR::hill_taxa(longSADinvas, q)

q <- 0 #Species richness 

neut_0 <- hillR::hill_taxa(finalSADneut, q)
comp_0 <- hillR::hill_taxa(finalSADComp, q)
pre_0 <- hillR::hill_taxa(preSADinvas, q)
post_0 <- hillR::hill_taxa(postSADinvas, q)
long_0 <- hillR::hill_taxa(longSADinvas, q)


q <- 1-(1e-5) #Shannon's Index 
neut_1 <- hillR::hill_taxa(finalSADneut, q)
comp_1 <- hillR::hill_taxa(finalSADComp, q)
pre_1 <- hillR::hill_taxa(preSADinvas, q)
post_1 <- hillR::hill_taxa(postSADinvas, q)
long_1 <- hillR::hill_taxa(longSADinvas, q)

# q = 0, get species richness 
# q can't be exactly one, but as q gets close to 1 use 1-1E^-5 (like 1)

##Compiling 
tab <- matrix(c(neut_2, comp_2, pre_2, post_2,long_2,
                neut_0, comp_0, pre_0, post_0, long_0,
                neut_1, comp_1, pre_1, post_1, long_1), ncol=3, nrow = 5)
colnames(tab) <- c('Simpsons Index','Species Richness', 'Shannons Index')
rownames(tab) <- c('Neutral Model','Competitive Model', 'Pre-Invasion', 'Post-Invasion', 'After Invasion')
tab <- as.table(tab)
tab <- round(tab, digits = 2)
tab

