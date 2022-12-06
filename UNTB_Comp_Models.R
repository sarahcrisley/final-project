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
#?Is this log-series? What does the untbSim simulate as the original pop?

# returns a data.frame, see `?getAbund` for more detials
head(neutAbund)



############################
###COMPETITIVE SIMULATION###
############################

S <- 50 # number of species
#? Why do we specify this variable here rather than in the simulation like above?

# make a competition matrix
a <- aa <- matrix(runif(S^2, 0, 5), nrow = S) #creates duplicate matrices
diag(aa) <- runif(S, 1, 10) #? Extracting the diagonal of the matrix aa? 


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
par(mfrow=c(2,1))
plot(sort(finalSADneut, decreasing = TRUE), log = 'y', 
     xlab = 'Rank', ylab = 'Abundance',
     main = 'Neutral SAD', col = 'blue',
     xlim = c(0,25))
plot(sort(finalSADComp, decreasing = TRUE), log = 'y', 
     xlab = 'Rank', ylab = 'Abundance',
     main = 'Competitive SAD', col = 'red',
     xlim = c(0,25))



##Time series for a specific species: neut
# look at population timeseries for species 10
sp10_neut <- neutAbund[neutAbund$spID == 10, ]

plot(sp10_neut$tstep, sp10_neut$abund, type = 'b', 
     xlab = 'Time step', ylab = 'Abundance')


##Time series for a specific species: comp
# look at population timeseries for species 10
sp10 <- compAbund[compAbund$spID == 10, ]

plot(sp10$tstep, sp10$abund, type = 'b', 
     xlab = 'Time step', ylab = 'Abundance')




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







