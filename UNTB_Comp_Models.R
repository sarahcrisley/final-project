#Ecolgy & Evolutionary Theory Final Project
#Name: Phoebe Jekielek and Sarah Risley


##Installing roleLite
devtools::install_github('role-model/roleLite')

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

# returns a data.frame, see `?getAbund` for more detials
head(neutAbund)

##   tstep spID abund
## 1     0    1   245
## 2     0    2   159
## 3     0    3   103
## 4     0    4    83
## 5     0    5    73
## 6     0    6    57





############################
###COMPETITIVE SIMULATION###
############################

S <- 50 # number of species

# make a competition matrix
a <- aa <- matrix(runif(S^2, 0, 5), nrow = S)
diag(aa) <- runif(S, 1, 10)

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

##   tstep spID abund
## 1     0    1   283
## 2     0    2   145
## 3     0    3   105
## 4     0    4    80
## 5     0    5    66
## 6     0    6    55






#######################
###SAD VISUALIZATION###
#######################

##Looking at the final SAD
finalSADComp <- compAbund$abund[compAbund$tstep == max(compAbund$tstep)]

plot(sort(finalSADComp, decreasing = TRUE), log = 'y', 
     xlab = 'Rank', ylab = 'Abundance')




##Time series for a specific species
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







