#Ecology & Evolutionary Theory Final Project
#Name: Phoebe Jekielek and Sarah Risley


##Installing roleLite
#devtools::install_github('role-model/roleLite')

##Libraries
library(roleLite)
library(vegan)
library(pika)
library(tidyverse)


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

# params set-up
Sm <- 50
sadMeta <- rlseries(Sm, 0.001)
J <- 1000
nu <- 0
m <- 0.1
niter <- 10000


## Neutral sim hill 
hillz <- replicate(500, # number of times to replicate
                   { # put all the code you need to replicate inside the { }
                     oneNeut <- untbSim(J = J, nu = nu, m = m, niter = niter, 
                                        initMeta = sadMeta)
                     
                     oneFinalSAD <- getAbund(oneNeut, tstep = nrow(oneNeut))
                     
                     # hill numbers
                     theseHill <- renyi(oneFinalSAD$abund, scales = 0:2, 
                                        hill = TRUE)
                     
                     # last line inside the {} should be the output you want
                     return(theseHill)
                   })

# Summarizing this raw simulation, mean and SD are good options
neut_hillMeans <- apply(hillz, 1, mean)
neut_hillSDs <- apply(hillz, 1, sd)



## Comp sim hill 
hillz_2 <- replicate(500, # number of times to replicate
                   { # put all the code you need to replicate inside the { }
                           oneComp <- compSim(J = J, Sm = Sm, Jm = 10000, 
                                              nu = nu,m = m, alpha = a, niter = niter) 
                           
                           oneFinalSAD <- getAbund(oneComp, tstep = nrow(oneComp))
                           
                           # hill numbers
                           theseHill <- renyi(oneFinalSAD$abund, scales = 0:2, 
                                              hill = TRUE)
                           
                           # last line inside the {} should be the output you want
                           return(theseHill)
                   })


# Summarizing this raw simulation, mean and SD are good options
comp_hillMeans <- apply(hillz_2, 1, mean)
comp_hillSDs <- apply(hillz_2, 1, sd)



## Invasion sim hill 

# All
hillz_3 <- replicate(500, # number of times to replicate
                     { # put all the code you need to replicate inside the { }
                             metaSAD <- pika::rlseries(S, 0.001) 
                             untbPre <- untbSim(J = J, nu = nu, m = m, niter = niter,
                                                initMeta = metaSAD) # pre invasion
                             
                             newMetaSAD <- metaSAD # modifying meta community
                             irare <- which.min(newMetaSAD)[1] # index of a rare species
                             newMetaSAD[irare] <- 2 * max(newMetaSAD) # amplifying rare x2
                             
                             untbPost <- untbSim(J = J, nu = nu, m = m, niter = niter, 
                                                 initLocal = untbPre[nrow(untbPre), ], 
                                                 initMeta = newMetaSAD) # post invasion
                             
                             abundPre <- getAbund(untbPre) # combining sims
                             abundPost <- getAbund(untbPost)
                             abundPost <- abundPost[abundPost$tstep != 0, ] # removing overlap
                             abundPost$tstep <- abundPost$tstep + max(abundPre$tstep) # timesteps post after pre
                             abundAll <- rbind(abundPre, abundPost) # combining
                             
                             oneFinalSAD <- abundAll$abund[abundAll$tstep == max(abundAll$tstep)]
                
                             
                             # hill numbers
                             theseHill <- renyi(oneFinalSAD, scales = 0:2, 
                                                hill = TRUE)
                             
                             # last line inside the {} should be the output you want
                             return(theseHill)
                            
                     })

# Summarizing this raw simulation, mean and SD are good options
allinvas_hillMeans <- apply(hillz_3, 1, mean)
allinvas_hillSDs <- apply(hillz_3, 1, sd)




# Pre
hillz_4 <- replicate(500, # number of times to replicate
                     { # put all the code you need to replicate inside the { }
                             metaSAD <- pika::rlseries(S, 0.001) 
                             untbPre <- untbSim(J = J, nu = nu, m = m, niter = niter,
                                                initMeta = metaSAD) # pre invasion
                             
                             newMetaSAD <- metaSAD # modifying meta community
                             irare <- which.min(newMetaSAD)[1] # index of a rare species
                             newMetaSAD[irare] <- 2 * max(newMetaSAD) # amplifying rare x2
                             
                             untbPost <- untbSim(J = J, nu = nu, m = m, niter = niter, 
                                                 initLocal = untbPre[nrow(untbPre), ], 
                                                 initMeta = newMetaSAD) # post invasion
                             
                             abundPre <- getAbund(untbPre) # combining sims
                             abundPost <- getAbund(untbPost)
                             abundPost <- abundPost[abundPost$tstep != 0, ] # removing overlap
                             abundPost$tstep <- abundPost$tstep + max(abundPre$tstep) # timesteps post after pre
                             abundAll <- rbind(abundPre, abundPost) # combining
                             
                             oneFinalSAD <- abundPre$abund[abundPre$tstep == max(abundPre$tstep)]
                             
                             
                             # hill numbers
                             theseHill <- renyi(oneFinalSAD, scales = 0:2, 
                                                hill = TRUE)
                             
                             # last line inside the {} should be the output you want
                             return(theseHill)
                             
                     })

# Summarizing this raw simulation, mean and SD are good options
preinvas_hillMeans <- apply(hillz_4, 1, mean)
preinvas_hillSDs <- apply(hillz_4, 1, sd)




# Short-term post
hillz_5 <- replicate(500, # number of times to replicate
                     { # put all the code you need to replicate inside the { }
                             metaSAD <- pika::rlseries(S, 0.001) 
                             untbPre <- untbSim(J = J, nu = nu, m = m, niter = niter,
                                                initMeta = metaSAD) # pre invasion
                             
                             newMetaSAD <- metaSAD # modifying meta community
                             irare <- which.min(newMetaSAD)[1] # index of a rare species
                             newMetaSAD[irare] <- 2 * max(newMetaSAD) # amplifying rare x2
                             
                             untbPost <- untbSim(J = J, nu = nu, m = m, niter = niter, 
                                                 initLocal = untbPre[nrow(untbPre), ], 
                                                 initMeta = newMetaSAD) # post invasion
                             
                             abundPre <- getAbund(untbPre) # combining sims
                             abundPost <- getAbund(untbPost)
                             abundPost <- abundPost[abundPost$tstep != 0, ] # removing overlap
                             abundPost$tstep <- abundPost$tstep + max(abundPre$tstep) # timesteps post after pre
                             abundAll <- rbind(abundPre, abundPost) # combining
                             
                             oneFinalSAD <- abundPost$abund[abundPost$tstep == 13000]
                             
                             
                             # hill numbers
                             theseHill <- renyi(oneFinalSAD, scales = 0:2, 
                                                hill = TRUE)
                             
                             # last line inside the {} should be the output you want
                             return(theseHill)
                             
                     })
# Summarizing this raw simulation, mean and SD are good options
postinvas_hillMeans <- apply(hillz_5, 1, mean)
postinvas_hillSDs <- apply(hillz_5, 1, sd)




# Long-Term Post
hillz_6 <- replicate(500, # number of times to replicate
                     { # put all the code you need to replicate inside the { }
                             metaSAD <- pika::rlseries(S, 0.001) 
                             untbPre <- untbSim(J = J, nu = nu, m = m, niter = niter,
                                                initMeta = metaSAD) # pre invasion
                             
                             newMetaSAD <- metaSAD # modifying meta community
                             irare <- which.min(newMetaSAD)[1] # index of a rare species
                             newMetaSAD[irare] <- 2 * max(newMetaSAD) # amplifying rare x2
                             
                             untbPost <- untbSim(J = J, nu = nu, m = m, niter = niter, 
                                                 initLocal = untbPre[nrow(untbPre), ], 
                                                 initMeta = newMetaSAD) # post invasion
                             
                             abundPre <- getAbund(untbPre) # combining sims
                             abundPost <- getAbund(untbPost)
                             abundPost <- abundPost[abundPost$tstep != 0, ] # removing overlap
                             abundPost$tstep <- abundPost$tstep + max(abundPre$tstep) # timesteps post after pre
                             abundAll <- rbind(abundPre, abundPost) # combining
                             
                             oneFinalSAD <- abundPost$abund[abundPost$tstep == max(abundPost$tstep)]
                             
                             
                             # hill numbers
                             theseHill <- renyi(oneFinalSAD, scales = 0:2, 
                                                hill = TRUE)
                             
                             # last line inside the {} should be the output you want
                             return(theseHill)
                             
                     })

# Summarizing this raw simulation, mean and SD are good options
longinvas_hillMeans <- apply(hillz_6, 1, mean)
longinvas_hillSDs <- apply(hillz_6, 1, sd)



## Compiling
hills_all <- rbind(neut_hillMeans,neut_hillSDs,
             comp_hillMeans, comp_hillSDs,
             preinvas_hillMeans, preinvas_hillSDs,
             postinvas_hillMeans, postinvas_hillSDs,
             longinvas_hillMeans,longinvas_hillSDs)

hills_all <-round(hills_all, digits = 2)

hills_all <- hills_all %>% 
        rename("0" = "Species Richness",
               "1" = "Shannons Index",
               "2" = "Inverse Simpson",
               "neut_hillMeans" = "Neutral Model",
               "neut_hillSDs" = "SD",
               "comp_hillMeans" = "Competitive Model",
               "comp_hillSDs" = "SD",
               "preinvas_hillMeans" = "Pre-Invasion Model",
               "preinvas_hillSDs" = "SD",
               "postinvas_hillMeans" = "Short-Term Post Invasion",
               "postinvas_hillSDs" = "SD",
               "longinvas_hillMeans" = "Long-Term Post Invasion",
               "longinvas_hillSDs" = "SD")

kable(hills_all, caption = "Mean Biodiversity Indices")
