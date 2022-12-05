S <- 50
m <- matrix(runif(S^2, 0, 0.01), nrow = S)
diag(m) <- runif(S, 0.01, 0.2)

n <- rep(10, S)

d0 <- ((m %*% n) * n)[, 1]
la <- 1 * max(d0)

ntstep <- 12

nmat <- matrix(0, nrow = ntstep, ncol = S)

nmat[1, ] <- n


for(i in 2:ntstep) {
    # birth or death?
    # death probabilities
    d <- ((m %*% n) * n)[, 1]
    
    # birth probabilities
    b <- la * n
    
    # `event = 0` is death and `event = 1` is birth
    event <- sample(0:1, size = 1, prob = c(sum(d), sum(b)))
    
    if(event == 0) { # someone dies
        # who dies?
        thisOneDead <- sample(1:S, size = 1, prob = d)
        
        # update populations
        n[thisOneDead] <- n[thisOneDead] - 1
    } else { # someone is born
        thisOneBorn <- sample(1:S, size = 1, prob = b)
        
        # update populations
        n[thisOneBorn] <- n[thisOneBorn] + 1
    }
    
    nmat[i, ] <- n
}

# timeseries
matplot(nmat, type = 'l', lty = 1, col = viridis::viridis(S))

# sad
x <- nmat[ntstep, ]
x <- sort(x[x > 0], TRUE)
plot(x, log = 'y')