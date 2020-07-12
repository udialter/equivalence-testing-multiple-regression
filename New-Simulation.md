New sim
================
Udi Alter
20/06/2020

``` r
# universal/general pre-settings
alpha <- .05 # Type I error rate
nsim <- 5000 # number of replications
delta <-.15 # SESOI

# preparing empty place-holders for later use in the function
betas <- vector(mode = "double", length = 5) # will hold beta weight values
berror <- vector(mode = "double", length = 5) # will hold standard error for each beta 

# preparing a p-value count table
pval <-matrix(0, nrow = nsim, ncol =20) 
colnames(pval)<-c('DB.pval.b1', 'DB.pval.b2', 'DB.pval.b3', 'DB.pval.b4', 'DB.pval.b5', 'TOST.pval.b1a', 'TOST.pval.b1b', 'TOST.pval.b2a', 'TOST.pval.b2b', 'TOST.pval.b3a','TOST.pval.b3b', 'TOST.pval.b4a', 'TOST.pval.b4b', 'TOST.pval.b5a', 'TOST.pval.b5b', 'AH.pval.b1', 'AH.pval.b2', 'AH.pval.b3', 'AH.pval.b4', 'AH.pval.b5')

# preparing the results intermediate row
results <- matrix(data = 0, nrow =1, ncol = 15) #number of tests concluding lack of association
colnames(results)<-c('DB.b1', 'DB.b2', 'DB.b3', 'DB.b4', 'DB.b5', 'TOST.b1', 'TOST.b2', 'TOST.b3', 'TOST.b4', 'TOST.b5', 'AH.b1', 'AH.b2', 'AH.b3', 'AH.b4', 'AH.b5')

# preparing the final results matrix
final <- matrix(data = 0, nrow =6, ncol = 15) #number of tests concluding lack of association
colnames(final)<-c('DB.b1', 'DB.b2', 'DB.b3', 'DB.b4', 'DB.b5', 'TOST.b1', 'TOST.b2', 'TOST.b3', 'TOST.b4', 'TOST.b5', 'AH.b1', 'AH.b2', 'AH.b3', 'AH.b4', 'AH.b5')
rownames(final) <- c("n = 50", "n = 75", "n = 100", "n = 250", "n = 500", "n = 1000")
```

``` r
simulate_me <- function(N)  { #begin function

  for(i in 1:nsim) {  # for each of the 5000 replications, perform:
                    sim_arguments <- list( # Setting fixed parameters to remain constant for every simulated dataset replication
                                          formula = y ~ 1 + x1 + x2 + x3 + x4 + x5,
                                          fixed = list(x1 = list(var_type = 'continuous', mean = 0, sd = 1),
                                                       x2 = list(var_type = 'continuous', mean = 0, sd = 1),
                                                       x3 = list(var_type = 'continuous', mean = 0, sd = 1),
                                                       x4 = list(var_type = 'continuous', mean = 0, sd = 1),
                                                       x5 = list(var_type = 'continuous', mean = 0, sd = 1)),
                                          error = list(variance = 1),
                                          sample_size = N,
                                          reg_weights = c(0, 0, 0.05, 0.1, 0.15, 0.2)
                                          )
                    # The actual simulation of data with the above fixed parameters
                    model <- simulate_fixed(data = NULL, sim_arguments) %>%
                    simulate_error(sim_arguments) %>%
                    generate_response(sim_arguments) %>% 
                    model_fit(sim_arguments) %>%
                    extract_coefficients()

                    df <- N-(5+1) # degrees of freedom, sample size minus the number of predictor plus 1, all condition will have 5 predictors.

                    for(b in 1:5) { # On each of the 5 predictors within each replication,conduct 3 tests: difference-based, two one-sided tests, and Anderson-Hauck
                                   betas[b] <- model$estimate[b+1] #Beta weights estimates extraction
                                   berror[b] <- model$std.error[b+1] #Standard error extraction per predictor
           
                                   ## traditional difference-based test 
                                   pval[i,b] <- model$p.value[b+1] # p value extraction per b
                                   ifelse(pval[i,b] >= alpha, results[1,b] <-results[1,b]+1, results[1,b] <- results[1,b]) # collecting non-sgnificant p values per b
            
                                   ## Schuirmann equivalence test TOST
                                   t1<- (betas[b]-(-delta))/(sqrt(berror[b]^2)) 
                                   t2<- (delta-betas[b])/(sqrt(berror[b]^2))
                                   p1<- 1-pt(t1, df) 
                                   p2<- 1-pt(t2, df)
                                   pval[i,2*b+4] <- p1 
                                   pval[i,2*b+5] <- p2
                                   ifelse(p1 < alpha & p2 < alpha, results[1,b+5] <-results[1,b+5]+1, results[1,b+5] <- results[1,b+5])
            
                                   ## Anderson Hauck procedure
                                   tah <- model$statistic[b+1]
                                   pah <- pt((abs(betas[b])-delta)/(sqrt(berror[b]^2)),df) - pt((-abs(betas[b])-delta)/(sqrt(berror[b]^2) ), df)
                                   pval[i,b+15] <- pah
                                   ifelse(pah < alpha, results[1,b+10] <-results[1,b+10]+1, results[1,b+10] <-   results[1,b+10])
            
                                  } # end of inner loop (1:5)
                    
                    } # end of outer loop (1:nsim)

  return(results/nsim)

} # end of function
```

# looping the function across all sample size conditions

``` r
sampsizes <- c(50, 75, 100, 250, 500, 1000)

for(n in seq_along(sampsizes)) {
                                final[n,] <- simulate_me(sampsizes[n])
                                }

final
```

    ##           DB.b1  DB.b2  DB.b3  DB.b4  DB.b5 TOST.b1 TOST.b2 TOST.b3 TOST.b4
    ## n = 50   0.9482 0.9364 0.9040 0.8400 0.7316  0.0000  0.0000  0.0000  0.0000
    ## n = 75   0.9434 0.9340 0.8780 0.7578 0.6266  0.0006  0.0004  0.0002  0.0000
    ## n = 100  0.9538 0.9266 0.8328 0.6930 0.5280  0.0058  0.0054  0.0036  0.0014
    ## n = 250  0.9512 0.8856 0.6586 0.3680 0.1208  0.5188  0.3974  0.1840  0.0562
    ## n = 500  0.9532 0.7944 0.3882 0.0822 0.0078  0.9080  0.7154  0.2890  0.0458
    ## n = 1000 0.9514 0.6510 0.1036 0.0042 0.0000  0.9988  0.9380  0.4734  0.0514
    ##          TOST.b5  AH.b1  AH.b2  AH.b3  AH.b4  AH.b5
    ## n = 50    0.0000 0.0860 0.0792 0.0622 0.0540 0.0306
    ## n = 75    0.0000 0.1100 0.1068 0.0772 0.0544 0.0272
    ## n = 100   0.0004 0.1514 0.1394 0.0900 0.0476 0.0230
    ## n = 250   0.0058 0.5288 0.4066 0.1896 0.0574 0.0060
    ## n = 500   0.0032 0.9080 0.7154 0.2890 0.0458 0.0032
    ## n = 1000  0.0004 0.9988 0.9380 0.4734 0.0514 0.0004
