## this script contains the full set of R commands suggested in the worksheet
## for the 'Parameterisation & ABC Workshop'

## if it is loaded using 'source', all functions will be available and the
## suggested plots will be drawn, but no other output will be visible;
## instead, it is suggested to manually copy the instructions into the R console

################################################################################
##             Commands for Step 1: Finding 95% Credible Intervals            ## 
################################################################################

setwd("/Users/Elske/Desktop/ABCWorkshop")

all.data <- read.table("Inputs/all_data.txt", header = TRUE)

full.priors <- read.table("Results/full_priors.txt", header = TRUE,
    stringsAsFactors = FALSE)
full.results <- read.table("Results/full_results.txt", header = TRUE,
    stringsAsFactors = FALSE)

head(full.priors)
head(full.results)

source("Code/R/ParameterEstimation.R")

full.abcEst <- create.abcEst(target = all.data, priors = full.priors,
    results = full.results, rate = 0.001)

summary(full.abcEst)
plot(full.abcEst)

################################################################################
##                    Commands for Step 2: Model Selection                    ##
################################################################################

simple.priors <- read.table("Results/simple_priors.txt", header = TRUE,
    stringsAsFactors = FALSE)
simple.results <- read.table("Results/simple_results.txt", header = TRUE,
    stringsAsFactors = FALSE)

source("Code/R/ModelSelection.R")

all.indexes <- c(full.priors[,1], simple.priors[,1])
all.results <- rbind(full.results, simple.results)

both.abcSel <- create.abcSel(target = all.data, indexes = all.indexes,
    results = all.results, rate = 0.001)
    
col.e1 <- c(1:27)col.e2 <- c(28:46)

e1.abcSel <- create.abcSel(target = all.data[ , col.e1], indexes = all.indexes,
    results = all.results[ , col.e1], rate = 0.001)
e2.abcSel <- create.abcSel(target = all.data[ , col.e2], indexes = all.indexes,
    results = all.results[ , col.e2], rate = 0.001)

summary(both.abcSel)
summary(e1.abcSel)
summary(e2.abcSel)

################################################################################
##                   Commands for Step 3: Posterior Checking                  ##
################################################################################

source("Code/R/PosteriorPlotting.R")

plot.postCheck(full.abcEst, draws = 100, rerun = FALSE)

## after this point in the practical, scripts require adjusting to your local
## setup to run; please see the worksheet

## this code was distributed during the 'Parameterisation & ABC Workshop',
## organised at the University of Bournemouth during a NERC course on
## 'Agent-Based Modelling', 25th - 29th January 2016.

## based on van der Vaart, Beaumont, Johnston & Sibly, 2015,
##      Ecological Modelling, 312, 180 - 190.

## code inspired by the R package 'abc':
## Csillery, Francois & Blum, 2012,
##      Methods in Ecology and Evolution, 3, 475 - 479.

## please use freely, but attribute!
## Elske van der Vaart (elskevdv@gmail.com), University of Reading