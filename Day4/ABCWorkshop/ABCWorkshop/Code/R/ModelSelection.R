################################################################################
##                   ACTUALLY DOING ABC PARAMETER ESTIMATION                  ## 
################################################################################

## create.abcSel() does ABC model selection, with the following arguments:
## target = the empirical data;
## indexes = a list specifying which rows of results belong to which model;
## results = the simulation results;
## rate = the acceptance rate, i.e., the proportion of runs accepted

## to use this function, store its output in an object,
## i.e, at the console, type:
## > my.abcSel <- create.abcSel(target = my.target, indexes = my.indexes,
##      results = my.results, rate = my.accept.rate)
## then look at the result by typing: > summary(my.abcEst)

create.abcSel <- function(target, indexes, results, rate) {
    ## store the original function call, so it can be viewed later
    call <- match.call()
    
    ## make sure the empirical data is formatted as a simple vector
    target <- unlist(target)
	
    ## calculate the standard deviation of each column of simulation results
    sim.sds <-  apply(X = results, MARGIN = 2, FUN = sd)
	
    ## scale the simulation results by their standard deviations
    ## (this to correct for different scales between different types of results)
    scaled.results <- results
    scaled.results[, sim.sds != 0] <- sweep(x = results[, sim.sds != 0],
        MARGIN = 2, STATS = sim.sds[sim.sds != 0], FUN = "/")
	
    ## scale the empirical data in the same way as the simulation results
    scaled.data <- target
    scaled.data[sim.sds != 0] <- target[sim.sds != 0] / sim.sds[sim.sds != 0]
		
    ## for each row, calculate the distance (or error) from each simulation
    ## result to the corresponding empirical data point
    errors <- sweep(x = scaled.results, MARGIN = 2, STATS = scaled.data,
        function(x, y) sqrt((y - x)^2))
	
    ## for each row, sum the distances between each simulation result and the
    ## corresponding empirical data point (giving the total error)
    summed.errors <- apply(X = errors, MARGIN = 1, FUN = sum)
    
    ## calculate the number of runs to accept given the acceptance rate
    number.to.accept <- ceiling(dim(errors)[1] * rate)
	
    ## calculate the max error to accept given the number of runs to accept
    error.to.accept <- sort(summed.errors)[number.to.accept]
	
    ## accept the runs with a total error less than the max acceptable error
    accepted <- (summed.errors <= error.to.accept)
    
    ## for each model, calculate the number of runs that were accepted
    ## unique(indexes) gives the number of different models
    models.accepted <- c()
    for (i in 1:length(unique(indexes))) {
		models.accepted <- c(models.accepted,
            length(indexes[accepted][indexes[accepted] == unique(indexes)[i]]))
    }
    
    ## assemble an object describing what's been done, and return it to the user
    outcome <- list(call = call, target = target, indexes = indexes,
        results = results, rate = rate, errors = summed.errors,
        accepted = accepted, models.accepted = models.accepted)
	class(outcome) <- "abcSel"
    
    return(outcome)
}

################################################################################
##                         PRINT & SUMMARY FUNCTIONS                          ## 
################################################################################

## this function prints a quick overview of an abcSel object
## see this overview by typing > print(name.of.object) or just > name.of.object

print.abcSel <- function(x, ...) {
    cat("call:\n")
    show(x$call)

    cat("\nattributes:\n")
    show(attributes(x))

    cat("# of sumstats:\t\t\t", length(x$target), sep = "")
    cat("\n# of models:\t\t\t\t", length(unique(x$indexes)), sep = "")
    cat("\n# of runs:\t\t\t", nrow(x$results), sep = "")
}

## this function computes the main results
## (i.e., the Bayes factors) of an abcSel object;
## see this summary by typing > summary(name.of.object) in the console 

summary.abcSel <- function(x, ...) {
	accepted <- x$models.accepted
	factors <- matrix(ncol = length(accepted), nrow = length(accepted))
	colnames(factors) <- paste(unique(x$indexes))
	rownames(factors) <- paste(unique(x$indexes))
	for (i in 1:length(accepted)) {
		for (j in 1:length(accepted)) {
			factors[i, j] <- round(accepted[i] / accepted[j], 2)
		}
	}
        
        outcome <- list(call = x$call, num.sumstats = length(x$target),
            num.models = length(unique(x$indexes)),
            num.runs = nrow(x$results), rate = x$rate, factors = factors)
	
        class(outcome) <- "summary.abcSel"
        outcome
}

## this function prints the summary of an abcSel object,
## as produced by the function summary.abcSel()

print.summary.abcSel <- function(x, ...) {
	cat("call:\n")
	show(x$call)

	cat("\acceptance rate:\t\t",  x$rate, sep = "")	
	cat("\n")
	cat("\n# of models:\t\t\t\t", x$num.models, sep = "")
	cat("\n# of sumstats:\t\t\t",  x$num.sumstats, sep = "")
	cat("\n# of runs:\t\t\t\t", x$num.runs, sep = "")
	
	cat("\n\n")
	cat("Bayes factors:\n")
	print(x$factors)
}

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