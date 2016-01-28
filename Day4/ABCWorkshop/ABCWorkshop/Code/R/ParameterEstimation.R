################################################################################
##                   ACTUALLY DOING ABC PARAMETER ESTIMATION                  ## 
################################################################################

## create.abcEst() does ABC parameter estimation, with the following arguments:
## target = the empirical data;
## priors = the parameters used to generate the simulation results;
## where the first column is a model index;
## results = the simulation results;
## rate = the acceptance rate, i.e., the proportion of runs accepted

## to use this function, store its output in an object,
## i.e, at the console, type:
## > my.abcEst <- create.abcEst(target = my.target, priors = my.priors,
##      results = my.results, rate = my.accept.rate)
## then look at the result by typing: > summary(my.abcEst)

create.abcEst <- function(target, priors, results, rate) {
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
    
    ## find the run with the smallest total error
    best <- priors[summed.errors == min(summed.errors),]
	
    ## calculate the median values of the accepted parameter values
    medians <- apply(X = priors[accepted, 2:ncol(priors)], MARGIN = 2,
    	FUN = median)
    med.values <- best
    med.values[2:ncol(priors)] <- medians
	
    ## assemble an object describing what's been done, and return it to the user
    outcome <- list(call = call, target = target, priors = priors,
        results = results, rate = rate, errors = summed.errors,
        accepted = accepted, best.values = best, med.values = med.values)
	class(outcome) <- "abcEst"
    
    return(outcome)
}

################################################################################
##                       PLOT, PRINT & SUMMARY FUNCTIONS                      ## 
################################################################################

## this function prints a quick overview of an abcEst object
## see this overview by typing > print(name.of.object) or just > name.of.object

print.abcEst <- function(x, ...) {
    cat("call:\n")
    show(x$call)

    cat("\nattributes:\n")
    show(attributes(x))

    cat("type of model:\t\t\t'", x$priors[1,1], "'", sep = "")
    cat("\n")
    cat("\n# of sumstats:\t\t\t", length(x$target), sep = "")
    cat("\n# of parameters:\t\t", length(x$priors[1,]) - 1, sep = "")
    cat("\n# of runs:\t\t\t", nrow(x$results), sep = "")
}

## this function computes the main results
## (i.e., the posterior parameter distributions) of an abEst object;
## see this summary by typing > summary(name.of.object) in the console 

summary.abcEst <- function(x, cred.int = 0.95, ...) {
    acc.params <- x$priors[x$accepted, 2:ncol(x$priors)]
    mins <- apply(acc.params, 2, min)
    lows <- apply(acc.params, 2, quantile, ((1 - cred.int) / 2))
    medians <- apply(acc.params, 2, quantile, 0.5)
    means <- apply(acc.params, 2, mean)
    highs <- apply(acc.params, 2, quantile, (1 - ((1 - cred.int) / 2)))
    maxes <- apply(acc.params, 2, max)
    sums <- rbind(mins, lows, medians, means, highs, maxes)
    rownames(sums) <- c("Min.:",
        paste((((1 - cred.int) / 2) * 100), "% Perc.:", sep = ""),
        "Median:", "Mean:",
        paste(((1 - (1 - cred.int) / 2) * 100), "% Perc.:", sep = ""), "Max.:")
	
    outcome <- list(call = x$call, model = x$priors[1,1],
        num.sumstats = length(x$target), num.params = ncol(x$priors) - 1,
        num.runs = nrow(x$results), rate = x$rate, sums = sums)        
	
    class(outcome) <- "summary.abcEst"
    outcome
}

## this function prints the summary of an abcEst object,
## as produced by the function summary.abcEst()

print.summary.abcEst <- function(x, digits = 3, ...) {
    cat("call:\n")
    show(x$call)

    cat("\ntype of model:\t\t'", x$model, "'", sep = "")
    cat("\nacceptance rate:\t\t",  x$rate, sep = "")	
    cat("\n")
    cat("\n# of sumstats:\t\t", x$num.sumstats, sep = "")
    cat("\n# of parameters:\t\t", x$num.params, sep = "")
    cat("\n# of runs:\t\t\t", x$num.runs, "\n\n", sep = "")
	
    cat("posteriors:\n")
    print(signif(x$sums, digits))
}

## this function plots the main results
## (i.e., the posterior parameter distributions) of an abcEst object
## see this plot by typing > plot(name.of.object) in the console 

plot.abcEst <- function(x, ...) {
    ## find the parameter values stored in the parameter estimation object
    ## ignore the first column, as it's model indexes 
    params <- x$priors[, c(2:ncol(x$priors))]
    num.params <- ncol(params)
	
    ## select the parameter values that were accepted
    posteriors <- params[x$accepted, ]
	
    ## calculate the mean of each parameter's prior values 
    prior.means <- colMeans(params)

    ## scale each parameter's priors and posteriors by the means of their priors
    scaled.priors <- sweep(x = params, MARGIN = 2, STATS = prior.means, "/")
    scaled.posts <- sweep(x = posteriors, MARGIN = 2, STATS = prior.means, "/")

    ## format the parameter names for attractive plotting;
    ## letters after "_"'s are subscripted
    par.names <- c()
	
    for (i in 1:num.params) {
	if (grepl("_", colnames(params)[i])) {
	    parts <- strsplit(colnames(params)[i], "_")
	    fixed <- paste("italic(", parts[[1]][1], "[", parts[[1]][2], "])",
                sep = "")
        } else {
	    fixed <- paste("italic(\"", colnames(params)[i], "\"[])", sep = "") 
	}
	par.names <- c(par.names, parse(text = fixed))
    }

    ## calculate the 95% credible intervals of the priors
    low.priors <- apply(scaled.priors, 2, quantile, c(0.025))
    top.priors <- apply(scaled.priors, 2, quantile, c(0.975))
	
    ## calculate the 95% credible intervals of the posteriors
    low.posts <- apply(scaled.posts, 2, quantile, c(0.025))
    top.posts <- apply(scaled.posts, 2, quantile, c(0.975))

    ## specify the x-coordinates for plotting
    x.cors <- c(0.75, seq(1.75, num.params, by = 1))

    ## start the right kind of graphic device for the current operating system
    graph.args <- list(width = 0.75 + ((num.params + 1) * 0.35), height = 2)
    switch(Sys.info()[1], 
        Windows = { do.call(windows, graph.args) },
        Linux = { do.call(x11, graph.args) },
        Darwin = { do.call(quartz, graph.args) },
        stop("Unknown operating system! Can't plot.")
    )
	
    ## set the plot margins and default line width
    par(mai=c(0.6, 0.6, 0.15, 0.15), lwd = 2)

    ## draw an empty plot
    plot(x.cors, colMeans(scaled.priors), type = "n", axes = FALSE,
        xlim = c(0.125 * (0.75 + num.params * 0.35),
        (num.params + 0.75 - (0.125 * (0.75 + num.params * 0.35)))),
        col = "grey90", ylim = c(min(scaled.priors), max(scaled.priors)))

    ## add the priors to the plot
    segments(x.cors, low.priors, x.cors, top.priors, col = "grey60", lty = 2)
    segments(x.cors - 0.25, low.priors, x.cors + 0.25, low.priors,
        col = "grey60", lty = 2)
    segments(x.cors - 0.25, top.priors, x.cors + 0.25, top.priors,
        col = "grey60", lty = 2)
    points(x.cors, colMeans(scaled.priors), pch = 21, bg = "grey60",
        col = "grey60")

    ## add the posteriors to the plot
    segments(x.cors, low.posts, x.cors, top.posts)
    segments(x.cors - 0.25, low.posts, x.cors + 0.25, low.posts)
    segments(x.cors - 0.25, top.posts, x.cors + 0.25, top.posts)
    points(x.cors, colMeans(scaled.posts), pch = 21, bg = "white")
	
    ## add a border, parameter names, axis titles and y-axis ticks to the plot
    box(lwd = 1.5)
    for (i in 1:num.params) {
	axis(1, at = x.cors[i], labels = par.names[i], cex.axis = 0.9,
            lwd.ticks = 2, mgp = c(3, 0.6, 0), las = 1, tck = -0.05)
    }
    title(xlab = "parameter", line = 1.5, cex.lab = 0.9)
    title(ylab = "scaled value", line = 1.75, cex.lab = 0.9)
    axis(2, cex.axis = 0.9, lwd.ticks = 2, mgp = c(3, 0.5, 0), las = 1,
        tck = -0.05)
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