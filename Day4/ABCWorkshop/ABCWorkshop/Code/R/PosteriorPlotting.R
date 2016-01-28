################################################################################
##                      POSTERIOR CHECKING FOR EARTHWORMS                     ## 
################################################################################

## plot.postCheck() plots a 'posterior predictive check' for an ABC analysis
## of Johnston et al.'s earthworm IBM [1], as done in the 'Parameterisation
## & ABC Workshop'; see the worksheet for full information

## it can do either a 'pseudo check', where it just replots accepted runs, or a
## true 'posterior predictive check', where it re-runs the underlying IBMs;
## to do this, adapt & source 'RunningNetLogo.R' first; see the worksheet

## plot.postCheck() takes the following arguments:
## abcEstObj = an object produced by create.abcEst() in ParameterEstimation.R;
## draws = the number of draws from the approximate
##      posterior parameter distribution to plot
## rerun = whether or not the underlying IBMs should actually be rerun

################################################################################
##                             THE PLOT FUNCTION                              ## 
################################################################################

plot.postCheck <- function(abcEstObj, draws, rerun = FALSE) {	
    ## specify which columns of the data and simulation results relate to
	## Experiment 1 & Experiment 2 respectively
	col.e1 <- c(1:27)
	col.e2 <- c(28:46)

    ## find the empirical data to fit as stored in the abcEst object
    target <- unlist(all.data)
    
    ## randomly select 'draws' of the accepted runs to plot
	ran.draws <- sort(sample(1:nrow(abcEstObj$priors[abcEstObj$accepted, ]),
        draws, replace = TRUE))
    
    if (rerun) { 
		## prepare some objects to store posterior checks in	
		results.best <- data.frame(matrix(nrow = draws, ncol = length(target)))
    		accepted.runs <- data.frame(matrix(nrow = draws, ncol = length(target)))

		## for as many times as specified by 'draws'...
		for (i in 1:draws) {
			## run the model with the 'best' parameter values identified by ABC
			results.best[i, col.e1] <- do.runExp1(abcEstObj$best.values)
			results.best[i, col.e2] <- do.runExp2(abcEstObj$best.values)
			
			## run the model with a random set of accepted parameter values
			accepted.runs[i, col.e1] <- do.runExp1(
				abcEstObj$priors[abcEstObj$accepted, ][ran.draws[i], ])
			accepted.runs[i, col.e2] <- do.runExp2(
				abcEstObj$priors[abcEstObj$accepted, ][ran.draws[i], ])
		}
		
		## take the average of the runs using the 'best' parameter values
		results.best <- colMeans(results.best)
		
    } else {
    		## store the results of the best-fitting accepted run	
    		results.best <- abcEstObj$results[abcEstObj$error == min(abcEstObj$error), ]
    	
    		## store the results of 'draws' of the accepted runs, selected randomly
        accepted.runs <- abcEstObj$results[abcEstObj$accepted, ][ran.draws, ]
    }
    
    ## compute the fit of the best-fitting parameter set		
    rsqs.best <- c()
    rsqs.best[1] <- 1 - (sum((target[col.e1] - results.best[col.e1])^2) /
        sum((target[col.e1] - mean(target[col.e1]))^2))
    rsqs.best[2] <- 1 - (sum((target[col.e2] - results.best[col.e2])^2) /
        sum((target[col.e2] - mean(target[col.e2]))^2))

    ## specify the plot annotations which differ between Exp 1 & Exp 2
    result.cols <- list(col.e1, col.e2)
    x.cors <- list(c(0:26)*7, c(0:18)*10)
    y.lims <- list(0.55, 60)
    y.labs <- list("mass (grams)", "cocoons")
    y.ticks <- list(c(0, 0.25, 0.5), c(0, 20, 40, 60))
    arrow.xcors <- list(0, seq(from = 0, to = 160, by = 20))
    arrow.ycors <- list(c(-0.145, -0.115), c(-15.81818, -12.54545))
    
    ## start the right kind of graphic device for the current operating system
    switch(Sys.info()[1], 
        Windows = { windows(width = 9, height = 3) },
        Linux = { x11(width = 9, height = 3) },
        Darwin = { quartz(width = 9, height = 3) },
        stop("Unknown OS! Can't plot.")
    )
	
    ## set the plot margins and default line width; divide the plot in 3 panels
    par(mai = c(0.7,0.6,0.15,0.1), lwd = 2, mfrow = c(1,3))
	
    ## create the legend in the first panel
    plot(x <- c(0, 1), y <- c(0, 1), type = "n", ann = FALSE, axes = FALSE,
        ylim = c(0,1), xlim = c(0,1))
    legend(0.06, 0.9,
        legend = c("empirical data", "best fitting run", "posterior check"),
        lwd = c(2, 4, 4, 4), pch = c(21, 20, 20), pt.cex = c(1.2, 0.001, 0.001),
        pt.bg = "white", col = c("black", "dimgrey", "grey"), cex = 1.6,
        bty = "n", seg.len = 1, title = expression(bold("legend")))
    arrows(0.115, 0.20, 0.115, 0.25, length = 0.05)
    text(0.44, 0.225, "food added", cex = 1.6)
    
    ## automatically plot the posterior checks for Exp 1 & Exp 2
    for (i in 1:2) {
        ## draw an empty plot
        plot(x <- x.cors[[i]], target[result.cols[[i]]], xlim = c(-1, 183),
            ylim = c(0, y.lims[[i]]), axes = FALSE, ann = FALSE, type = "n")
        
        ## add a border, axis titles and axes ticks to the plot
        box(lwd = 1.5)
        title(xlab = "time (days)", line = 3.85, cex.lab = 1.8)
        title(ylab = y.labs[[i]], line = 2.75, cex.lab = 1.8)
        axis(1, at = seq(0, 180, by = 60), cex.axis = 1.8, mgp = c(3, 1.3, 0))
        axis(2, at = y.ticks[[i]], cex.axis = 1.8, mgp = c(3, 0.8, 0))

        ## add arrows at the times where the earthworms were fed
        for (j in 1:length(arrow.xcors[[i]])) {
            arrows(arrow.xcors[[i]][j], arrow.ycors[[i]][1],
                arrow.xcors[[i]][j], arrow.ycors[[i]][2], xpd = TRUE,
                length = 0.05)
        }
        
        for (j in 1:draws) {
            lines(x.cors[[i]], accepted.runs[j, result.cols[[i]]],
                lwd = 4, col = rgb(0.75, 0.75, 0.75, 0.2), type = "l")
        }
        lines(x <- x.cors[[i]], results.best[result.cols[[i]]],
            cex = 1.5, type = "l", col = "dimgrey", lwd = 4)
        lines(x <- x.cors[[i]], target[result.cols[[i]]], cex = 1.2,
            type = "o", bg = "white", pch = 21)
        legend("topleft", legend = paste("Exp", i), text.font = 2, cex = 1.6,
            bty = "n", inset = c(-0.1, -0.015))
        legend("bottomright", legend = substitute(expression(paste(RSQ, " ",
            italic(R)^2)), list(RSQ = round(rsqs.best[i], 2)))[[2]],
            lty = 1,lwd = 4, col = "dimgrey", bty = "n", cex = 1.5,
            seg.len = 0.4, x.intersp = 0.6, y.intersp = 1.2, inset = c(0.02, 0))
    }
}

