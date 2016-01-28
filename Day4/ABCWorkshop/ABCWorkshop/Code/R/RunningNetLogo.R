################################################################################
##                        RE-RUNNING THE NETLOGO MODELS                       ## 
################################################################################

## this script makes it possible to run Johnston et al.'s earthworm IBM [1]
## directly from R, in both its "Simple" and "Full" versions, as done in the
## 'Parameterisation & ABC Workshop' - after adapting the paths in the first
## section; see the worksheet for further information

################################################################################
##                         LOADING THE MODELS INTO R                          ## 
################################################################################

Sys.setenv(NOAWT = 1)
library(RNetLogo)

## setting the path to the NetLogo application
nl.path <- "/Applications/NetLogo 5.2"

## loading the NetLogo models for Experiment 1
exp1_full.netlogo <- "exp1_full"
NLStart(nl.path, gui = FALSE, nl.obj = exp1_full.netlogo)
NLLoadModel("/Users/Elske/Desktop/ABCWorkshop/Code/NetLogo/Exp1_Full.nlogo",
	nl.obj = exp1_full.netlogo)

exp1_simple.netlogo <- "exp1_simple"
NLStart(nl.path, gui = FALSE, nl.obj = exp1_simple.netlogo)
NLLoadModel("/Users/Elske/Desktop/ABCWorkshop/Code/NetLogo/Exp1_Simple.nlogo",
	nl.obj = exp1_simple.netlogo)
		
## loading the NetLogo models for Experiment 2
exp2_full.netlogo <- "exp2_full"
NLStart(nl.path, gui = FALSE, nl.obj = exp2_full.netlogo)
NLLoadModel("/Users/Elske/Desktop/ABCWorkshop/Code/NetLogo/Exp2_Full.nlogo",
	nl.obj = exp2_full.netlogo)

exp2_simple.netlogo <- "exp2_simple"
NLStart(nl.path, gui = FALSE, nl.obj = exp2_simple.netlogo)
NLLoadModel("/Users/Elske/Desktop/ABCWorkshop/Code/NetLogo/Exp2_Simple.nlogo",
	nl.obj = exp2_simple.netlogo)

################################################################################
##                             RUNNING THE MODELS                             ## 
################################################################################

## this function returns Johnston et al.'s (2014) original literature values
## version = "Full" or "Simple"

make.litValues <- function(version) {
	values <- data.frame(
		Model = version,
		B_0 = 967,
		E = 0.25,
		E_c = 7,
		E_f = 10.6,
		E_s = 3.6,
		IGm = 0.7,
		h = 3.5,
		M_b = 0.011,
		M_c = 0.015,
		M_p = 0.25,
		M_m = 0.5,
		r_B = 0.177,
		r_m = 0.182,
		s = 0.0004
	)
	
	if (version == "Simple") {
		values$IGm[1] <- 0.15
		values <- values[, !colnames(values) %in% c("E_f", "h", "s")]
	}
	
	values
}

## these functions use RNetLogo to run the earthworm IBMs in NetLogo
## settings = a row of parameter values, where the first column is a
## model index, either "Full" or "Simple"

do.runExp1 <- function(settings) {
	instance <- make.instance("Exp1", settings[1])
	set.parameters(settings, instance)
	NLCommand("setup-interface", nl.obj = instance)
	NLCommand("setup", nl.obj = instance)
	birth <- NLReport("mean-mass", nl.obj = instance)
	mass <- unlist(NLDoReport(iterations = 26, command = "go-7",
		reporter = "mean-mass", nl.obj = instance))
	result <- c(birth, mass)
	result
}

do.runExp2 <- function(settings) {
	instance <- make.instance("Exp2", settings[1])
	set.parameters(settings, instance)
	NLCommand("setup-interface", nl.obj = instance)
	NLCommand("setup", nl.obj = instance)
	NLCommand("go-25", nl.obj = instance)
	start <- c(NLReport("sum-hatchlings", nl.obj = instance))
	outcome <- unlist(NLDoReport(iterations=18, command="go-10",
		reporter = c("sum-hatchlings"), nl.obj = instance))
	result <- c(start, outcome)
	result
}

## this function returns the right NetLogo instance to run,
## given the experiment & version required
## experiment = "Exp1" or "Exp2", version = "Full" or "Simple"

make.instance <- function(experiment, version) {
	instance <- NULL
	if (experiment == "Exp1") {
		if (version == "Full") { instance = exp1_full.netlogo }
		if (version == "Simple") { instance = exp1_simple.netlogo }
	} else if (experiment == "Exp2") {
		if (version == "Full") { instance = exp2_full.netlogo }
		if (version == "Simple") { instance = exp2_simple.netlogo }
	}
	instance
}

## this function loads a set of parameter values into a NetLogo instance
## settings = a row of parameter values; instance = the NetLogo instance
## to load into, produced by the function 'make.instance'

set.parameters <- function(settings, instance) {
	NLCommand(paste("set B_0",  settings$B_0), nl.obj = instance)
	NLCommand(paste("set activation_energy", settings$E), nl.obj = instance)
	NLCommand(paste("set energy_tissue", settings$E_c), nl.obj = instance)
	if ("E_f" %in% colnames(settings)) {
		NLCommand(paste("set energy_food", settings$E_f), nl.obj = instance)
	}
	NLCommand(paste("set energy_synthesis", settings$E_s), nl.obj = instance)	
	NLCommand(paste("set max_ingestion_rate", settings$IGm), nl.obj = instance)
	if ("h" %in% colnames(settings)) {
		NLCommand(paste("set half_saturation_coef", settings$h), nl.obj = instance)
	}
	NLCommand(paste("set mass_birth", settings$M_b), nl.obj = instance)
	NLCommand(paste("set mass_cocoon", settings$M_c), nl.obj =  instance)
	NLCommand(paste("set mass_sexual_maturity", settings$M_p), nl.obj = instance)
	NLCommand(paste("set mass_maximum", settings$M_m), nl.obj = instance)
	NLCommand(paste("set growth_constant", settings$r_B), nl.obj = instance)
	NLCommand(paste("set max_reproduction_rate", settings$r_m), nl.obj = instance)
	if ("s" %in% colnames(settings)) {
		NLCommand(paste("set speed", settings$s), nl.obj = instance)
	}
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