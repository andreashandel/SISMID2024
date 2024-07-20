library('DSAIRM') # within-host models
library('DSAIDE') # between-host models
library('pracma') # for trapezoid integration function

###########################
# within-host model
###########################

# vector of treatment times
txvec = c(seq(0.1,2,length=10),30)

# vector to hold total virus load
Vtot = rep(0,length(txvec))

# run within-host model for different timing of treatment
for (i in 1:length(txvec))
{
	sim1 <-	simulate_virusandtx_ode(
		U = 10000000,	I = 0, V = 1,
		n = 0, dU = 0,
		dI = 1, dV = 2,
		b = 0.03, p = 0.001, g = 0, e = 0.95,
		f = 0,
		tstart = 0, tfinal = 30, dt = 0.1,
		steadystate = FALSE,
		txstart = txvec[i]
	)
	# plot virus load, just as diagnostic during coding process
	plot(sim1$ts$time,sim1$ts$V,type='l')
	
	# total virus load
	Vtot[i] = log(pracma::trapz(sim1$ts$time,sim1$ts$V))	
}

# plot virus load as function of treatment start
# just as a check to make sure things look ok
plot(txvec,Vtot,type='p')

###########################
# population level model
###########################

# different values for transmission rate b
b0 = 0.002
Vmax = max(Vtot)
bvec = b0 * head(Vtot,-1)/Vmax

# vector that will hold total number of infected
Itot = rep(0,length(bvec))

# run between-host model for different transmission values (virus load)
for (i in 1:length(bvec))
{
	sim2 <-	simulate_SIR_model_ode(
			S = 1000, I = 1, R = 0,
			b = bvec[i],
			g = 1,
			tstart = 0,
			tfinal = 100,
			dt = 0.1
		)
	Itot[i] = head(sim2$ts$S,1) - tail(sim2$ts$S,1)
}

# plot results
plot(head(txvec,-1), Itot, type = 'p', xlab = 'treatment start', ylab = 'outbreak size')

