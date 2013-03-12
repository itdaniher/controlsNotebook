# coding=utf-8

from scipy.signal import *
from pylab import *

# It uses scipy.integrate.odeint, which in turn uses the Fortran library
# LSODA.  Like any good solver, LSODA is an adaptive solver--it adjusts its
# step size to be as large as possible while keeping estimates of the error
# bounded.  For the problem you are solving, with initial condition 0, the
# exact solution is initially exactly 0.  This is such a nice smooth solution
# that the solver's step size quickly grows--so big, in fact, that it skips
# right over your pulse and never sees it.

taus = [0, 0.01, 0.1, 0.5, 1.0]

tauN = 0

adjustprops = dict(left=0.1, bottom=0.1, right=0.97, top=0.93, hspace=0.5)   
f = figure(figsize=(8,8), dpi=300)
f.subplots_adjust(**adjustprops) 
impulsePlot = f.add_subplot(4, 1, 1)
title(unicode("τ  num = "+str(taus[0])+"..."+str(taus[-1])+"; τ denom = "+str(tauN)))
stepPlot = f.add_subplot(4, 1, 2, sharex=impulsePlot)
magPlot = f.add_subplot(4, 1, 3)
phasePlot = f.add_subplot(4, 1, 4, sharex=magPlot)

for tau in taus:

	tauD = tau

	system = lti([1], polymul([tauD, 1], [1, 1, 1]))
	
	w, mag, phase = system.bode()
	tImpulse, youtImpulse = system.impulse()
	tStep, youtStep = system.step()
	
	impulsePlot.plot(tImpulse,youtImpulse)#, label="impulse")
	ylabel("amplitude")

	stepPlot.plot(tStep, youtStep)#, label="step")
	xlabel("seconds")
	ylabel("amplitude")
	
	#title("zeroes of "+str(system.zeros)+"; poles of "+str(system.poles), size="small")
	magPlot.semilogx(w, mag)#, label="magnitude")
	ylabel("amplitude")
	
	phasePlot.semilogx(w, phase)#, label="phase")
	xlabel("radians per second")
	ylabel("degrees")
	

setp(magPlot.get_xticklabels(), visible=False)
setp(impulsePlot.get_xticklabels(), visible=False)
f.savefig("figure.png")
