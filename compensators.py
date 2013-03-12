# coding=utf-8

import sympy
from sympy.abc import s#, j
from pylab import *
from scipy.signal import *

j = sympy.I

blacks = lambda L: L/(1+L)

ct = 10000

def eeval(expression, w):
	""" evaluate a sympy expression at omega. return magnitude, phase."""
	e = expression.subs(s, j*w)
	return float(sympy.Abs(e)), float(sympy.arg(e)/sympy.pi*180)

def e2sys(expression):
	""" basic helper function that accepts a sympy expression, expands it, 
	attempts to simplify it, and returns a scipy-type LTI system equivalent. """
	expression = expression.expand()
	expression = expression.cancel()
	n = sympy.Poly(sympy.numer(expression), s).all_coeffs()
	d = sympy.Poly(sympy.denom(expression), s).all_coeffs()
	n = map(float, n)
	d = map(float, d)
	return lti(n, d)

def sys2e(system):
	""" basic helper function that accepts an instance of the scipy.signal.lti class
	and returns a sympy expression given by the ratio of the expanded numerator and denomenator"""
	den = sum([x*s**y for y,x in enumerate(system.den[::-1])])
	num = sum([x*s**y for y,x in enumerate(system.num[::-1])])
	return sympy.factor(num/den)

def phaseMargin(system):
	w, mag, phase = bode(system, n=ct)
	crossingPoint = argmin(abs(mag)) # point where magnitude is closest to 0dB
	return {"w_c": w[crossingPoint], "p_m": 180+phase[crossingPoint]}

def lagCompensate(system, target):
	w, mag, phase = bode(system, n=ct)
	phaseTarget = (-180 + (target - 6)) # find location of new crossing goal; extra '-6' is from 10/wc rule
	print phaseTarget
	phaseTargetLocation = argmin(abs(phaseTarget - phase)) 
	tau = 10/w[phaseTargetLocation]
	for alpha in logspace(0, 2, ct): 
		# equation for basic lag compensator
		G_c = (tau*s+1)/(alpha*tau*s+1)
		if abs(phaseMargin(e2sys(G_p*G_c))[1]-target) < 0.1:
			break
	return {"alpha": alpha, "tau": tau}

def leadCompensate(system, target):
	# alpha = 10 for 55 degrees phase margin
	w, mag, phase = bode(system, n=ct)
	phase = array([i if i < 0 else -(360-i) for i in phase])
	alpha = 10
	phaseTarget = -(180 + (55 - target))
	targetIndex = argmin(abs(phase-phaseTarget))
	w_c = w[targetIndex]
	tau = 1 / (sqrt(alpha) * w_c)
	K_l = 1
	G_c = K_l * (alpha*tau*s+1)/(tau*s+1)
	# get magnitude of loop transfer function L(s) at w_c
	K_l = 1/ eeval(G_c*sys2e(system), w_c)[0] 
	# equation for basic lead compensator
	G_c = K_l * (alpha*tau*s+1)/(tau*s+1)
	return {"K_l": K_l, "alpha": alpha, "tau":tau}, G_c

def characterize(system):
	""" draw impulse, step, bode for a given scipy.signal.lti instance. """
	adjustprops = dict(left=0.1, bottom=0.1, right=0.97, top=0.93, hspace=0.5)
	f = figure()
	f.subplots_adjust(**adjustprops)
	impulsePlot = f.add_subplot(4, 1, 1)
	stepPlot = f.add_subplot(4, 1, 2, sharex=impulsePlot)
	magPlot = f.add_subplot(4, 1, 3)
	phasePlot = f.add_subplot(4, 1, 4, sharex=magPlot)
	w, mag, phase = bode(system)
	phase = [i if i < 0 else -(360-i) for i in phase]
	tImpulse, youtImpulse = impulse2(system)
	tStep, youtStep = step2(system)
	impulsePlot.plot(tImpulse,youtImpulse)
	impulsePlot.set_title("impulse")
	impulsePlot.set_ylabel("amplitude")
	stepPlot.plot(tStep, youtStep)
	stepPlot.set_title("step")
	stepPlot.set_xlabel("seconds")
	stepPlot.set_ylabel("amplitude")
	magPlot.semilogx(w, mag)
	magPlot.set_title("magnitude")
	magPlot.set_ylabel("amplitude")
	phasePlot.semilogx(w, phase)
	phasePlot.set_title("phase")
	phasePlot.set_xlabel("radians per second")
	phasePlot.set_ylabel("degrees")
	setp(magPlot.get_xticklabels(), visible=False)
	setp(impulsePlot.get_xticklabels(), visible=False)
	show()

G_p = 10/(s*(s+1)*(0.1*s+1))

system = e2sys(G_p)

print phaseMargin(system)
# should be tau=16.9, alpha=14 for lag compensation
#print lagCompensate(system,50)
# should be tau =  0.086, K_l = 0.47, alpha = 10
#print leadCompensate(system, 50)

G_c = leadCompensate(system, 50)[-1]

print phaseMargin(e2sys(G_p*G_c))

characterize(e2sys(G_p*G_c))
