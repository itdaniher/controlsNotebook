# coding=utf-8

import sympy
from sympy.abc import s
from pylab import *
import numpy
import scipy.signal as signal
j = sympy.I

ct = 1000

decibels = lambda lin: 20*numpy.log10(norm(lin))

def eeval(expression, w):
	""" evaluate a sympy expression at omega. return magnitude, phase."""
	num, den = e2sys(expression)
	y = numpy.polyval(num, 1j*w) / numpy.polyval(den, 1j*w)
	phase = numpy.arctan2(y.imag, y.real) * 180.0 / numpy.pi
	mag = abs(y)
	return mag, phase

def bode(expression, n = 10):
	freqs = signal.findfreqs(e2sys(expression)[0], e2sys(expression)[1], n)
	magnitude = numpy.array([])
	phase = numpy.array([])
	for freq in freqs:
		(m, p) = eeval(expression, freq)
		magnitude = numpy.append(magnitude, m)
		if p >= 0:
			p = -360+p
		phase = numpy.append(phase, p)
	magnitude = array(map(decibels, magnitude))
	return freqs, magnitude, phase

def e2sys(expression):
	""" basic helper function that accepts a sympy expression, expands it, 
	attempts to simplify it, and returns a numerator and denomenator pair for the instantiation of a scipy
	LTI system object. """
	expression = expression.expand()
	expression = expression.cancel()
	n = sympy.Poly(sympy.numer(expression), s).all_coeffs()
	d = sympy.Poly(sympy.denom(expression), s).all_coeffs()
	n = map(float, n)
	d = map(float, d)
	return (n, d)

def sys2e(system):
	""" basic helper function that accepts an instance of the scipy.signal.lti class
	and returns a sympy expression given by the ratio of the expanded numerator and denomenator"""
	den = sum([x*s**y for y,x in enumerate(system.den[::-1])])
	num = sum([x*s**y for y,x in enumerate(system.num[::-1])])
	return sympy.factor(num/den)

def phaseMargin(expression):
	w, mag, phase = bode(expression, n=ct)
	crossingPoint = argmin(abs(mag)) # point where magnitude is closest to 0dB
	return {"w_c": w[crossingPoint], "p_m": phase[crossingPoint]+180}

def steadyStepError(expression):
	""" use sympy's symbolic solver at 0 to determine steady state error to a step. """
	e = expression.subs(s, j*0)
	fvt = float(sympy.Abs(e))
	return {"sse": 1/(1+fvt)}

def reducedGainCompensate(expression, target):
	""" find the uncompensated system's magnitude at the right phase to give the right margin, normalize against it"""
	w, mag, phase = bode(expression, n=ct)
	phaseTarget = (-180 + target)
	targetIndex = argmin(abs(phaseTarget - phase))
	magnitude = mag[targetIndex]
	magnitude = 10**(magnitude/20)
	return {"K": 1/magnitude}, 1/magnitude

def dominatePoleCompensate(expression, target):
	""" add a pole, call reducedGainCompensate"""
	K = reducedGainCompensate(1/s*expression, target)[0]['K']
	return {"K": K}, K/s

def lagCompensate(expression, target):
	w, mag, phase = bode(expression, n=ct)
	phaseTarget = -180 + (target + 6) # find location of new crossing goal; extra '-6' is from 10/wc rule
	targetIndex = argmin(abs(phaseTarget - phase))
	w_c = w[targetIndex]
	tau = 10/w_c
	# equation for basic lag compensator
	alpha = eeval(expression, w_c)[0]
	G_c = (tau*s+1)/(alpha*tau*s+1)
	return {"alpha": alpha, "tau": tau}, G_c

def leadCompensate(expression, target):
	# alpha = 10 for 55 degrees phase margin
	w, mag, phase = bode(expression, n=ct)
	alpha = 10
	phaseTarget = -(180 + (55 - target))
	targetIndex = argmin(abs(phase-phaseTarget))
	w_c = w[targetIndex]
	tau = 1 / (sqrt(alpha) * w_c)
	G_c = (alpha*tau*s+1)/(tau*s+1)
	# get magnitude of loop transfer function L(s) at w_c
	K_l = 1/ eeval(G_c*expression, w_c)[0] 
	# equation for basic lead compensator
	G_c = K_l * (alpha*tau*s+1)/(tau*s+1)
	return {"K_l": K_l, "alpha": alpha, "tau":tau}, G_c

def drawBode(expression, f=figure(), color="k", labeled=""):
	""" draw bode plot for a given scipy.signal.lti instance. """
	adjustprops = dict(left=0.1, bottom=0.1, right=0.97, top=0.93, hspace=0.2)
	f.subplots_adjust(**adjustprops)
	magPlot = f.add_subplot(2, 1, 1)
	phasePlot = f.add_subplot(2, 1, 2, sharex=magPlot)
	w, mag, phase = bode(expression)
	phase = [i if i < 0 else -(360-i) for i in phase]
	magPlot.semilogx(w, mag, label=labeled, color=color)
	magPlot.set_title("magnitude")
	magPlot.set_ylabel("amplitude (dB)")
	phasePlot.semilogx(w, phase, label=labeled, color=color)
	phasePlot.set_title("phase")
	phasePlot.set_xlabel("radians per second")
	phasePlot.set_ylabel("degrees")
	setp(magPlot.get_xticklabels(), visible=False)
	legend(loc="best")

if __name__ == "__main__":
	print "before :"
	#G_p = 100/((s+1)*(0.1*s+1)*(0.01*s+1))
	G_p = 10/(s*(s+1)*(0.1*s+1))
	print phaseMargin(G_p)
	print steadyStepError(G_p)
	drawBode(G_p, color="k", labeled="uncompensated")
	print "leadCompensate"
	(coeffs, G_c) = leadCompensate(G_p, 50)
	print coeffs
	print phaseMargin(G_p*G_c)
	print steadyStepError(G_p*G_c)
	drawBode(G_p*G_c, color="b", labeled="lead")
	print "lagCompensate"
	(coeffs, G_c) = lagCompensate(G_p, 50)
	print coeffs
	print phaseMargin(G_p*G_c)
	print steadyStepError(G_p*G_c)
	drawBode(G_p*G_c, color="r", labeled="lag")
	print "reducedGain"
	(coeffs, G_c) = dominatePoleCompensate(G_p, 50)
	print coeffs
	print phaseMargin(G_p*G_c)
	print steadyStepError(G_p*G_c)
	drawBode(G_p*G_c, color="g", labeled="dominatePole")
	show()
