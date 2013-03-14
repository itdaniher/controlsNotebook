# coding=utf-8

import sympy
from sympy.abc import s
from pylab import *
from scipy.signal import *

j = sympy.I

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

def steadyStepError(system):
	fvt = eeval(sys2e(system), 0)[0]
	return {"sse": 1/(1+fvt)}

def reducedGainCompensate(system, target):
	w, mag, phase = bode(system, n=ct)
	phase = [i if i < 0 else -(360-i) for i in phase]
	phase = array(phase)
	phaseTarget = (-180 + target)
	print phaseMargin(system)
	print phaseTarget
	targetIndex = argmin(abs(phaseTarget - phase))
	print mag[targetIndex]
	magnitude = mag[targetIndex]
	return {"K": 1/magnitude}, 1/magnitude

def dominatePoleCompensate(system, target):
	w, mag, phase = bode(e2sys(sys2e(system)*1/s), n=ct)
	phase = [i if i < 0 else -(360-i) for i in phase]
	phase = array(phase)
	phaseTarget = (-180 + target)
	targetIndex = argmin(abs(phaseTarget - phase))
	print abs(phaseTarget-phase)
	magnitude = mag[targetIndex]
	return {"K": 1/magnitude}, (1/magnitude)*1/s

def lagCompensate(system, target):
	w, mag, phase = bode(system, n=ct)
	phase = [i if i < 0 else -(360-i) for i in phase]
	phase = array(phase)
	phaseTarget = (-180 + (target - 6)) # find location of new crossing goal; extra '-6' is from 10/wc rule
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

def characterize(system, f=figure(), color="k", labeled=""):
	""" draw bode plot for a given scipy.signal.lti instance. """
	adjustprops = dict(left=0.1, bottom=0.1, right=0.97, top=0.93, hspace=0.2)
	f.subplots_adjust(**adjustprops)
	magPlot = f.add_subplot(2, 1, 1)
	phasePlot = f.add_subplot(2, 1, 2, sharex=magPlot)
	w, mag, phase = bode(system)
	phase = [i if i < 0 else -(360-i) for i in phase]
	tImpulse, youtImpulse = impulse2(system)
	tStep, youtStep = step2(system)
	magPlot.semilogx(w, mag, label=labeled, color=color)
	magPlot.set_title("magnitude")
	magPlot.set_ylabel("amplitude")
	phasePlot.semilogx(w, phase, label=labeled, color=color)
	phasePlot.set_title("phase")
	phasePlot.set_xlabel("radians per second")
	phasePlot.set_ylabel("degrees")
	setp(magPlot.get_xticklabels(), visible=False)
	legend(loc="best")


print "before :"
G_p = 100/((s+1)*(0.1*s+1)*(0.01*s+1))
print phaseMargin(e2sys(G_p))
print steadyStepError(e2sys(G_p))
characterize(e2sys(G_p), color="k", labeled="uncompensated")
print "leadCompensate"
(coeffs, G_c) = leadCompensate(e2sys(G_p), 45)
print coeffs
print phaseMargin(e2sys(G_p*G_c))
print steadyStepError(e2sys(G_p*G_c))
characterize(e2sys(G_p*G_c), color="b", labeled="lead")
print "reducedGain"
(coeffs, G_c) = reducedGainCompensate(e2sys(G_p), 45)
print coeffs
print phaseMargin(e2sys(G_p*G_c))
print steadyStepError(e2sys(G_p*G_c))
characterize(e2sys(G_p*G_c), color='r', labeled="reducedGain")
print "majorPole"
(coeffs, G_c) = dominatePoleCompensate(e2sys(G_p), 45)
print coeffs
print phaseMargin(e2sys(G_p*G_c))
print steadyStepError(e2sys(G_p*G_c))
characterize(e2sys(G_p*G_c), color='g', labeled="majorPole")

show()
