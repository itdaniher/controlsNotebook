from __future__ import division
import connectClient
import time
from scipy.interpolate import *
from pylab import *
import numpy

cee = connectClient.CEE()
# set a sample rate of 40,000 samples per second
cee.setSampleRate(40)
# for four quadrant operation, use a "zero" point of 2v5
zero = 2.5

def getCounters():
	# get four bytes from xmega onboard cee
	data = cee.dev.ctrl_transfer(0x80 | 0x40, 0x10, 0, 0, 4)
	# ticks are time increments of units samples
	sampleCounter = data[1] << 8 | data[0]
	position = data[3] << 8 | data[2]
	return sampleCounter, position

def getRPS(dt = .1):
	# catch overflow of 16b timer/counters by comparing two subsequent values
	def catchOverflow(a_v, b_v):
		if (a_v > (2**16)*(3/4)) and (b_v < (2**16)*(1/4)):
			b_v += 2**16
		return a_v, b_v
	# getCounters
	a_t, a_s = getCounters()
	# wait
	time.sleep(dt)
	# getCounters
	b_t, b_s = getCounters()
	# clean
	a_t, b_t = catchOverflow(a_t, b_t)
	a_s, b_s = catchOverflow(a_s, b_s)
	try:
		assert b_s > a_s
	except:
		return getRPS(dt)
	# normalize position data
	rotations = (b_t-a_t)/1088
	# normalize time data
	duration = (b_s-a_s)*cee.devInfo['sampleTime']
	# rps = rotations / duration
	rps = rotations / duration
	return rps


def getDCSample(v, dt = .1):
	ss = cee.setOutputConstant('b', 'v', zero+v)['startSample']
	(v, i) = cee.getInput('b', resample=dt, count=1, start = ss+int(dt/cee.devInfo['sampleTime']))
	v = v - zero
	dw = getRPS(dt)
	return {'v':v, 'i':i, 'dw':dw}

def plotTwoAxes(x, y1, y2, xlabel="x", y1label="y1", y2label="y2"):
	f = figure()
	hold(True)
	ax1 = f.add_subplot(111)
	ax1.plot(x, y1, 'r.')
	ax1.set_xlabel(xlabel)
	ax1.set_ylabel(y1label, color="r")
	for tl in ax1.get_yticklabels():
		tl.set_color('r')
	ax2 = ax1.twinx()
	ax2.plot(x, y2, 'b.')
	ax2.set_ylabel(y2label, color="b")
	for tl in ax2.get_yticklabels():
		tl.set_color('b')
	return f

def DCAnalysis():
	# ensure "zero" point is where we want it to be
	cee.setOutputConstant('a', 'v', zero)
	# set reasonable timestep for steady state analysis
	dt = 0.2
	# set to max negative in light of future sampling
	cee.setOutputConstant('b', 'v', 0)
	data = []
	# go through len 50 list of voltages, get voltage, current, and rotational velocity for each voltage
	for v in linspace(-zero, zero, 50):
		data.append(getDCSample(v, dt))
	# cleanup
	v = [d['v'] for d in data]
	i = [d['i'] for d in data]
	dw = [d['dw'] for d in data]
	v = array(v)
	i = array(i)
	dw = array(dw)
	f = figure()
	hold(True)
	ax1 = f.add_subplot(111)
	ax1.plot(v, dw, 'r.', label="dw/dt data (rps)")
	legend(loc="best")
	ax1.set_xlabel("voltage (v)")
	ax1.set_ylabel("rotations per second", color="r")
	for tl in ax1.get_yticklabels():
		tl.set_color('r')
	ax2 = ax1.twinx()
	ax2.plot(v, i, 'bo', label="measured current (mA)")
	ax2.set_ylabel("current draw (mA)", color="b")
	for tl in ax2.get_yticklabels():
		tl.set_color('b')
	resistance = polyfit(v, i, 1)
	ax2.plot(v, polyval(resistance, v), '-', label="resistance fit")
	legend(loc='best')
	f.savefig("DCanalysis.png")	
	print resistance[0], " ohms"
	i = -i
	k_e = polyfit(v - i*resistance[0], dw, 1)[0]
	print k_e, "rps per volt"
	k_e = polyfit(v, dw, 1)[0]
	print k_e, "rps per volt"
	legend(loc='best')
	#show()
	return data

def ACAnalysis():
	# make sure zero is actually zero
	cee.setOutputConstant('a', 'v', zero)
	# total observable behavior should span 4 tau
	# tau is somewhat arbitrarily chosen constant for what seemed to ecapsulate the interesting parts
	tau = .25
	# value of 10, in this case, 10mA
	v = 10
	# calculate how many samples are contained in four timesteps
	sampleCt = int(4*tau/cee.devInfo['sampleTime'])
	# set a step from 0mA to 10mA to happen at tau/4, measure until 4tau
	# "ss" is the integer value of the starting sample
	ss = cee.setOutputArbitrary('b', 'i', [0, tau/4, tau/4, 4*tau], [0, 0, +v, +v], repeat=0)['startSample']
	# instantiate empty array
	data = []
	while True:
		# inner loop, simply call getCounters, shove the data into the array, break if the sample count is more than our target end point
		data.append(getCounters())
		datum = data[-1]
		if datum[1] > ss+sampleCt:
			break
	# get 'sampleCt' samples into lists "v" and "i" with no resampling, starting at the sample point the arbitrary waveform started
	(v, i) = cee.getInput('b', resample=0, count=sampleCt, start=ss)
	# normalize to "zero"
	v = array(v) - zero
	# generate array of "sampleCt" sample indexes
	s = arange(ss, sampleCt+ss)
	# "t" or ticks is the 2nd element in data
	t = [d[1] for d in data]
	# "w" or omega is the 1st element
	w = [d[0] for d in data]
	# plot motor voltage and current on the same plot
	plotTwoAxes(s, v, i, "samples", "voltage", "current").show()
	# generate an abbreviated set of times on the continuum from the first measured sample count to the last
	x_f = linspace(t[0], t[-1], 100)
	# fit a spline to our rotations over time data
	fit = UnivariateSpline(t, w)
	# show quality of fit
	figure()
	plot(x_f, fit(x_f), label="univariate spline fit")
	xlabel("time (samples)")
	ylabel("position")
	title("position over time")
	plot(t, w, '.', label="data")
	legend(loc="best")
	# use spline as low-jitter source of rotational data capable of being numerically integrated
	figure()
	plot(t[1::], diff(w)/diff(t), '.', label='numerically differentiated data')
	plot(x_f[1::], diff(fit(x_f))/diff(x_f), '-', label='derivative of interpolated and dt-normalized w')
	rpsps = diff(fit(x_f), 2)/(diff(x_f)[1::])
	semilogy(x_f[2::], rpsps, '-', label='second derivative of interpolated and dt-normalized w')
	xlabel("time (samples)")
	ylabel("data")
	legend(loc='best')
	show()
