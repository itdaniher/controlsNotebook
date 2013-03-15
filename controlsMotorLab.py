from __future__ import division
import connectClient
import time
from scipy.interpolate import *
from pylab import *
import numpy

cee = connectClient.CEE()
cee.setSampleRate(40)
zero = 2.5
_unTwos = lambda x, bitlen: x-(1<<bitlen) if (x&(1<<(bitlen-1))) else x
_chunk = lambda l, x: [l[i:i+x] for i in xrange(0, len(l), x)] 
_flatten = lambda l: list(itertools.chain(*[[x] if type(x) not in [list] else x for x in l]))

def getCounters():
	data = cee.dev.ctrl_transfer(0x80 | 0x40, 0x10, 0, 0, 4)
	ticks = data[1] << 8 | data[0]
	samples = data[3] << 8 | data[2]
	return ticks, samples

def getRPS(dt = .1):
	def catchOverflow(a_v, b_v):
		if (a_v > (2**16)*(3/4)) and (b_v < (2**16)*(1/4)):
			b_v += 2**16
		return a_v, b_v
	a_t, a_s = getCounters()
	a_ss = cee.setOutputConstant('a', 'v', zero)['startSample']
	time.sleep(dt)
	b_t, b_s = getCounters()
	a_t, b_t = catchOverflow(a_t, b_t)
	a_s, b_s = catchOverflow(a_s, b_s)
	try:
		assert b_s > a_s
	except:
		return getRPS(dt)
	rotations = (b_t-a_t)/1088
	duration = (b_s-a_s)*cee.devInfo['sampleTime']
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
	cee.setOutputConstant('a', 'v', zero)
	dt = 0.2
	cee.setOutputConstant('b', 'v', 0)
	data = []
	for v in linspace(-zero, zero, 50):
		data.append(getDCSample(v, dt))
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
	cee.setOutputConstant('a', 'v', zero)
	tau4 = .5
	v = 10
	sampleCt = int(2*tau4/cee.devInfo['sampleTime'])
	ss = cee.setOutputArbitrary('b', 'i', [0, tau4/4, tau4/4, 2*tau4], [0, 0, +v, +v], repeat=-1)['startSample']
	data = []
	while True:
		data.append(getCounters())
		datum = data[-1]
		if datum[1] > ss+sampleCt:
			break
	(v, i) = cee.getInput('b', resample=0, count=sampleCt, start=ss)
	v = array(v) - zero
	s = arange(ss, sampleCt+ss)
	t = [d[1] for d in data]
	w = [d[0] for d in data]
	plotTwoAxes(s, v, i, "samples", "voltage", "current").show()
	x_f = linspace(t[0], t[-1], 100)
	fit = UnivariateSpline(t, w)
	figure()
	plot(x_f, fit(x_f), label="univariate spline fit")
	xlabel("time (samples)")
	ylabel("position")
	title("position over time")
	plot(t, w, '.', label="data")
	legend(loc="best")
	figure()
	plot(t[1::], diff(w)/diff(t), '.', label='numerically differentiated data')
	plot(x_f[1::], diff(fit(x_f))/diff(x_f), '-', label='derivative of interpolated and dt-normalized w')
	rpsps = diff(fit(x_f), 2)/(diff(x_f)[1::])
	semilogy(x_f[2::], rpsps, '-', label='second derivative of interpolated and dt-normalized w')
	xlabel("time (samples)")
	ylabel("data")
	legend(loc='best')
	show()

if __name__ == "__main__":
	print DCAnalysis()
