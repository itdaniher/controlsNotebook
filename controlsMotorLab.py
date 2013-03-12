from __future__ import division
import connectClient
import time

from pylab import *
import numpy

cee = connectClient.CEE()

_unTwos = lambda x, bitlen: x-(1<<bitlen) if (x&(1<<(bitlen-1))) else x
_chunk = lambda l, x: [l[i:i+x] for i in xrange(0, len(l), x)] 
_flatten = lambda l: list(itertools.chain(*[[x] if type(x) not in [list] else x for x in l]))

def getCounters():
	data = cee.dev.ctrl_transfer(0x80 | 0x40, 0x10, 0, 0, 4)
	ticks = _unTwos(data[1] << 8 | data[0], 16)
	ticks = data[1] << 8 | data[0]
	samples = data[3] << 8 | data[2]
	return ticks, samples

def getRPS():
	def catchOverflow(a_v, b_v):
		if (a_v > (2**16)*(3/4)) and (b_v < (2**16)*(1/4)):
			b_v += 2**16
		return a_v, b_v
	a_t, a_s = getCounters()
	time.sleep(.01)
	b_t, b_s = getCounters()
	a_t, b_t = catchOverflow(a_t, b_t)
	a_s, b_s = catchOverflow(a_s, b_s)
	rps = ((b_t-a_t)/1024)/(b_s-a_s)*1.0/cee.devInfo['sampleTime']
	return rps

def setMotorVoltage(v):
	if v > 0:
		cee.setOutputConstant('b', 'v', abs(v))
		return cee.setOutputConstant('a', 'v', 0.0)
	if v < 0:
		cee.setOutputConstant('b', 'v', 0)
		return cee.setOutputConstant('a', 'v', abs(v))
	else:
		cee.setOutputConstant('b', 'v', 0)
		return cee.setOutputConstant('a', 'v', 0)
	
def setMotorCurrent(i):
	if i > 0:
		cee.setOutputConstant('b', 'i', i)
		cee.setOutputConstant('a', 'v', 0)
	if i < 0:
		cee.setOutputConstant('b', 'v', 0)
		cee.setOutputConstant('a', 'i', i)
	else:
		cee.setOutputConstant('a', 'd', 0)
		cee.setOutputConstant('b', 'v', 0)

print setMotorVoltage(-2)

while True:
	print getRPS()
