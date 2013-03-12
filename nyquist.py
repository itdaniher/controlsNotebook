import scipy.signal as signal
from pylab import *
hold(True)

# L(s) = K/s
# K/(K+s)

for i in logspace(-10,10):
	w, h = signal.freqresp(((i), (1, i)))
	plot(h.real, h.imag, 'b', label=str(i))
	plot(h.real, -h.imag, 'r', label=str(i))
show()
