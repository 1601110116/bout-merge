
restore, filename='rms_68.dat'
restore, filename='rms_132.dat'

x68 = (indgen(68) - 1) / 68.
x132 = (indgen(132) - 1) /132.

plot, x68, rms_68(34,*) / max(rms_68(34,*))

