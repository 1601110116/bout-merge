
p = collect(path='data', var='P')
b = fft(p, dimension=2)
moment_xyzt, b, rms=rms
siz = size(rms)
x = 514
y = 32
plot, deriv(alog(rms(x,0,*)))
;n0 = collect(path='data',var='n0')
print, deriv(alog(rms(x,0,*))); * sqrt(n0(x,y))
print, min(deriv(alog(rms(x,0,*))))

