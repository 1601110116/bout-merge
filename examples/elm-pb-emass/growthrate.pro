
p = collect(path='data', var='P')
moment_xyzt, p, rms=rms
siz = size(rms)
z = siz(3)
x = where(rms(*,32,z-1) eq max(rms(*,32,z-1)))
y = where(rms(x,*,z-1) eq max(rms(x,*,z-1)))
plot, deriv(alog(rms(x,y,*)))
;n0 = collect(path='data',var='n0')
print, deriv(alog(rms(x,y,*))); * sqrt(n0(x,y))
print, min(deriv(alog(rms(x,y,*))))

