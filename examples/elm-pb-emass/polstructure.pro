
restore, f='p.idl.dat'

p10 = (p(512,*,8,10))
p10 = p10 / max(abs(p10))
p50 = (p(512,*,8,50))
p50 = p50 / max(abs(p50))

z10 = p(512,32,*,10)
z10 = z10 / max(abs(z10))
z50 = p(512,32,*,50)
z50 = z50 / max(abs(z50))

b10 = abs(b(512,*,8,10))
b10 = b10 / max(abs(b10))
b50 = abs(b(512,*,8,50))
b50 = b50 / max(abs(b50))


p1 = plot(p10, color='black', name='t=10', xtitle='y', ytitle='p', font_size=18)
p2 = plot(p50, color='red', name='t=50', /overplot)
l1 = legend(target=[p1,p2], font_size=18)

p3 = plot(b10, color='black', name='t=10', xtitle='FFT(y)', ytitle='p', font_size=18)
p4 = plot(b50, color='red', name='t=50', /overplot)
l2 = legend(target=[p3, p4], font_size=18)

p5 = plot(z10, color='black', name='t=10', xtitle='z', ytitle='p', font_size=18)
p6 = plot(z50, color='red', name='t=50', /overplot)
l3 = legend(target=[p5,p6], font_size=18)

p7 = plot(deriv(alog(rmsp(512,32,*))), color='black', name='before FFT', xtitle='t', ytitle='$\gamma$', font_size=18)
p8 = plot(deriv(alog(rmsb(512,0,*))), color='red', name='after FFT', /overplot)
l4 = legend(target=[p7,p8], font_size=18)

