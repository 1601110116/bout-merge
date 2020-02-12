g=file_import('outputx5y40.nc')
siz = size(g.rxy)
p1=plot(g.rxy(0,*), g.zxy(0,*))
for i=1, siz[1]-1 do p1=plot(g.rxy(i,*), g.zxy(i,*),/overplot)

h=file_import('outputx9y40v3.nc')
siz = size(h.rxy)
for i=0, siz[1]-1 do p2=plot(h.rxy(i,*), h.zxy(i,*), color='red', /overplot)

for i=0, 4 do print, (h.rxy(i+2,20)-h.rxy(4,20))/(g.rxy(i,20)-g.rxy(2,20))

