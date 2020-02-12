



; file_list in the dmp.nc files
;phi0 J0 P0 B0 Dphi0 U0 V0 Ti0 Te0 N0 eta vexb_x vexb_y vexb_z vbtild_x
;Di_couple kaii_couple kaie_couple jpar phi U P Psi

print,'Collecting gfile '
g=file_import('cbm18_8_y064_x516_090309.nc')

path='data'

print,'Collecting data and set Path =',path

print,'Collecting parameters '

tbar=collect(path=path,var='tbar')
t=collect(path=path,var='t_array')
bbar=collect(path=path,var='bbar')

print,'Defining some  parameters '
PI=3.1415926
mu0=4.0e-7*PI

punit= bbar*bbar/2./mu0

p0=g.pressure

print,'Collecting P'
p=collect(path=path,var='P')
moment_xyzt,p,rms=rmsp,dc=dcp

sdcp=size(dcp)
sp0=size(p0)

;;reform the 2D p0 to 3D array of size same as dcp array
p03=rebin(reform(p0,sp0[1],sp0[2],1),sdcp[1],sdcp[2],sdcp[3])
dcpr=dcp*punit     ;; change to the SI unit
ptlt=dcpr+p03      ;; total pressure

;print,'Collecting jpar'
;jpar=collect(path=path,var='jpar')
;print,'Collecting U'
;u=collect(path=path,var='U')

;print,'Collecting Psi'
;psi=collect(path=path,var='Psi')
;print,'Collecting Di_couple'
;di=collect(path=path,var='Di_couple')
;print,'Collecting Chii_couple'
;chii=collect(path=path,var='kaii_couple')
;print,'Collecting Chie_couple'
;chie=collect(path=path,var='kaie_couple')

jy0=32
nt0=0
dnt=20
fchars=2
fthick=2
fcharthick=2
  window,1
;string(t,m,s,format='(i02,":",i02,":",i02)'
f=ptlt
ftitle=string('P0+P.DC at different times for dnt=',STRTRIM(dnt))
;ftitle=STRJOIN('P0+P.DC at different times for dnt=',string(dnt))

fxtitle='jx'
  plot,f[*,jy0,nt0],chars=fchars,charthick=fcharthick,thick=fthick,title=ftitle,xtitle=fxtitle, $
       /xst
  oplot,f[*,jy0,nt0+dnt],thick=fthick,color=3
  oplot,f[*,jy0,nt0+2*dnt],thick=fthick,color=4
  oplot,f[*,jy0,nt0+3*dnt],thick=fthick,color=5
  oplot,f[*,jy0,nt0+4*dnt],thick=fthick,color=6
  oplot,f[*,jy0,nt0+5*dnt],thick=fthick,color=2

