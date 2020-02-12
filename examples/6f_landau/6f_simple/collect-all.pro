ni=collect(path='data',var='ni')
moment_xyzt,ni,rms=rmsni,dc=dcni
ti=collect(path='data',var='ti')
moment_xyzt,ti,rms=rmsti,dc=dcti
te=collect(path='data',var='te')
moment_xyzt,te,rms=rmste,dc=dcte
ps=collect(path='data',var='psi')
moment_xyzt,ps,rms=rmsps,dc=dcps
ph=collect(path='data',var='phi')
moment_xyzt,ph,rms=rmsph,dc=dcph
p=collect(path='data',var='p')
moment_xyzt,p,rms=rmsp,dc=dcp
jp=collect(path='data',var='jpar')
moment_xyzt,jp,rms=rmsjp,dc=dcjp
u=collect(path='data',var='u')
moment_xyzt,u,rms=rmsu,dc=dcu
vp=collect(path='data',var='vipar')
moment_xyzt,vp,rms=rmsvp,dc=dcvp

g=file_import("cbm18_dens8.grid_nx68ny64.nc")
;g=file_import("east077741.03500_psi080to105_x260y64.nc")

xx=41;190
yy=32;25
gr=deriv(alog(rmsp[xx,yy,*]))
psn=(g.psixy[*,yy]-g.psi_axis)/(g.psi_bndry-g.psi_axis)

window,1
@plot-mode

window,2
@plot-mode2
