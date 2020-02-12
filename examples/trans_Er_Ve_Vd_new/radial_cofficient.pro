
g=file_import('CMod_1100223012_1150_260x64y_0.9psi_v1_diff_source.bout.nc')
filename = "CMod_1100223012_1150_260x64y_0.9psi_v1_diff_source.bout.nc"
mu = 4.e-7*3.1415926
Ni_x=1.
Te_x= 10.
Ti_x= 10.
density_unit = 1.e20                          ; 1/m^3
ee = 1.6022e-19                               ; elementary charge 

g11=g.Rxy*g.Rxy*g.Bpxy*g.Bpxy/g.rmag/g.rmag/g.bmag/g.bmag
J=(g.hthe/g.rmag)/(g.Bpxy/g.bmag)
lbar=g.rmag
unit_psi=g.rmag*g.rmag*g.bmag
jy=38
rxy=g.Rxy[*,jy]/lbar
bpxy=g.Bpxy[*,jy]/g.bmag
x=(g.psixy[*,jy]-g.psi_axis)/(g.psi_bndry-g.psi_axis)
xreal=g.psixy[*,jy]/unit_psi

NX=g.NX
NY=g.NY
jx=g.IXSEPS1
yup=g.JYSEPS1_1
ydown=g.JYSEPS2_2

ni=g.Niexp/Ni_x
ti=g.Tiexp/Ti_x
te=g.Teexp/Te_x

dnidx=deriv(xreal,ni[*,jy])
dtidx=deriv(xreal,ti[*,jy])
dtedx=deriv(xreal,te[*,jy])

print,'Gradient', '    Core','     separatrix','       Edge'
print,'dni : ', dnidx[0],'  ---',dnidx[jx],'  ---',dnidx[NX-1]
print,'ti_no_source : ', dtidx[0],'  ---',dtidx[jx],'  ---',dtidx[NX-1]
print,'te_no_source : ', dtedx[0],'  ---',dtedx[jx],'  ---',dtedx[NX-1]

print,'value', '    Ni','     Ti','       Te'
print,'x=0 value: ', ni[0,jy],'  ---',ti[0,jy],'  ---',te[0,jy]
print,'x=-1 value: ', ni[-1,jy],'  ---',ti[-1,jy],'  ---',te[-1,jy]

diff0 = (dnidx[0]*J[0,jy]*g11[0,jy])/(dnidx*J[*,jy]*g11[*,jy])
diff0=abs(diff0)
chi0 = dtidx[0]*J[0,jy]*g11[*,jy]/(dtidx*J[*,jy]*g11[*,jy])
che0 = dtedx[0]*J[0,jy]*g11[*,jy]/(dtedx*J[*,jy]*g11[*,jy])

print,'Setting PFR coefficients'
diff_pf=diff0
chi0_pf=chi0
che0_pf=che0
for i=jx-3,jx do diff_pf[i]=diff0[jx]
for i=jx-3,jx do chi0_pf[i]=chi0[jx]
for i=jx-3,jx do che0_pf[i]=che0[jx]
for i=0,NX-jx-1 do diff_pf[jx-3-i]=diff0[jx+i]
for i=0,NX-jx-1 do chi0_pf[jx-3-i]=chi0[jx+i]
for i=0,NX-jx-1 do che0_pf[jx-3-i]=che0[jx+i]
for i=0,2.*jx-NX-2 do diff_pf[i]=diff0[-1]
for i=0,2.*jx-NX-2 do chi0_pf[i]=chi0[-1]
for i=0,2.*jx-NX-2 do che0_pf[i]=che0[-1]

V_conv = diff0/ni[*,jy]*dnidx
window,1
plot,V_conv,chars=2,title='velorcity V_conv'

term_conv_e = te[*,jy]*deriv(xreal,V_conv)+3./2.*V_conv*dtedx
term_conv_i = ti[*,jy]*deriv(xreal,V_conv)+3./2.*V_conv*dtidx
Si=fltarr(NX)
Se=fltarr(NX)
for i=0,NX-1 do Si[i] = 0.0
for i=0,NX-1 do Se[i] = 0.0

for i=1,NX-1 do Si[i] = Si[i-1]-(term_conv_i[i]+term_conv_i[i-1])*(xreal[i]-xreal[i-1])/2.
for i=1,NX-1 do Se[i] = Se[i-1]-(term_conv_e[i]+term_conv_e[i-1])*(xreal[i]-xreal[i-1])/2.

window,2
plot,Si,chars=2,title='Source'
oplot,Se,color=2

chi = dtidx[0]*J[0,jy]*g11[0,jy]/(dtidx*J[*,jy]*g11[*,jy]) + Si*J[0,jy]*g11[0,jy]/(dtidx*J[*,jy]*g11[*,jy])
che = dtedx[0]*J[0,jy]*g11[*,jy]/(dtedx*J[*,jy]*g11[*,jy]) + Se*J[0,jy]*g11[0,jy]/(dtedx*J[*,jy]*g11[*,jy])
;for i=240,NX-1 do chi[i]=chi[240]
;for i=240,NX-1 do che[i]=che[240]
print,'Setting PFR coefficients with diffusion terms'
chi_pf=chi
che_pf=che
for i=jx-3,jx do chi_pf[i]=chi[jx]
for i=jx-3,jx do che_pf[i]=che[jx]
for i=0,NX-jx-1 do chi_pf[jx-3-i]=chi[jx+i]
for i=0,NX-jx-1 do che_pf[jx-3-i]=che[jx+i]
for i=0,2.*jx-NX-2 do chi_pf[i]=chi[-1]
for i=0,2.*jx-NX-2 do che_pf[i]=che[-1]

diff=fltarr(NX,NY)
chi_i=fltarr(NX,NY)
chi_e=fltarr(NX,NY)

;print,'calculate heat coefficient without diffusion terms'
;for i= 0,NY-1 do diff[*,i]=diff0
;for i= 0,NY-1 do chi_i[*,i]=chi0
;for i= 0,NY-1 do chi_e[*,i]=che0
;for i= 0,yup do diff[*,i]=diff_pf
;for i= 0,yup do chi_i[*,i]=chi0_pf
;for i= 0,yup do chi_e[*,i]=che0_pf
;for i= ydown+1,NY-1 do diff[*,i]=diff_pf
;for i= ydown+1,NY-1 do chi_i[*,i]=chi0_pf
;for i= ydown+1,NY-1 do chi_e[*,i]=che0_pf

print,'calculate heat coefficient with diffusion terms'
for i= 0,NY-1 do diff[*,i]=diff0
for i= 0,NY-1 do chi_i[*,i]=chi
for i= 0,NY-1 do chi_e[*,i]=che
for i= 0,yup do diff[*,i]=diff_pf
for i= 0,yup do chi_i[*,i]=chi_pf
for i= 0,yup do chi_e[*,i]=che_pf
for i= ydown+1,NY-1 do diff[*,i]=diff_pf
for i= ydown+1,NY-1 do chi_i[*,i]=chi_pf
for i= ydown+1,NY-1 do chi_e[*,i]=che_pf

print,'coefficient', '  separatrix'
print,'diff0 :  ', diff0[jx]
print,'chi :  ', chi[jx]
print,'che :  ', che[jx]

window,1
surface,diff,az=-20,chars=3,title='$D_{\perp}$'
window,2
surface,chi_i,az=-20,chars=3,title='$\chi_{i\perp}$'
window,3
surface,chi_e,az=-20,chars=3,title='$\chi_{e\perp}$'

window,4
plot,diff0,chars=2,title='perpendicular transport coefficient D'
window,5
plot,chi0,chars=2,title='chi_i'
oplot,chi,color=2
window,6
plot,che0,chars=2,title='chi_e'
oplot,che,color=2


handle = file_open(filename, /write)
s = file_write(handle, 'diff',diff)
s = file_write(handle, 'chi_i',chi_i)
s = file_write(handle, 'chi_e',chi_e)
file_close, handle
