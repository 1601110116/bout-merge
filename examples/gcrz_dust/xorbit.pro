;g=file_import("data/g056129.05550_2_x260y64_psi090to106_fitTe_decresePF_gaussianY_RF_thermal.nc")
;g=file_import("data/g056129.05550_2_x260y64_psi090to106_thermalSBC.nc")
;g=file_import("data/g056129.05550_2_x260y64_psi070to106_fitTe_decreasePF_gaussianY.nc")
;g=file_import("data/CMod_1100223012_1150_260x64y_0.9psi_v1_diff_Vi_para.bout.nc")
g=file_import("data/east056129.05550_psi085to105_x260y64_Vi_phi.nc")

;r1=collect(path='data',var='r1')
;z1=collect(path='data',var='z1')
;energy1=collect(path='data',var='energy1')
;Pzeta1=collect(path='data',var='pzeta1')
td=collect(path='data',var='td')
rd=collect(path='data',var='rd')
;bb=collect(path='data',var='bb')
phi=collect(path='data',var='phi')
r=collect(path='data',var='r')
z=collect(path='data',var='z')
vpara=collect(path='data',var='vpara')
;pzeta=collect(path='data',var='pzeta')
;region=collect(path='data',var='region')
;energy=collect(path='data',var='energy')
time1=collect(path='data',var='time1')
time2=collect(path='data',var='time2')
ti1=collect(path='data',var='ti1');
;te1=collect(path='data',var='te1');
ni1=collect(path='data',var='ni1');
;ne1=collect(path='data',var='ne1');
;test1=collect(path='data',var='test1');
;test2=collect(path='data',var='test2');
;test3=collect(path='data',var='test3');
;test4=collect(path='data',var='test4');
;test5=collect(path='data',var='test5');
;test6=collect(path='data',var='test6');
;test7=collect(path='data',var='test7');
;test8=collect(path='data',var='test8');
test11=collect(path='data',var='test11');
test12=collect(path='data',var='test12');
test13=collect(path='data',var='test13');
test14=collect(path='data',var='test14');
test15=collect(path='data',var='test15');
test16=collect(path='data',var='test16');
test17=collect(path='data',var='test17');
test18=collect(path='data',var='test18');
;test21=collect(path='data',var='test21');
;test22=collect(path='data',var='test22');
;test23=collect(path='data',var='test23');
;test24=collect(path='data',var='test24');
;test25=collect(path='data',var='test25');
;test26=collect(path='data',var='test26');
;test27=collect(path='data',var='test27');
;test28=collect(path='data',var='test28');
;ti0=collect(path='data',var='ti0');
;ni0=collect(path='data',var='ni0');
vp=collect(path='data',var='vp');
;kusai=collect(path='data',var='kusai');
;kusaia=collect(path='data',var='kusaia');
;kusaie=collect(path='data',var='kusaie');
;gammaa=collect(path='data',var='gammaa');
;gammae=collect(path='data',var='gammae');
meanre=collect(path='data',var='meanre');
deltath=collect(path='data',var='deltath');
deltasec=collect(path='data',var='deltasec');
vd=collect(path='data',var='vd');
aa=collect(path='data',var='aa');
zz=collect(path='data',var='zz');
;g.rxy=g.rxy/g.rmag
;g.zxy=g.zxy/g.rmag
r=r*g.rmag
z=z*g.rmag

nstep=3000
;print,energy1(1:nstep)/energy1(1)
;print,pzeta1(1:nstep)/pzeta1(1)
;!p.charsize=2
plot,g.rxy(259,*),g.zxy(259,*),psym=4,color=0,background=1,/iso
oplot,g.rxy(0,*),g.zxy(0,*),psym=4,color=0
oplot,g.rxy(*,0),g.zxy(*,0),psym=4,color=0
oplot,g.rxy(*,63),g.zxy(*,63),psym=4,color=0
oplot,g.rxy(g.IXSEPS1-1,*),g.zxy(g.IXSEPS1-1,*),psym=3,color=0
oplot,r(1:nstep-1),z(1:nstep-1),thick=3,color=2


;oplot,energy1(1:nstep-1),yrange=[-0.1,1.3],thick=3,charsize=2.0
;oplot,pzeta1(1:nstep-1),thick=3
