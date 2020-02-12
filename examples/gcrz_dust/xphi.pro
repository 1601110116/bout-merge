;g=file_import("data/g056129.05550_2_x260y64_psi090to106_fitTe_decresePF_gaussianY.nc")
;g=file_import("data/g056129.05550_2_x260y64_psi090to106_fitTe_decresePF_gaussianY_RF_thermal.nc")
;g=file_import("data/g056129.05550_2_x260y64_psi070to106_fitTe_decreasePF_gaussianY.nc")
;g=file_import("data/g056129.05550_2_x260y64_psi070to106_fitTe_decreasePF_gaussianY.nc")
;g=file_import("data/g056129.05550_2_x260y64_psi090to106_thermalSBC.nc")
;g=file_import("data/g056129.05550_2_x260y64_psi090to106_thermalSBC_Li.nc")
;g=file_import("data/east67238.grd.nc")
g=file_import("data/east056129.05550_psi085to105_x260y64_Vi_phi.nc")
;g=file_import("data/CMod_1100223012_1150_260x64y_0.9psi_v1_diff_Vi_para.bout.nc")
;openr,lun,'data/xphi.dat',/get_lun
;openr,lun,'data/xphi_nocore2.dat',/get_lun
;openr,lun,'data/xphi_090106RFSBC.dat',/get_lun
openr,lun,'data/errortest.dat',/get_lun
;openr,lun,'data/xphi_090106SBC_Li_constantcore.dat',/get_lun
;openr,lun,'data/xphi_090106SBC_Li_analytical3_fine.dat',/get_lun
data=fltarr(5,g.nx,g.ny)
print,data(0:4,0,0)
print,data(0:4,g.nx-1,g.ny-1)

readf,lun,data
; (*gRxy)(i,j), (*gZxy)(i,j), (*gphish0)(i,j),local value,  local error

print,data(0:4,0,0)
print,data(0:4,1,0)
print,data(0:4,g.nx-1,g.ny-1)
print,data(0:4,100,60)

maxphi=max(abs(data(2,*,*)))
print,"The maximum local err and its x,y index are:"
print, max(abs(data(4,*,*)),index),index mod g.nx, index/g.nx;
print,"The minimum local err and its x,y index are:"
print, min(abs(data(4,*,*)),index),index mod g.nx, index/g.nx;
surface, data(4,*,*),charsize=2

globalerr=(data(2,*,*)-data(3,*,*))/maxphi
print,"The maximum global err and its x,y index are:"
print, max(abs(globalerr),index),index mod g.nx, index/g.nx;
print,"The minimum global err and its x,y index are:"
print, min(abs(globalerr),index),index mod g.nx, index/g.nx;

;surface, (data(2,*,*)-data(3,*,*))/maxphi,charsize=3

