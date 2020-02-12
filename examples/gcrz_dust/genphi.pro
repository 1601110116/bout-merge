pro genphi
;g=file_import("data/g056129.05550_2_x260y64_psi090to106_thermalSBC.nc")
;g=file_import("data/east67238.grd.nc")
;g=file_import("data/CMod_1100223012_1150_260x64y_0.9psi_v1_diff_Vi_para.bout.nc")
g=file_import("data/east056129.05550_psi085to105_x260y64_Vi_phi.nc")

for i=0,g.IXSEPS1-1 do begin
endfor

openw,lun,'data/prodata.dat',/get_lun
for j=0,g.ny-1 do begin
	for i=0,g.nx-1 do begin      
        printf,lun,g.rxy(i,j),g.zxy(i,j),g.tiexp(i,j),g.niexp(i,j),g.vi_par(i,j),g.phi(i,j),format='(6e15.6)'
	endfor
endfor
free_lun,lun

end
