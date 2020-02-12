



; file_list in the dmp.nc files
;phi0 J0 P0 B0 Dphi0 U0 V0 Ti0 Te0 N0 eta vexb_x vexb_y vexb_z vbtild_x
;Di_couple kaii_couple kaie_couple jpar phi U P Psi

path='data'

print,'Collecting data and set Path =',path

print,'Collecting P'
p=collect(path=path,var='P')
print,'Collecting jpar'
jpar=collect(path=path,var='jpar')
print,'Collecting U'
u=collect(path=path,var='U')

print,'Collecting Psi'
psi=collect(path=path,var='Psi')
print,'Collecting Di_couple'
di=collect(path=path,var='Di_couple')
print,'Collecting Chii_couple'
chii=collect(path=path,var='kaii_couple')
print,'Collecting Chie_couple'
chie=collect(path=path,var='kaie_couple')
