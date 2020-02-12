
;;file_list in trans_neut module
;J psixy aveY_J aveY_g11 aveY_g11J eta_spitzer DDX_Ni DDX_Ti DDX_Te 
;Diffc_ni_step chic_i_step chic_e_step Diffc_ni_Hmode chic_i_Hmode chic_e_Hmode Diff_perp_couple chi_i_perp_couple
;chi_e_perp_couple Diff_ni_perp chi_i_perp chi_e_perp kappa_Te kappa_Ti Diffc_nn_perp 
;nu_ionz nu_CX nu_ionz_n nu_CX_n nu_diss nu_rec Si_p Si_CX S_diss S_rec 
; Grad_par_pei Grad_par_pn Grad_par_logNn tau_ei q_se_kappa q_si_kappa Gamma_ni_Diffc
;Gamma_nn_Diffc Flux_Ni_aveY Flux_Ti_aveY Flux_Te_aveY Ni Vi Te Ti Nn Tn Vn Nm Vmx Vmy Vmz

print,'Collecting gfile '
g=file_import('cbm18_8_y064_x516_090309.nc')

path='data'

density_unit = 1.e20   ;;; 1./m^3
ee=1.6e-19             ;; SI unit electron charge
punit=density_unit*ee  ;; SI unit, Pascal  

print,'Collecting data and set Path =',path

print,'Collecting parameters '

tebar=collect(path=path,var='te_x')
lbar=collect(path=path,var='lbar')
tbar=collect(path=path,var='tbar')
t=collect(path=path,var='t_array')
print,'Running time is [ms]'
print,t*tbar*1.e3   ;;in unit ms

diunit=lbar*lbar/tbar
print,'unit_diffusion coefficient is', diunit,'m^2/s'

var='ni'
print,'collecting data ', var
f=collect(path=path,var=var)
ni=f
var='ti'
print,'collecting data ', var
f=collect(path=path,var=var)
ti=f

var='te'
print,'collecting data ', var
f=collect(path=path,var=var)
te=f
var='vi'
print,'collecting data ', var
f=collect(path=path,var=var)
vi=f
var='Diff_perp_couple'
print,'collecting data ', var
f=collect(path=path,var=var)
di=f
var='chi_i_perp_couple'
print,'collecting data ', var
f=collect(path=path,var=var)
chi=f
var='chi_e_perp_couple'
print,'collecting data ', var
f=collect(path=path,var=var)
che=f

dir=di*diunit
chir=chi*diunit
cher=che*diunit

pei=ni*(ti+te)*tebar
peir=pei*punit



safe_colors,/first
