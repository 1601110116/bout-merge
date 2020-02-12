FUNCTION Edge_Collisionality, 	n_e=neIn, n_Norm=n_Norm, T_e=TeIn, T_Norm=T_Norm, T_i=Tiin, 	$
				Z_eff=Zeff, mu=mu_in, R_o=Ro, a=aIn, q=qin, BpBt=BpBt_in, 	$
				B_o=Bin, Omega_ci = Omega_in, chi=chiIn
;
; Purpose:
;
;	This function computes various quantitiies related to the edge collisionality
;	and returns them in a structure.
;
; Input Keywords:
;
;	n_e	A GKVs1D object containing the electron density in units of cm^-3
;		vs. radius.  The radial coordinate should be in units of cm.
;		n_e MUST be provided. NOT OPTIONAL!
;
;	n_Norm	In the event that the GKVs1D object containing n_e is not in units
;		of cm^-3, set n_Norm to a constant such that n_e*n_Norm yields 
;		the electtron density in the correct units.  Defaults to 1.0.
;		(Optional)
;
;	T_e	A GKVs1D object containing the electron temperature in units of eV
;		vs. radius.  The radial coordinate should be in units of cm.
;		T_e MUST be provided.  NOT OPTIONAL!
;
;	T_norm	In the event that the GKVs1D object containing T_e is not in units
;		of eV, set T_Norm to a constant such that T_e*T_Norm yields 
;		the electtron temperature in the correct units.  Defaults to 1.0.
;		(Optional)
;
;	T_i	A GKVs1D object containint the ion temperature in the SAME units in
;		which T_e is provided, such that T_i*T_Norm yields the ion temperature
;		in eV vs. radius.  The radial coordinate should be in units of cm.
;		Defaults to T_e.  (Optional)
;
;	Z_eff	The effective charge.  Defaults to 1.0. (Optional)
;
;	mu	The ion mass in units of the atomic mass.  Defaults to 2.0.
;		(Optional)
;
;	R_o	The major radius of the device in cm.  Defaults to 100.
;
;	a	The minor radius of the device in cm.  Devaults to 20.
;		(Optional)
;
;	q	The safety factor at the plasma edge.  This is used to estimate the 
;		The connection length, R_o*q.  Defaults to 3.0.  (Optional)
;
;	BpBt	Alternatively, the user may specify the ratio of the poloidal field
;		to the toroidal field at the outboard midplane.  The connection
;		length is then estimated as a/BpBt.  (Optional)
;
;	B_o	The magnitude of the mangetic field in Gauss.  Defaults to 10^4.
;		(Optional).
;
;   Omega_ci	Alternatively, you can specify the ion cyclotron frequency, which will
;		be used to compute B_o.  (Optional)
;
;	chi	The perpendicular heat conductivity in m^2/s.  This is used to estimate
;		The energy confinement time against perpendicular conduction losses only.
;		Defaults to 1 m^2/s. (Optional)
;
;
;  Output:
;
;	This function returns a structure containing:
;
;		nuStar_e	The ratio of the neoclassical effective electron
;				collision frequency to the transit frequency
;				vs. radius.  Calculated after Sauter et al,
;				Phys. Plasmas <6> 7, 2834 (July 1999).
;
;		nuStar_i	The ratio of the neoclassical effective ion
;				collision frequency to the transit frequency
;				vs. radius.  Calculated after Sauter et al,
;				Phys. Plasmas <6> 7, 2834 (July 1999).
;
;		l_mfp_e		The electron mean-free-path vs. radius.
;
;		l_mfp_i		The ion mean-free path vs. radius.
;
;		L_resistive	The characteristic resistive length vs. radius.
;
;		Normed_L_resistive	Ratio of the resistive length to the
;					connection length.
;
;		nu_Energy	The energy exchange rate vs. radius.
;
;		tau_ll		The electron energy confinement time against
;				heat conduction along the field line.
;
;		tau_perp	The electron energy confinement time against
;				heat conduction across the field.
;
;		Normed_Tau_ll	Product of tau_ll and the energy exchange time.
;
;		Normed_tau_perp	Produce of tau_perp and the energy exchange time.
;
;		kappa_n		(1/n_e)(dn_e/dx) vs. radius.
;
;		kappa_T		(1/T_e)(dT_e/dx) vs. radius.
;
;		kappa_Ti	(1/T_i)(dT_i/dx) vs. radius.
;
;		Normed_Kappa_n	(rho_p/n_e)(dn_e/dx), where rho_p is the local ion
;				poloidal gyroradius.
;
;		Normed_Kappa_T	(rho_p/T_e)(dT_e/dx),where rho_p is the local ion
;				poloidal gyroradius.
;
;		Normed_Kappa_Ti	(rho_p/T_i)(dT_i/dx),where rho_p is the local ion
;				poloidal gyroradius.
;
;		j_Bootstrap	The neoclassical bootstrap current with collisional
;				corrections as per Sauter et al, Phys. Plasmas <6> 7, 
;				2834 (July 1999).  Units are arbitrary!  Just for 
;				profile information.
;
;	
;
;  Written by W.M. Nevins
;	10/23/01
;  Revised by W.M. Nevins
;	10/27/02
;  
n_e= neIn	; values assumed input in units of /cm^3
T  = TeIn	; values assumed input in units of eV
IF(N_ELEMENTS(Tiin) EQ 0) THEN BEGIN
	Ti=Tein
ENDIF ELSE BEGIN
	Ti=T
ENDELSE

Z_eff = 1.0
IF(N_ELEMENTS(Zeff) GT 0) THEN Z_eff=Zeff

mu = 2.
IF(N_ELEMENTS(muin) GT 0) THEN mu=muin

R_o=100.
IF(N_ELEMENTS(Ro) GT 0) THEN R_o=Ro

a=20
IF(N_ELEMENTS(aIn) GT 0) THEN a=aIn

AspectRatio=R_o/a
epsilon = a/R_o

q=3.0
IF(N_ELEMENTS(qin) GT 0) THEN q=qin

BpBt = epsilon/q
IF(N_ELEMENTS(BpB_in) EQ 1) THEN BpBt=BpBt_in

connection_Length = a/BpBT


Omega_ci = 9.58e7/mu		; Rad/s
B_o = 1.e4	; Gauss
IF(N_ELEMENTS(Bin) GT 0) THEN BEGIN
	B_o=Bin
	Omega_ci = 9.58e3*B_o/mu
ENDIF ELSE BEGIN
	IF(N_ELEMENTS(Omega_In) GT 0) THEN BEGIN
		Omega_ci = Omega_In
		B_o = mu*Omega_In/9.58e3
	ENDIF
ENDELSE

chi=1.e4			; Perp thermal conductivity in cm^2/s
IF(N_ELEMENTS(chiIn) GT 0) THEN chi=1.e4*chiIn	; (assume input chi is in m^2/s, and convert to cm^2/s


T   -> Get, values=Tptr		; eV
Ti  -> Get, values=TiPtr	; ev
n_e -> Get, values=nptr		; in units of /cm^3

nValues = (*nPtr)
Tvalues = *Tptr
TiValues= *TiPtr

IF(N_ELEMENTS(n_Norm) EQ 1) THEN nValues=nValues*n_Norm
IF(N_ELEMENTS(T_Norm) EQ 1) THEN BEGIN
	Tvalues = Tvalues*T_Norm
	TiValues = TiValues*T_Norm
ENDIF

Ln_Lambda_e = 24. - ALOG( SQRT(nValues)/Tvalues )
Ln_Lambda_i = 23. - ALOG( SQRT(nValues*Z_eff/TiValues)/TiValues )

nu_e = 2.91e-6*(0.5*(1+Z_eff))*(nValues)*Ln_Lambda_e/( Tvalues^1.5)
nu_i = 4.80e-8*(       Z_eff )*(nValues)*Ln_Lambda_i/(TiValues^1.5)/SQRT(mu)

v_te = 4.19e7*SQRT(Tvalues)	; cm/s
v_ti = 9.79e5*SQRT(TiValues/mu)	; cm/s
c_s  = 9.79e5*SQRT(TValues/mu)	; cm/s

l_mfp_e = v_te/nu_e		; cm
l_mfp_i = v_ti/nu_i		; cm

l_mfp_e_Obj = n_e -> MAKECOPY(/novalues, /noerrorbars)
l_mfp_i_Obj = n_e -> MAKECOPY(/novalues, /noerrorbars)

l_mfp_e_Obj -> Set, values=PTR_NEW(l_mfp_e)
l_mfp_i_Obj -> Set, values=PTR_NEW(l_mfp_i)

l_mfp_e_Obj -> Set, title="!12l!X!S!U(e)!R!Dmfp!N", mnemonic='l_mfp_e', units='cm'
l_mfp_i_Obj -> Set, title="!12l!X!S!U(i)!R!Dmfp!N", mnemonic='l_mfp_i', units='cm'

nu_Energy = 3.2e-9*(nValues)*Z_eff*Ln_Lambda_e/( mu*Tvalues*SQRT(Tvalues) )	; Hz

nu_Energy_Obj = n_e -> MAKECOPY(/novalues, /noerrorbars)
nu_Energy_Obj -> set, values=PTR_NEW(nu_Energy)
nu_Energy_Obj -> Set, title='!4m!X!DE!N', mnemonic='nu_Energy', units='Hz'

tau_e = 3.44e5*Tvalues*SQRT(Tvalues)/(nValues)/Ln_Lambda_e		; s
kappa_ll = 3.2*v_te*v_te*tau_e						; cm^2/s
tau_ll = (49./12)*(connection_length)^2/kappa_ll					; s
tau_ll_alt = (3./7.)*connection_length/c_s				; s
tall_ll = tau_ll + tau_ll_alt						; s

tau_ll_Obj = n_e -> MAKECOPY(/novalues, /noerrorbars)
tau_ll_Obj -> Set, values=PTR_NEW(tau_ll)
tau_ll_Obj -> Set, title='!4s!S!D!9#!R!U!X(E)!N', mnemonic='tau_ll', units='s'

Normed_tau_ll_Obj = tau_ll_Obj -> times(nu_Energy_Obj)
Normed_tau_ll_Obj -> Set, title='!4s!S!D!9#!R!U!X(E)!N'+'!4m!X!DE!N', mnemonic='tau_ll_normed', units=''

;nuStar_e = q*(100.*R_o)*nu_e/v_te
;nuStar_i = q*(100.*R-o)*nu+i/v_ti

; following is after Sauter et al, Phys. Plasmas <6> 7, 2834 (July 1999) -- see pg. 2838

nuStar_e = 6.921e-20*connection_Length*(nValues*1.e6)*Z_eff*Ln_Lambda_e/(Tvalues^2*epsilon^1.5)
nuStar_i = 4.900e-20*connection_Length*(nValues*1.e6)*Z_eff*Ln_Lambda_i/(TiValues^2*epsilon^1.5)

nuStar_e_Obj = n_e -> MAKECOPY(/novalues, /noerrorbars)
nuStar_i_Obj = n_e -> MAKECOPY(/novalues, /noerrorbars)

nuStar_e_Obj -> SET, title='!4m!X!S!U*!R!De!N', mnemonic='nuStar_e', values=PTR_NEW(nuStar_e), units=''
nuStar_i_Obj -> SET, title='!4m!X!S!U*!R!Di!N', mnemonic='nuStar_i', values=PTR_NEW(nuStar_i), units=''

rho_i = 1.02e2*SQRT(mu*TiValues)/(B_o)	; cm
rho_p = rho_i/BpBt			; cm

dndx = n_e -> dbyd('x')
temp = dndx -> over(n_e)
temp -> Get, values=tempPtr
dLogNdxValues = *tempPtr
dndx -> trash
kappa_n = temp -> times(-1.)
kappa_n -> set, title='!4j!X!Dn!N', mnemonic='kappa_n', units='1/cm'

Normed_kappa_n = temp -> times(-rho_p)
temp -> trash
Normed_kappa_n -> set, title='!4j!X!Dn!N', mnemonic='kappa_n', units='1/!4q!X!Dp!N'

dTdx = T -> dbyd('x')
temp1 = dTdx -> over(T)
temp1 -> Get, values=temp1Ptr
dLogTdxValues = *temp1Ptr
dTdx -> trash
kappa_T = temp1 -> times(-1.)
kappa_T -> set, title='!4j!X!S!DT!R!U(e)!N', mnemonic='kappa_T', units='1/cm'

Normed_kappa_T = temp1 -> times(-rho_p)
Normed_kappa_T -> set, title='!4j!X!S!DT!R!U(e)!N', mnemonic='kappa_T', units='1/!4q!X!Dp!N'

dTidx = Ti -> dbyd('x')
temp3 = dTidx -> over(Ti)
temp3 -> Get, values=temp3Ptr
dLogTidxValues = *temp3Ptr
dTidx -> trash
kappa_Ti = temp3 -> times(-1.)
kappa_Ti -> set, title='!4j!X!S!DT!R!U(i)!N', mnemonic='kappa_Ti', units='1/cm'

Normed_kappa_Ti = temp3 -> times(-rho_p)
Normed_kappa_Ti -> set, title='!4j!X!S!DT!R!U(i)!N', mnemonic='kappa_Ti', units='1/!4q!X!Dp!N'

L_resistive = SQRT( L_mfp_e*SQRT(-mu*1836*R_o/(dLogNdxValues+dLogTdxValues)) )
L_resistive_Obj = n_e -> MAKECOPY(/novalues, /noerrorbars)
L_resistive_Obj -> set, values=PTR_NEW(L_resistive)
L_resistive_Obj -> set, title='L!S!D!9#!X!R!U(resistive)!N'
L_resistive_Obj -> set, mnemonic='L_resistive', units='cm'

Normed_L_resistive_Obj = L_resistive_Obj -> over(connection_Length)
Normed_L_resistive_Obj -> set,title='L!S!D!9#!X!R!U(resistive)!N/L!S!D!9#!X!R!U(connection)!N'
Normed_L_resistive_Obj -> set, mnemonic='Normed_L_resistive', units=''

rhoValues = 102*SQRT(mu*Tvalues)/B_o
normedValuesPtr = Normed_L_resistive_Obj -> GetValues()
normedValues = *normedValuesPtr
perpValues = SQRT(2.)*!PI*rhoValues/normedValues
Perp_Resistive_Obj = L_Resistive_Obj -> MakeCopy()
Perp_Resistive_Obj -> set, values=PTR_NEW(perpValues), title='L!S!D!9x!X!R!U(resistive)!N', units='cm', mnemonic='L_0'


d2Tdx = T -> d2byd('x')
temp2 = d2Tdx -> Over(T)
d2Tdx -> trash

temp2 -> Get, values=temp2Ptr
temp1 -> Get, values=temp1Ptr
dTvalues = *temp1Ptr
d2Tvalues= *temp2Ptr
temp1 -> trash
temp2 -> trash
temp3 -> trash

tau_perp = 1./( chi*(ABS(d2Tvalues) + dTvalues^2) )
tau_perp_Obj = n_e -> MakeCopy(/NoValues, /NoErrorBars)
tau_perp_Obj -> set, values=PTR_NEW(tau_perp)
tau_perp_Obj -> set, title='!4s!S!D!9x!X!R!U(E)!N', mnemonic='tau_perp', units='s'

normed_tau_perp_Obj = tau_perp_Obj -> times(nu_Energy)
normed_tau_perp_Obj -> set, title='!4s!S!D!9x!X!R!U(E)!N'+'!4m!X!DE!N', mnemonic='tau_perp_normed', units='s'

;
; profile of j_bootstrap following Sauter et al
;
f_t = SQRT(epsilon)

Z = Z_eff
f_31= f_t/( 1. + (1.-0.1*f_t)*SQRT(nuStar_e) + 0.5*(1.-f_t)*nuStar_e/Z )
L_31 = (1.+1.4/(Z+1.))*f_31 - 1.9/(Z+1.)*f_31^2 + 0.3/(Z+1.)*f_31^3 + 0.2/(Z+1.)*f_31^4

f_32_ee = f_t/( 1. + 0.26*(1-f_t)*SQRT(nustar_e) + 0.18*(1.-0.37*f_t)*nustar_e/SQRT(Z) )

L_32_ee = (0.05+0.62*Z)/( Z*(1+0.44*Z) )*(f_32_ee - f_32_ee^4) 				$
	+ 1./(1.+0.22*Z)*[ f_32_ee^2 - f_32_ee^4 - 1.2*(f_32_ee^3 - f_32_ee^4) ]	$
	+ 1.2/(1.+0.5*Z)*f_32_ee^4

f_32_ei = f_t/( 1. + (1.+0.6*f_t)*SQRT(nustar_e) + 0.85*(1.-0.37*f_t)*nustar_e*(1.+Z) )

L_32_ei = -(0.56+1.93*Z)/( Z*(1+0.44*Z) )*(f_32_ei - f_32_ei^4) 			$
	+ 4.95/(1.+2.48*Z)*[ f_32_ei^2 - f_32_ei^4 - 0.55*(f_32_ei^3 - f_32_ei^4) ]	$
	- 1.2/(1+0.5*Z)*f_32_ei^4
	
L_32 = L_32_ee + L_32_ei
	
f_34 = f_t/( 1. + (1.-0.1*f_t)*SQRT(nustar_e) + 0.85*(1.-0.37*f_t)*nustar_e*(1.+Z) )
L_34 = (1.+1.4/(Z+1.))*f_34 - 1.9/(Z+1.)*f_34^2 + 0.3/(Z+1.)*f_34^3 + 0.2/(Z+1.)*f_34^4

alpha_0 = - 1.17*(1-f_t)/( 1. - 0.22*f_t - 0.19*f_t^2 )

alpha 	= (   (alpha_0 + 0.25*(1.-f_t^2)*SQRT(nustar_i))/( 1 + 0.5*SQRT(nustar_i) )	$
	    - 0.315*nustar_i^2*f_t^6  )/( 1. + 0.15*nustar_i^2*f_t^6 )

p_e = 1.e6*nValues*Tvalues*1.6022e-19	; electron pressure in MPa

A_1 = dLogNdxValues + dLogTdxValues + (TiValues/Z)/Tvalues*(dLogNdxValues + dLogTidxValues)

A_2 = dLogTdxValues

A_2_i = dLogTidxValues

A_4 = alpha*(TiValues/(Z*Tvalues))*A_2_i

jProfile = - (p_e/BpBt)*(  1.e4/( B_o*SQRT(1.+BpBT*2) )  )*( L_31*A_1 + L_32*A_2 + L_34*A_4 )*1.e-4

jProfile_Obj = n_e -> MakeCopy(/NoValues, /NoErrorBars)
jProfile_Obj -> Set, values=PTR_NEW(jProfile), title='j!Dbootstrap!N', mnemonic='J_bootstrap', units='MA/m!U2!N'




output = {	NAME		:	"EdgeCollisionality"	, 	$
		nuStar_e	:	nuStar_e_Obj		, 	$
		nuStar_i	:	nuStar_i_Obj		, 	$
		l_mfp_e		:	l_mfp_e_Obj		,	$
		l_mfp_i		:	l_mfp_i_Obj		,	$
		L_resistive	:	L_resistive_Obj		,	$
		Normed_L_resistive:	Normed_L_resistive_Obj	,	$
		L_perp_resistive:	Perp_resistive_Obj	,	$
		nu_Energy	:	nu_Energy_Obj		,	$
		tau_ll		:	tau_ll_Obj		,	$
		tau_perp	:	tau_perp_Obj		,	$
		Normed_Tau_ll	:	Normed_tau_ll_Obj	,	$
		Normed_tau_perp	:	normed_tau_perp_Obj	,	$
		kappa_n		:	kappa_n			,	$
		kappa_T		:	kappa_T			,	$
		kappa_Ti	:	kappa_Ti		,	$
		Normed_Kappa_n	:	Normed_Kappa_n		,	$
		Normed_Kappa_T	:	Normed_Kappa_T		,	$
		Normed_Kappa_Ti	:	Normed_Kappa_Ti		,	$
		j_Bootstrap	:	jProfile_Obj			}

RETURN, output

END  ;****** Edge_Collisionality ******;