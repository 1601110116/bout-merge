;+ ; NAME: ;     MTANH_FUNCT;; PURPOSE: ;     Calculates a modified hyperbolic tangent function for fitting the density;     across the separatrix derived from the Thomson data;; CALLING SEQUENCE: ;     MTANH_FUNCT,X,A,F [,PDER];; INPUT PARAMETERS: ;     X : Elevation of point where ;     A : Function coefficients for the fit function;; OPTIONAL INPUT PARAMETERS: ;     PDER : Specify this variable if the partial derivatives of the fit ;            function are desired.;; KEYWORDS: ;     NONE;; OUTPUTS: ;     F :  Function values at the X locations;     PDER:  Partial derivative values at the X locations;; COMMON BLOCKS: ;     NONE;; SIDE EFFECTS: ;     NONE;; RESTRICTIONS:;     NONE;; PROCEDURE: ;     Calulate F = A(0) - A(1)*TANH(Z);        where Z=(X-A(2))/A(3);;     Calculate the partial derivatives, PDER, if specified on input;; CODE TYPE: modeling, analysis;; CODE SUBJECT:  edge, rf, transport, equilibrium;; EASE OF USE: can be used with existing documentation;; OPERATING SYSTEMS:  UNIX of all flavors;; EXTERNAL CALLS:  NONE;; RESPONSIBLE PERSON: Ray Jong;	; DATE OF LAST MODIFICATION:  02/17/99;; MODIFICATION HISTORY:;     Created by Gary D. Porter, LLNL;;       MODIFIED TO ADD LINEAR TERM 14 MAR 97 GDP;       Improvements suggested by Wolfgang Suttrop (12 May 97);       Name changed from gslts_funct to tanh_funct 17 FEB 99 RAJ;;-	PRO    MTANH_FUNCT,X,A,F,PDER    old_quiet=!quiet    !quiet=1    z=(x-a(2))/a(3)        PTAN = TANH(Z)        qtan=-1.*ptan        rtan=z*exp(-z)/(exp(z)+exp(-z))    f=a(0) - a(1)*PTAN -a(1)*a(4)*rtan    epsilon=1.e-9    df=-a(1)*((1-tanh(z)*tanh(z))+a(4)*rtan*(1-z-z*tanh(z))/z)        IF N_PARAMS(0) LE 3 THEN begin          !quiet=old_quiet          RETURN ;NEED PARTIAL?        endif         pder = [  [REPLICATE(1.,N_ELEMENTS(X))], $                  [qtan-a(4)*rtan], $                  [-df/(a(3)+epsilon)], $                  [-z*df/(a(3)+epsilon)], $                  [-a(1)*rtan]  ]        RETURNEND