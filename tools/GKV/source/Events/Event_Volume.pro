;; *****************************************************************************************************************; ******************************************     Copyright Notice     *********************************************; *                                                                                                               *; *  This work was produced at the University of California, Lawrence Livermore National Laboratory (UC LLNL)     *; *  under contract no. W-7405-ENG-48 (Contract 48) between the U.S. Department of Energy (DOE) and The Regents   *; *  of the University of California (University) for the operation of UC LLNL. Copyright is reserved to the      *; *  University for purposes of controlled dissemination, commercialization through formal licensing, or other    *; *  disposition under terms of Contract 48; DOE policies, regulations and orders; and U.S. statutes. The rights  *; *  of the Federal Government are reserved under Contract 48 subject to the restrictions agreed upon by the DOE  * ; *  and University as allowed under DOE Acquisition Letter 97-1.                                                 *; *                                                                                                               *; *****************************************************************************************************************;; *****************************************************************************************************************; **********************************************     DISCLAIMER     ***********************************************; *                                                                                                               *; *  This work was prepared as an account of work sponsored by an agency of the United States Government.         *; *  Neither the United States Government nor the University of California nor any of their employees, makes      *; *  any warranty, express or implied, or assumes any liability or responsibility for the accuracy, completeness, *; *  or usefulness  *of any information, apparatus, product, or process disclosed, or represents that its use     *; *  would not infringe privately-owned rights.  Reference herein to any specific commercial products, process,   *; *  or service by trade name, trademark, manufacturer or otherwise does not necessarily constitute or imply its  *; *  endorsement, recommendation, or favoring by the United States Government or the University of California.    *; *  The views and opinions of authors expressed herein do not necessarily state or reflect those of the United   *; *  States Government or the University of California, and shall not be used for advertising or product          *; *  endorsement purposes.                                                                                        *; *                                                                                                               *; *****************************************************************************************************************;Function Curvature, data;; computes curvature of "data" surface;if(Query_Real(data) ne 1) then begin	PRINT, "ERROR in Curvature:  Argument is not real"	return, 0endifinfo=size(data)if (info(0) ne 2) then begin	PRINT, "ERROR in Curvature:  Expect 2-D array as argument"	return, 0endif;; Define kernal of 9-point difference operator;k1 = -1.0d/12.k2 = -1.0d/6.kernel = [ [k1,k2,k1], [k2, 1.0d, k2], [k1,k2,k1] ];; Convolve kernal with data to compute curvature;result=CONVOL(data, kernel, /EDGE_WRAP)return, resultendPro Find_Peak, array, map, max_volume, fraction = f, compute_volume=cv, p_skip=ps, debug=d;; This proceedure finds all maximia in ARRAY, and integrates to find volume of peaks surrounding each maxima.; The surrounding region is that region "downhill" from the maximia extending to the boarder with the neighboring maxima.; On return each element of the long integer array MAP contains the number of the maxima associated with that cell.; The 0th element of MAX_VOLUME contains the total volume of unmarked cells, ; while the ith element contain the volume associated with the ith maximia.;compute_volume=4if KEYWORD_SET(cv) then compute_volume=cvfraction=0.995if KEYWORD_set(f)  then fraction=fp_skip=100if KEYWORD_SET(ps) then p_skip=psinfo=size(array)if(info(0) ne 2) then begin	PRINT, "ERROR in Find_Peak: Array dimension is ", info(0)	returnendifnx=info(1)ny=info(2)n_total=nx*nymap=LONARR(nx, ny)Total_Volume=Total(array);; Make linear arrays of displacements;x_offsets=[-1,-1,-1, 0, 0, 1, 1, 1]y_offsets=[-1, 0, 1,-1, 1,-1, 0, 1];; find all maximia using only array operatiions for greater speed;for ii=0,7 do map = map + ( array gt shift(array ,x_offsets(ii), y_offsets(ii)) )map = LONG(map gt 7) ;; MAP elements corresponding to maxima of ARRAY are now marked with 1. ; All other elements of MAP are zero.;n_maxima=Total(map)Result=HISTOGRAM(map, binsize=1, min=1, max=1, REVERSE_INDICES=R)n_max = LONG(R(1) - R(0))if( n_max ne n_maxima) then begin	PRINT, "ERROR in Find_Peak: REVERSE_INDICES array does not have expected structure"	returnendiffor max_number=1L, n_max do begin	index = R[max_number+1]	map[index] = max_numberendfor;; MAP elements corresponding to maxima of ARRAY are now marked with ; a unique Max_Number -- an integer between 1 and n_max.; All other elements are still zero;; Compute number of marked cells, n_marked;old_marked=0Ln_marked = TOTAL(map gt 0)iteration_number = 0L;if KEYWORD_SET(d) then PRINT, "Find_peak finished marking Maxima. N_max = ", n_maxstart_of_loop: iteration_number = iteration_number+1L;; Find all neighbors of marked cells,; mark with Max_Number of local maxima ; (which is stored in non-zero elements of MAP);;	if KEYWORD_SET(d) then PRINT, format='(A43,$)', "Compleated iteration over neighbors number "for ii=0,7 do begin	map = map + (map eq 0)*( array lt shift(array,x_offsets[ii],y_offsets[ii]) )*shift(map,x_offsets[ii],y_offsets[ii]);	if KEYWORD_SET(d) then PRINT, format='(T1, i2)', iiendforif(iteration_number gt 1000) then begin	PRINT, format='(A61, I3)', "ERROR in Find_Peaks: too many iterations, iteration_number = ", iteration_number	returnendifold_marked = n_markedn_marked = TOTAL(map gt 0)if KEYWORD_SET(d) then begin	PRINT, format='(A27, I8)', "Find_Peak iteration number ", iteration_number	PRINT, format='(A27, I8)', "       Total cells marked =", n_marked	PRINT, format='(A27, I8)', "Total # of cells in array =", n_totalendifif(n_marked gt fraction*n_total) then goto, finished_loop  ; Stop iteration when "fraction" of cells are marked (defaults to 99.5%)if (n_marked gt old_marked ) then goto, start_of_loop      ; Stop iteration when algorithm fails to find any more cells to markfinished_loop: $if KEYWORD_SET(d) then PRINT, "Computing volumes of each marked region (this can take a while...)'if KEYWORD_SET(d) then PRINT, format='(A20, I6 )', "Number of regions =", n_max;; Compute Volumes of each marked region;Max_Volume=dblarr(n_max+1);if KEYWORD_SET(d) then PRINT, format='(A16, i1)', "compute_volume= ", compute_volumeif (compute_volume eq 1) then begin  ; Use logical array operations 	for max_number=0L, n_max do begin		max_volume(max_number) = TOTAL( (map eq max_number)*array )		if ( KEYWORD_SET(d) and (P_skip*(max_number/P_skip)  eq max_number) ) then PRINT, format='(A12, I6)', "max_number =", max_number	endforendif if (compute_volume eq 2) then begin  ; Use HISTOGRAM method	result= HISTOGRAM(map, binsize=1, min=0, max=n_max, REVERSE_INDICES=R)	for max_number=0L, n_max do begin		if( R[max_number] lt R[max_number+1]-1L) then begin			index=R[ R[max_number]:R[max_number+1]-1L ]			max_volume(max_number)=total(array[index])		endif		if ( KEYWORD_SET(d) and (P_skip*(max_number/P_skip)  eq max_number) ) then PRINT, format='(A12, I6)', "max_number =", max_number	endforendif if( compute_volume eq 3) then begin  ; Use Where method	for max_number=0L, n_max do begin		index=WHERE(map eq max_number, count)		if (count ne 0) then max_volume(max_number)=total(array[index])		if ( KEYWORD_SET(d) and (P_skip*(max_number/P_skip)  eq max_number) ) then PRINT, format='(A12, I6)', "max_number =", max_number 	endforendifif (compute_volume eq 4) then begin  ; Use pure brute force...(which, oddly enough, seems to give rather good performance for VERY lorge arrays!	for i=0L, nx-1L do begin		for j=0L, ny-1 do begin			max_volume(map(i,j))=max_volume(map(i,j)) + array(i,j)		endfor		if ( KEYWORD_SET(d) and (P_skip*(i/P_skip)  eq i) )  then PRINT, format='(A4, I6)', "nx =", i	endforendif present_volume=TOTAL(max_volume)if KEYWORD_SET(d) then PRINT, format='(/, A16, e10.2)',    "Volume marked = ", present_volumeif KEYWORD_SET(d) then PRINT, format='(A16, e10.2)',    " Total volume = ", total_volumeif KEYWORD_SET(d) then PRINT, format='(A16, F6.2, A1)', "      Percent =", 100*present_volume/Total_volume, "%"returnend			FUNCTION FilterDown,	data, C=curvature, max_its=its, err_tol=errtol, $			  	rms_err=rms_errs, max_err=max_errs, Debug=d;;  This function filters the input array, DATA, removing local peaks.  ;  The filtered array is always LESS THAN OR EQUAL TO the input.;  max_its=1if KEYWORD_SET(its) then max_its=itsC=0.if KEYWORD_SET(curvature) then C=curvatureerr_tol=1.e-4if KEYWORD_SET(errtol) then err_tol=errtolinfo=SIZE(data)if( info(0) ne 2 ) then begin	PRINT, "ERROR in FIlterDown:  Input array dimension is ", info(0)	return, 0endifnx=info(1)nt=info(2)result=dataresult1=FLTARR(nx, nt);; Make kernel for 9-point averaging stencil;k1=1.0d/12.k2=1.0d/6.kernel = [ [k1,k2,k1], [k2, 0,k2], [k1,k2,k1] ];rms_errs=1.0				;  Initialize rms_errs variable...max_errs=1.0				;  Initialize max_errs;for n=1L, max_its do begin	result1=0.*result1						;  reset temporary array to zero	result1=CONVOL(result, kernel, /EDGE_WRAP) + C	;  average neighboring values & add curvature	result1= result1 < result					;  Constrain filter to only reduce RESULT	err=SQRT(TOTAL((Result-Result1)^2)/Total(Result^2))	max_err=MAX(result-result1)*nx*nt/TOTAL(result)	rms_errs=[rms_errs, err]					;  concatenate err into "errs" array	max_errs=[max_errs, max_err]				;  concatenateMax_err into "max_errs" array	if KEYWORD_SET(d) then PRINT, "Iteration # ", n, "     rms_err = ", err, "     max_err = ", max_err	result=result1	if(err lt err_tol) then goto, doneendfordone: $return, resultend	function Boundary_Time, 	data, Front=ff, Back=bb,	$ 					Constant=CC, Constant_Gradient=CG, Zero=ZZ, Double=ddinfo=size(data)if (info(0) ne 2) then begin	PRINT, "ERROR in Boundary_Time:  Expected 2-D array as as input"	return, 0endifif( Query_real(data) ne 1) then begin	PRINT, "ERROR in Boundary_Time:  Input array is not real"	return, 0endifnx=info(1)nt=info(2)if KEYWORD_SET(CC) then goto, Constant_valuesif KEYWORD_SET(zz) then goto, zero_curvature;; Default--"Constant Gradient: boundary:;; Values in boundary layer are extrapolate from interior points ; assuming constant �/�t;border = 2.0*data(*,0) - data(*,1)if KEYWORD_SET(BB) then border = 2.0*data(*,nt-1) - data(*,nt-2)return, border; ; "Constant" boundary:;; Values in boundary layer are equal to those in the first column of data ;Constant_values:	$border = data(*,0)if KEYWORD_SET(BB) then border = data(*,nt-1)return, border;; "Zero" boundary:;; "zero-curvature" extrapolation of "data" into boundary columns; at the front (default) or back of data. ; NOTE:  it is NOT POSSIBLE to remove the "even-odd" component of; the "curvature" on the boundary with this extrapolation for the ; kernal used in FILTER_DOWN, as the 1, 2, 1 weighting annihilates; any even/odd component in the boundary layer.;zero_curvature:	$zeros=dblarr(nx)temp=dblarr(nx,4)temp(*,  0) = zerostemp(*,1:3) = data(*,0:2)if KEYWORD_SET(bb) then begin	temp(*,0:2) = data(*,nt-3:nt-1)	temp(*,  2) = zerosendiftemp1=Curvature(temp)source=temp1(*,1)if KEYWORD_SET(bb) then source=temp1(*,2)c_k = fft(source, -1, Double=dd)k = 2.0d*!Dpi*dindgen(nx)/(nx)x_k = 6.0d*c_k/(1.0d + COS(k))x_k(nx/2) = 0.0border=fft(x_k, 1, Double=dd);; check result;temp2 = dblarr(nx,4)temp2(*,0) = bordertemp2(*,1:3) = data(*,0:2)if KEYWORD_SET(bb) then begin	temp2(*,0:2) = data(*,nt-3:nt-1)	temp2(*,  3) = borderendiftemp3=curvature(temp2)test=temp3(*,1)if KEYWORD_SET(bb) then test=temp3(*,2)print, max(test)return, borderendFunction Scale, array, result, map, max_vol, Curvature,	$			rms_err=rms_err, max_err=max_err,		$						err_tol=errtol, max_its=max_its, Debug=d;;  Applies filter to input ARRAY, to separate events;  at scale defined by Curvature.  These events are;  stored in RESULT, while the residual (ARRAY - RESULT);  is returned.;;  At input,;;	ARRAY		contains 2-D array of data to be filtered.;;	CURVATURE	contains the curvature of the surface to be fit;			to ARRAY;;	ERR_TOL	contains the maximum RMS error to be allowed;			(provided no more than MAX_ITS iterations are ;			 required to achieve this error);;	MAXITS	contains the maximum number of iterations;			allowed in attempting to achieve an RMS error;			less than ERR_TOL;;  On return, ;;	RESULT	a 2-D array (same dimensions as ARRAY) containing;			features at scale defined by CURVATURE;;	MAP		a 2-D integer array (same dimensions as ARRAY) maping;			location of "events" (essentially local maxima) found;			by FIND_PEAK within the RESULT array.;;	MAX_VOL	a 1-D array containing volume of each local maximia found;			in RESULT.  MAX_VOL(0) contains the volume not associated;			with any maxima;	;	RMS_ERR	will contain decreasing sequence of RMS;			deviations across grid for each iteration;			required to filter ARRAY to specified tolerance;;	MAX_ERR	will contain decreasing sequence of maximum;			deviations across grid for each iteration;			required to filter ARRAY to specified tolerance;;info=size(array)if(info(0) ne 2) then begin	PRINT, "ERROR in Scale: Array dimension is ", info(0)	return, 0endiferr_tol=1.e-4if KEYWORD_SET(errtol) then err_tol=errtolnx=info(1)ny=info(2)nt=nyn_total=nx*ny;; Add boundary layer (in "y", or time direction).; Default is constant �/�t.;left_border =Boundary_Time(array, /Front)right_border=Boundary_time(array, /Back )temp=dblarr(nx, ny+2)temp(*,0)=left_bordertemp(*,nt+1)=right_bordertemp(*,1:nt)=array;; compute filtered array; result1=FilterDown(temp, C=Curvature, max_its=max_its, err_tol=err_tol, $			  rms_err=rms_err, max_err=max_err, Debug=d)temp1 = temp - result1diff=temp1(*,1:nt)result=result1(*,1:nt);;  Count events of largest scale;result1(*,   0)=0.0d					; Set boundary layers to zero to prevent "events" result1(*,nt+1)=0.0d					; from propagating across periodic boundary in time.Find_Peak, result1, map, max_vol, Debug=dreturn, diffendFUNCTION Event_Volume, ARRAY, SPAN=s, nScales=n_scales, C0=cc, FEATURES=features, $		   		  MAPS=maps, RMS_ERRORS=rms_errors, MAX_ERRORS=max_errors, $	        		  ERR_TOL=errtol, MAX_ITS=maxits, DEBUG=d;;  This Function returns a structure containing:;;	MAX_VOL	An array containing as its elements the volume of each "event" found in ARRAY.  ;			Note, however, tham EVENT_Volume(0) contains the total volume that was not ;			"marked" by Find_Peak. If this volulme issignificant, then validity of resulting ;			PDF is in doubt.;;	N_Events	N_Events is an array containing the "address" within MAX_VOL;			of the "Events" found at each scale.  This array used by PDFs ;			to locate the Event_Volumes corresponding to each scale in the Max_Vol array;			and construct the partial probability distribution functions at each scale.  ;			The final element of N_Events contains the total number of iterations over scales.;			;  ;; KEYWORDS:;;	nScales		On input, 'nScales' determines the maximum number of space/time ;			'scales' in the filter.  Defaults to 7. (Optional);;	SPAN		If set on input, SPAN determines the rate at which FILTER_DOWN;			decreases the spatial scale on each iteration.;			Default is to determine SPAN from nScales. (Optional);;	CO		Determines the initail "curvature" of the filtered surface;			in units of (1/L_x^2).;			(defaults to 0.1);;	FEATURES	On return, FEATURES is a 3-D array containing a sequence of ;			(2-D) collections of features from ARRAY at each scale;			in the iteration over scales.;;	MAPS		On return, MAPS is a 3-D array containing a sequence of (2-D) maps of;			"events" extrated from the "features" at the corresponding level.;;	RMS_ERRORS	On return, RMS_ERRORS is a 1-D array containing a record of ;			the relative RMS error in the determination of the filtered  ;			images of ARRAY at each iteration.  This information is ;			packed with the RMS_ERRORS(0) equal to the number of iterations ;			required for constructing the first filtered image, followed by;			RMS_ERRORS(0) rms error values.  RMS_ERRORS( RMS_ERRORS(0) ) then;			contains the number of iterations required for ocnstructing the ;			second filtered image, etc.;;	MAX_ERRORS	On return, MAX_ERRORS is a 1-D array analogous to RMS_ERRORS.;			MAX_ERRORS contains the maximum error at each iteration.;;	ERR_TOL		On input, ERR_TOL contains the maximum tollerable relative ;			error in any filtered image of ARRAY (defaults to 1.0 e-4);;	MAX_ITS		On input, MAX_ITS contains the maximum number of iterations ;			(in units of NX*NT) in producing any filtered image of ARRAY.  ;			FILTER_DOWN will only perform MAX_TIS iterations even if the ;			RMS error is not reduced below ERR_TOL. (Defaults to 2);;  Check validity of input ARRAY;IF( QUERY_REAL(ARRAY) ne 1) then begin	Print, "ERROR in Event_Volume:  ARRAY is not real"	return, 0endifinfo=SIZE(array)if(info(0) ne 2) then begin	PRINT, "ERROR in Event_Volume: Array dimension is ", info(0)	return, 0endifnx=info(1)nt=info(2)n_total=nx*ntC0=0.1if( KEYWORD_SET(cc) ) then C0=cccurve=C0*MAX(array)/(nx^2);; Compute default SPAN;left_bdy = Boundary_Time(array, /Front)	;  Get boundary cells to load at  right_bdy= Boundary_Time(array,  /Back)	;  front and back in  time.temp=fltarr(nx, nt+2)temp(*,0)    = left_bdy				;  Copy front boundary cells into TEMP			temp(*,1:nt) = array				;  Copy ARRAY into TEMPtemp(*,nt+1) = right_bdy			;  Copy back boundary cells into TEMPmax_curve = Max(Curvature(temp))		;  Compute maximum value of curvature in input ARRAYnScales = 7.0IF(N_ELEMENTS(n_Scales) EQ 1) THEN nScales=n_Scales > 1.span= exp(aLog(max_curve/Curve)/nScales)	;  Choose SPAN for 6 iterations to smallest scaleif( KEYWORD_SET(s) ) then span=sprint, "Event_Volume:  Minimum Curvature = ", Curveprint, "Event_Volume:  Maximum Curvature = ", Max_curveprint, "Event_Volume:  Span = ", spanerr_tol=1.0e-3if( KEYWORD_SET(errtol) ) then err_tol=errtolmax_its=2*(nx+nt)if( KEYWORD_SET(maxits) ) then max_its=maxits;;  Initialize variables for iteration over scales;total_vol=TOTAL(array);;  Filter input ARRAY to find features at largest scale;resid=Scale(array, feature, map, max_vol, curve, 	$		rms_err=rms_err, max_err=max_err,		$		err_tol=err_tol, max_its=max_its, debug=d)features=feature				; initialize FEATURES arraymaps=map						; initialize MAPS arrayn_events=N_ELEMENTS(MAX_VOL)-1		; initialize N_EVENTS arrayrms_errors=N_ELEMENTS(rms_err)		; initialize RMS_ERRORS arrayrms_errors=[rms_errors, rms_err]	; set first element of RMS_ERRORS arraymax_errors=N_ELEMENTS(max_err)		; initialize MAX_ERRORS arraymax_errors=[rms_errors, max_err]rel_resid_vol=TOTAL(resid)/total_volit=0;;  Begin iteration over scales;WHILE (rel_resid_vol gt err_tol) do begin	it=it+1	curve= span*curve				;  Increase curvature of filtered surface	max_its = ( max_its/SQRT(span) ) > 1	;  Fewer iterations required for smaller features	;  	resid1 =Scale(resid, feature, map, max_vol1, curve, 	$				rms_err=rms_err, max_err=max_err,		$				max_its=max_its, err_tol=err_tol, debug=d	)	;	features=[[[features]], [[feature]]]	;  Save features at current scale	maps=[[[maps]], [[map]]]			;  Save map of events at current scale	max_vol(0)=max_vol(0)+max_vol1(0)	;  Pack unmarked volulme into max_vol(0)	max_vol=[max_vol, max_vol1(1:*)]	;  Add volumes of events at current scale to MAX_VOL								;  Concatenate # if events at current scale into	n_events=[n_events, N_ELEMENTS(MAX_VOL1)-1 ]	; N_EVENTS array.	;	;  add error sequences from call to SCALE to RMS_ERRORS and MAX_ERRORS	;	rms_errors=[rms_errors, N_ELEMENTS(rms_err), rms_err]		max_errors=[max_errors, N_ELEMENTS(max_err), max_err]		;	;  prepare for next iteration	;	resid=resid1	rel_resid_vol=TOTAL(resid)/total_vol		print, "Iternation # = ",it, "     curve = ", curve, "     rel_resid_vol = ", rel_resid_vol	if(it eq 10) then goto, done		;  10 interations should be plenty!ENDWHILE;; Finished with iteration over event scales.  Return structure containing information needed to ; construct and plot probability distribution function (PDF).;done: n_events = [n_events, it]			; Concatenate total number of iterations into final								; element of n_events array.Event_Vol = { Max_Vol:Max_Vol, n_events:n_events }								return, Event_volend