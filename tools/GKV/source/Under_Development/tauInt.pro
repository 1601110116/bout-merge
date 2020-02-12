FUNCTION GKVs2D::AVG_t, lw=ll, tauCorrINT=tci, std=stdd, skewness=skwn, kurtosis=kurt
;
; Compute the average of 'self' over final (2nd) independent
; variable, and return the result as a GKVs1D object.
;
;  Written by W.M. Nevins
;	6/19/01
;
temp = self -> MakeCopy()
temp -> Restrict
irange = temp.grid1.irange
imin = irange[0]
imax = irange[1]
nx = imax - imin + 1
firstSlice = temp -> slice(axis=1, index=imin)
IF( (Query_Integer(ll) + Query_Real(ll)) EQ 0) THEN BEGIN
	firstCorrs = firstSlice -> xcorr()
	fullWidth = firstCorrs -> FullWidth()
	tRange = self.grid2.range
	tauHalf = fullWidth/2.
	tTot = tRange[1] - tRange[0]
	lw = SQRT(tauHalf*tTot)
ENDIF ELSE BEGIN
	lw = ll
ENDELSE
;
; set up arrays for output values
;
avg = FLTARR(nx)
avgPM = FLTARR(nx)
tauCorrInt = FLTARR(nx)
std = FLTARR(nx)
tauCorrInt = FLTARR(nx)
skewness = FLTARR(nx)
kurtosis = FLTARR(nx)
;
; Get data for first time slice
;
stats = firstSlice -> STATS(lw=lw)
;
; Store it in appropriate array
;
avg[0] = stats.avg
avgPM[0] = stats.avgPM
tauCorrInt[0] = stats.tauCorrInt
std[0] = stats.std
tauCorrInt[0] = stats.tauCorrInt
skewness[0] = stats.skewness
kurtosis[0] = stats.kurtosis
;
; Start loop over first independent variable
;
FOR i=imin+1, imax DO BEGIN
	thisSlice = temp -> Slice(axis=1, index=i)
	stats = thisSlice -> stats(lw=lw)
	thisSlice -> trash
	avg[i] = stats.avg
	avgPM[i] = stats.avgPM
	tauCorrInt[i] = stats.tauCorrInt
	std[i] = stats.std
	tauCorrInt[i] = stats.tauCorrInt
	skewness[i] = stats.skewness
	kurtosis[i] = stats.kurtosis
ENDFOR
;
; Make output structures
;
avgStr = {GKVs1D}
FOR i=0, N_TAGS({GKVsd}) DO avgStr.(i) = temp.(i)
avgStr.title = '<' + temp.title + '>'
avgStr.mnemonic = temp.mnemonic + '_avg'
indices = ['*']
avgStr.indices = PTR_NEW(indices)
avgStr.values = PTR_NEW(avg)
avgStr.ErrorBars = PTR_NEW(avgPM)
vMin = GKVsd_MIN(avg, MAX=vMax)
avgStr.Vrange = [vMin, vMax]
grid = GKVsd_GridCopy(temp.grid1)
avgStr.grid1 = grid
avgObj = OBJ_NEW('GKVs1D', avgStr)

temp -> trash
RETURN, avgObj

END ; ****** FUNCTION GKVs2D::AVG_t ****** ;

FUNCTION GKVs1D::TauInt, MaxSteps=maxSteps, iLagMax=iLagMaxIn, LagMax=LagMaxIn
;
;
; Find length of (packed) correlation function
;
selfCorrs = self -> xcorr(/pack, /norm, lw=0, /NoAvg)
CorrFcn = FLOAT(*selfCorrs.values)
nCorrs = N_ELEMENTS(CorrFcn)
iLagMin=0L
iLagMax =nCorrs/2L
;
; Compute dTau
;
selfGridValues = *self.grid1.values
dTau = selfGridValues[1] - selfGridValues[0]
;
; Compute maximum lag window
;
IF(N_ELEMENTS(ilagMaxIn) EQ 1) THEN iLagmax = LONG( iLagMaxIn < iLagMax )
IF(N_ELEMENTS(LagMaxIn ) EQ 1) THEN iLagMax = LONG( (LagMaxIn/dTau + 1L)  < iLagMax)
;
; Compute the 'skip' for loop over values of the lag window
;
iskip = 1L
IF(N_ELEMENTS(maxSteps) EQ 1) THEN BEGIN
	iskip = iLagMax/maxSteps > 1L
ENDIF
;
; Compute tauInt, tauMax vs. length of lag window
;
nLags = (iLagMax)/iSkip + 1L
tauInt = FLTARR(nLags)
tauMax = FLTARR(nLags)
tauOne = FLTARR(nLags)
lagValues = FLTARR(nLags)
tauInt[0] = 0.
tauMax[0] = 0.
tauOne[0] = 0.
index = -1L
FOR ilw=1L, ilagMax-1L, iskip DO BEGIN
	index = index+1L
	temp = (LagWindow_Fcn(nCorrs, ilw)*CorrFcn)[ncorrs/2-ilw:ncorrs/2+ilw]
	temp = SHIFT(temp, -ilw)
	int = TOTAL(temp, /Cumulative)
	tauInt[index] = int[2*ilw]*dTau
	tauMax[index] = 2.*MAX(int[0:ilw], iMax)*dTau
	tauOne[index] = (iMax+.5)*dTau
	lagValues[index] = dTau*ilw
ENDFOR
;
; Make output objects
;
tauIntStr = {GKVs1D}
tauMaxStr = {GKVs1D}
tauOneStr = {GKVs1D}
FOR i=0, N_TAGS({GKVsd})-1 DO tauIntStr.(i) = self.(i)
FOR i=0, N_TAGS({GKVsd})-1 DO tauMaxStr.(i) = self.(i)
FOR i=0, N_TAGS({GKVsd})-1 DO tauOneStr.(i) = self.(i)

tauIntStr.title = '!4s!X!DINT!N {' + self.title + '}'
tauMaxStr.title = '!4s!X!DMAX!N {' + self.title + '}'
tauOneStr.title = '!4s!X!D1!N {' + self.title + '}'

tauIntStr.mnemonic = 'tau_int'
tauMaxStr.mnemonic = 'tau_max'
tauOneStr.mnemonic = 'tau_1'

indices = ['*']
tauIntStr.indices = PTR_NEW(indices)
tauMaxStr.indices = PTR_NEW(indices)
tauOneStr.indices = PTR_NEW(indices)

tauIntStr.values = PTR_NEW(tauInt)
tauMaxStr.values = PTR_NEW(tauMax)
tauOneStr.values = PTR_NEW(tauOne)

vmax = MAX(tauInt)
tauIntStr.vrange = [0., vmax]
vMax = MAX(tauMax)
tauMaxStr.vrange = [0., vmax]
vMax = MAX(tauOne)
tauOneStr.vrange = [0., vmax]

tauIntStr.units = self.Grid1.units
tauMaxStr.units = self.Grid1.units
tauOneStr.units = self.Grid1.units
;
; Make output grid structures
;
tauIntGrid = {Grid}
tauIntGrid.values = PTR_NEW(lagValues)
tauIntGrid.range = [0., dTau*(iLagMax-1L)]
tauIntGrid.irange = [0, nLags-1L]
tauIntGrid.units = self.Grid1.units
tauIntGrid.title = '!4s!X!Dlag!N'
tauIntGrid.mnemonic = 'tau_lag'

tauMaxGrid = tauIntGrid
lagValuesMax = lagValues
tauMaxGrid.values = PTR_NEW(lagValuesMax)

tauOneGrid = tauIntGrid
lagValuesOne = lagValues
tauOneGrid.Values = PTR_NEW(lagValuesOne)
;
; Load grid structures
;
tauIntStr.grid1 = tauIntGrid
taumaxStr.grid1 = tauMaxGrid
tauOneStr.grid1 = tauOneGrid
;
; Create output objects
;
tauIntObj = OBJ_NEW('GKVs1D', tauIntStr)
tauMaxObj = OBJ_NEW('GKVs1D', tauMaxStr)
tauOneObj = OBJ_NEW('GKVs1D', tauOneStr)
;
; Create output structure
;
output = { Name:'tauInt', tauInt:tauIntObj, tauMax:tauMaxObj, tau1:tauOneObj }
RETURN, output
END  ; ****** GKVs1D::TauInt ****** ;