PRO PERIOD_v5, other_args

; call like bash -i -c 'idl -e period_v3 -args' FILE_PATH, GTI_FILE_PATH, OUTPUT_FILE_PATH, obsID, OFAC, HIFAC, energy_range, obs_no

;OFAC,HIFAC,FILE_PATH,PX,PY,NOUT,JMAX,PROB
;
; Computes a Lomb-Scargle periodogram (power spectrum)
;
; I/O Variables:
;
; X (input) -- times
; Y (input) -- fluxes
; OFAC (input) -- Lomb-Scargle oversampling factor (try 1-8)
; HIFAC (input) -- Lomb-Scargle high frequency factor (try 1)
; PX (output) -- frequencies
; PY (output) -- Power spectrum
; NOUT (output) -- Number of frequencies/periods
; JMAX (output) -- index of best frequenc/(period 
; PROB (output) -- false alarm probability of strongest frequency/period
;                  based on white noise assumption.
;

;hardcoded or passed parameters section
args = command_line_args()
FILE_PATH = args[0] 
GTI_FILE_PATH = args[1] 
OUTPUT_FILE_PATH = args[2]
obsID = args[3]
OFAC = args[4]
HIFAC = args[5]
energy_range = args[6]
obs_no = args[7]



time_dat = mrdfits(FILE_PATH, 1, hdr, columns=['TIME'])
flux_dat = mrdfits(FILE_PATH, 1, hdr, columns=['RATE'])
frmtm = sxpar(hdr, 'TIMEDEL')

print, "frametime passed ", frmtm


;if frmtm eq 73 then frmtm = .0734
 ;     if frmtm eq 199 then frmtm = .1991
 ;     if frmtm eq 50 then exit
 ;     if frmtm eq 0 then frmtm = 2.6
	;  if frmtm eq 48 then exit
;print, "this time frame time: ", frmtm 

 
temp_X = time_dat.time
temp_Y = flux_dat.rate

;sort out GTIs 
num_lines = file_lines(GTI_FILE_PATH)
en = ''
openr, snow, GTI_FILE_PATH, /get_lun
  for i = 0, num_lines-1 do begin
      readf, snow, en
      bts = strsplit(en, /extract)
	    temp = where(temp_X gt bts[0] and temp_X lt bts[1])
	    if i eq 0 then (cumulative_X = temp) else (cumulative_X = [cumulative_X,temp])
      if i eq 0 then (cumulative_Y = temp) else (cumulative_Y = [cumulative_Y,temp])
  endfor

      
X = temp_X(cumulative_X)
Y = temp_Y(cumulative_Y)

;added to get rid of bright sources that take ages to process
if (mean(Y) gt 0.5) then begin
  exit
endif

rebinner,X,Y,frmtm,0.01,newX,newY
X=newX
Y=newY

;put into test nfft code --- remove 

;openw, rebinnerX, '/home/UTU/andmas/Desktop/rebin_resultsX.txt', /get_lun
;openw, rebinnerY, '/home/UTU/andmas/Desktop/rebin_resultsY.txt', /get_lun
;printf, rebinnerX, X, format = '(f0.12)' 
;printf, rebinnerY, Y, format = '(f0.12)'
;free_lun, rebinnerX
;free_lun, rebinnerY
;close, /all


N=N_ELEMENTS(X)
NMAX=500000L
TWOPID=2.0*!PI
HALFD=DOUBLE(0) & TWOD=DOUBLE(0)
HALFD=0.5 & TWOD=2.0
np=round(0.5*ofac*hifac*n)
PX=DBLARR(NP)
PY=DBLARR(NP)
WPR=FLTARR(NMAX) & WPI=FLTARR(NMAX) & WR=FLTARR(NMAX) & WI=FLTARR(NMAX)
NOUT=0.5*OFAC*HIFAC*N
IF(NOUT GT NP) THEN PRINT,'OUTPUT ARRAY TOO SHORT'
VAR=STDEV(Y,AVE)^2
XMIN=MIN(X) & XMAX=MAX(X)
XDIF=XMAX-XMIN
XAVE=0.5*(XMAX+XMIN)
PYMAX=0.0
PNOW=1.0/(XDIF*OFAC)
;
ARG=TWOPID*((X-XAVE)*PNOW)
WPR=-TWOD*SIN(HALFD*ARG)^2
WPI=SIN(ARG)
WR=COS(ARG)
WI=WPI
;
old=0
i=long(0)
WHILE I lt NOUT DO BEGIN
  PX(I)=PNOW
  SUMSH=0.0
  SUMC=0.0
  SUMSH=TOTAL(WI*WR)
  SUMC=TOTAL((WR-WI)*(WR+WI))
  WTAU=0.5*ATAN(2.0*SUMSH/SUMC)
  SWTAU=SIN(WTAU)
  CWTAU=COS(WTAU)
  SUMS=0.0
  SUMC=0.0
  SUMSY=0.0
  SUMCY=0.0
  SS=WI*CWTAU-WR*SWTAU
  CC=WR*CWTAU+WI*SWTAU
  SUMS=TOTAL(SS^2)
  SUMC=TOTAL(CC^2)
  YY=Y-AVE
  SUMSY=TOTAL(YY*SS)
  SUMCY=TOTAL(YY*CC)
  WTEMP=WR
  WR=(WR*WPR-WI*WPI)+WR
  WI=(WI*WPR+WTEMP*WPI)+WI
  PY(I)=0.5*(SUMCY^2/SUMC+SUMSY^2/SUMS)/VAR
  IF (PY(I) GE PYMAX) THEN BEGIN
    if(I gt 10) then begin
      PYMAX=PY(I)
      JMAX=I
    endif
  ENDIF
  PNOW=PNOW+1.0/(OFAC*XDIF)
;  
  i=i+1
 
ENDwhile
EXPY=EXP(-PYMAX)
EFFM=2.0*NOUT/OFAC
PROB=EFFM*EXPY
IF(PROB GT 0.01) THEN PROB=1.0-(1.0-EXPY)^EFFM
II=WHERE(PX GT 0.0)
PX=PX(II)
PY=PY(II)


;insert code to write results to file
fname_freq = OUTPUT_FILE_PATH + obsID + "_" + energy_range + '_obs_no_' + obs_no +  '_freq.dat'
fname_pow = OUTPUT_FILE_PATH + obsID  + "_" + energy_range + '_obs_no_' + obs_no + '_power.dat'
fname_rest = OUTPUT_FILE_PATH + obsID  + "_" + energy_range + '_obs_no_' + obs_no + '_restresults.dat'
fname_tbin = OUTPUT_FILE_PATH + obsID + "_" + energy_range + '_obs_no_' + obs_no +  '_tbin.dat'
fname_fbin = OUTPUT_FILE_PATH + obsID + "_" + energy_range + '_obs_no_' + obs_no +  '_fbin.dat'
openw, freq, fname_freq, /get_lun
openw, pow, fname_pow, /get_lun
openw, rest, fname_rest, /get_lun
openw, tbin, fname_tbin, /get_lun
openw, fbin, fname_fbin, /get_lun
printf, freq, PX, format = '(f0.18)' 
printf, pow, PY, format = '(f0.18)' 
printf, rest, NOUT, format = '(f0.18)'
printf, rest, JMAX, format = '(f0.18)' 
printf, rest, PROB, format = '(f0.18)'
printf, tbin, X, format = '(f0.18)'
printf, fbin, Y, format = '(f0.18)' 
print, "reached end of writes" 
free_lun, freq
free_lun, pow
free_lun, rest
free_lun, tbin
free_lun, fbin


RETURN

END
  
