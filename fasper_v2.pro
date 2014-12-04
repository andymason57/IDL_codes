pro SPREAD, y, yy, n, x, m
;
;   Given an array yy(0:n-1), extirpolate (spread) a value y into
;   m actual array elements that best approximate the "fictional"
;   (i.e., possible noninteger) array element number x.  The weights
;   used are coefficients of the Lagrange interpolating polynomial

         nfac=[0,1,1,2,6,24,120,720,5040,40320,362880]

         if (m gt 10) then message,'factorial table too small in spread'
         ix=long(x)
         if x eq float(ix) then yy(ix)=yy(ix)+y $
         else begin

              ilo= ( long(x-0.5*m+1.0) > 0 ) < (n-m)
              ihi=ilo+m-1
              nden=nfac(m)
              fac=x-ilo
              for j=ilo+1,ihi do fac = fac*(x-j)
              yy(ihi) = yy(ihi) + y*fac/(nden*(x-ihi))
              for j=ihi-1,ilo,-1 do begin
                   nden=(nden/(j+1-ilo))*(j-ihi)
                   yy(j) = yy(j) + y*fac/(nden*(x-j))
              endfor
         endelse
end
;
function sign,a,b
;
i1=where(b lt 0.0,COMPLEMENT=i2)
if(i1(0) ne -1) then a(i1)=-abs(a(i1))
if(i2(0) ne -1) then a(i2)=abs(a(i2))
f=a
return,f
end




pro FASPER_v2, other_args

;wk1, wk2, nout, jmax, prob

;FILE_PATH, GTI_FILE_PATH, OUTPUT_FILE_PATH, obsID, ofac, hifac, energy_range, detid




args = command_line_args()
FILE_PATH = args[0] 
GTI_FILE_PATH = args[1] 
OUTPUT_FILE_PATH = args[2]
obsID = args[3]
OFAC = args[4]
HIFAC = args[5]
energy_range = args[6]
detid = args[7]

;convert ofac and hifac from strings to integers
ofac = ofac+0
hifac = hifac+0



         MACC = 4  ;Number of interpolation points per 1/4 cycle
                   ;of highest frequency

; added by andy m. assign x and y 
time_dat = mrdfits(FILE_PATH, 1, hdr, columns=['TIME'])
flux_dat = mrdfits(FILE_PATH, 1, hdr, columns=['RATE'])
frmtm = sxpar(hdr, 'TIMEDEL')

print, "frametime passed ", frmtm

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

;   Check dimensions of input arrays
         n    = n_elements( x )
         sz   = size( y )
         if n ne sz(1) then message, 'Incompatible arrays.'

         nout=0.5*ofac*hifac*n
         print, nout
         nfreqt=long(ofac*hifac*n*MACC)     ;Size the FFT as next power
         nfreq=64L                          ;of 2 above nfreqt.

         while (nfreq lt nfreqt) do nfreq = ISHFT( nfreq,1 )
         ndim=long(ISHFT( nfreq,1 ))

         var=(stdev(y,ave))^2.              ;Compute the mean, variance
                                            ;and range of the data.
         xmin=MIN( x, MAX=xmax )
         xdif=xmax-xmin

         wk1 =complexarr(ndim)            ;Extirpolate the data into
         wk2 =complexarr(ndim)            ;the workspaces.
         fac=ndim/(xdif*ofac)
         fndim=ndim
         ck  =(x-xmin)*fac MOD fndim
         ckk =2.0*ck MOD fndim

         for j=0L,n - 1 do begin
              SPREAD, y(j)-ave,wk1,ndim,ck(j),MACC
              SPREAD, 1.0,wk2,ndim,ckk(j),MACC
         endfor

         wk1  = FFT( wk1,1, /OVERWRITE )    ;Take the Fast Fourier Transforms.
         wk2  = FFT( wk2,1, /OVERWRITE )

         wk1  = wk1(1:nout)
         wk2  = wk2(1:nout)
         rwk1 = double( wk1 ) & iwk1 = imaginary( wk1 )
         rwk2 = double( wk2 ) & iwk2 = imaginary( wk2 )
         df=1.0/(xdif*ofac)

         hypo2 = 2.0 * abs( wk2 )                  ;Compute the Lomb value for each
         hc2wt= rwk2/hypo2               ;frequency.
         hs2wt= iwk2/hypo2

         cwt  = sqrt(0.5+hc2wt)
         swt  = SIGN(sqrt(0.5-hc2wt),hs2wt)
         den  = 0.5*n+hc2wt*rwk2+hs2wt*iwk2
         cterm= (cwt*rwk1+swt*iwk1)^2./den
         sterm= (cwt*iwk1-swt*rwk1)^2./(n-den)

         wk1  = df*(findgen(nout)+1.)
         wk2  = (cterm+sterm)/(2.0*var)
         
         num_elements_wk1 = n_elements(wk1)
         num_elements_wk2 = n_elements(wk2)
         
         print, num_elements_wk1
         print, num_elements_wk2

         ; slice first 10 frequencies - see if high power red noise if not - leave LC alone, otherwise ignore first 10 freq
         redn_wk1 = wk1(0:9)
         redn_wk2 = wk2(0:9)
         
         red_pmax = MAX(redn_wk2, red_jmax)
         print, 'red_jmax', redn_wk2(red_jmax) 
         
         if redn_wk2(red_jmax) gt 10 then begin
          ;slice wk1 and wk2 
          new_wk1 = wk1(39:num_elements_wk1-1) 
          new_wk2 = wk2(39:num_elements_wk2-1)
          pmax = MAX( new_wk2, jmax )
          print, "red noise present"
         jmax = jmax + 40
         endif else begin
           pmax = MAX(wk2, jmax)
           print, 'red noise not present'
         endelse
                  
         expy =exp(-pmax)                   ;Estimate significance of largest
         effm =2.0*(nout)/ofac              ;peak value.
         prob =effm*expy
         if (prob gt 0.01) then prob=1.0-(1.0-expy)^effm
         
; added by andy m to output results

;insert code to write results to file
fname_freq = OUTPUT_FILE_PATH + obsID + "_" + detid + "_" + energy_range + '_freq.dat'
fname_pow = OUTPUT_FILE_PATH + obsID + "_" + detid + "_" + energy_range + '_power.dat'
fname_rest = OUTPUT_FILE_PATH + obsID  + "_" + detid + "_" + energy_range + '_restresults.dat'
fname_tbin = OUTPUT_FILE_PATH + obsID + "_" + detid + "_" + energy_range + '_tbin.dat'
fname_fbin = OUTPUT_FILE_PATH + obsID + "_" + detid + "_" + energy_range + '_fbin.dat'
openw, freq, fname_freq, /get_lun
openw, pow, fname_pow, /get_lun
openw, rest, fname_rest, /get_lun
openw, tbin, fname_tbin, /get_lun
openw, fbin, fname_fbin, /get_lun
printf, freq, wk1, format = '(f0.18)' 
printf, pow, wk2, format = '(f0.18)' 
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
         

end

;+
; NAME:
;        FASPER
;
; PURPOSE:
;        Given abscissas x (which need not be equally spaced) and ordinates
;        y, and given a desired oversampling factor ofac (a typical value
;        being 4 or larger). this routine creates an array wk1 with a
;        sequence of nout increasing frequencies (not angular frequencies)
;        up to hifac times the "average" Nyquist frequency, and creates
;        an array wk2 with the values of the Lomb normalized periodogram at
;        those frequencies.  The arrays x and y are not altered.  This
;        routine also returns jmax such that wk2(jmax) is the maximum
;        element in wk2, and prob, an estimate of the significance of that
;        maximum against the hypothesis of random noise. A small value of prob
;        indicates that a significant periodic signal is present.
;
; CATEGORY:
;        Numerical Recipes routines.
;
; CALLING SEQUENCE:
;
;        FASPER, X, Y, Ofac, Hifac, Wk1, Wk2 [, Nout, Jmax, Prob]
;
; INPUTS:
;        X:       Abscissas array, (e.g. an array of times).
;
;        Y:       Ordinates array, (e.g. an array of corresponding counts).
;
;        Ofac:     Oversampling factor.
;
;        Hifac:    Hifac * "average" Nyquist frequency = highest frequency
;                  for which values of the Lomb normalized periodogram will
;                  be calculated.
;
; OUTPUTS:
;        Wk1:      An array of Lomb periodogram frequencies.
;
;        Wk2:      An array of corresponding values of the Lomb periodogram.
;
; OPTIONAL OUTPUTS:
;        Nout:     The dimension of Wk1 and Wk2, i.e. the number of frequencies
;                  calculated.
;
;        Jmax:     The array index corresponding to the MAX( Wk2 ).
;
;        Prob:     The False Alarm Probability (FAP) of the largest value of
;                  of Lomb normalized periodogram.
;
; MODIFICATION HISTORY:
;        Written by:    Han Wen, December 1994. (Adapted from a Numerical
;                       Recipes routine with the same name).
;        15-AUG-1996    Andrew Lee, Modified spread and fasper to use arrays
;                       starting at 0 instead of 1, and fixed some bugs where
;                       int was used instead of long.
;-

; (C) Copr. 1986-92 Numerical Recipes Software
