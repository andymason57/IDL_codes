PRO alt_chimax, other_args

args = command_line_args()

;epfold,intime,inrate,raterr=inerror,pstart=pstart,pstop=pstop, $
;           nbins=nbins,sampling=sampl,chatty=chatty,chierg=chierg, $
;           maxchierg=maxchierg,gti=gti,persig=persig,linear=linear,        $
;           trial_period=trial_period,fitchi=fitchi, tolerance=tolerance, $
;           gap=gap,pdot=pdot,pddot=pddot,time0=time0,fchatty=fchatty,mjd=mjd
;           
;           
;+
; NAME:
;             epfold
;
;
; PURPOSE: 
;
; CATEGORY: 
;             timing tools
;
;
; CALLING SEQUENCE:
;             epfold,time,rate,raterr=raterr,pstart=pstart,pstop=pstop, 
;                    nbins=nbins,sampling=sampling,/chisq
;                       
; 
; INPUTS:
;             time : a vector containing the time in arbitary units 
;                    (usually to be assumed seconds)
;             rate : a vector containing the countrate; if not given,
;                    the times are assumed to come from individual events
;             pstart:   lowest period to be considered
;             pstop:    highest period to be considered
;
; OPTIONAL INPUTS:
;             gti: at the moment for event data only: gti[2,*] is an
;                  array with the good times from the event selection,
;                  gti[0,*] are the start times, gti[1,*] the stop
;                  times. Obviously, use the same time system for time
;                  and gti (i.e., also barycenter gtis!). Needed for
;                  the computation of the exposure time of the phase
;                  bins, if not given, the time of the first and the
;                  last event is used instead.
;             gap: array containing bad bins (e.g., gaps in the
;                  lightcurve). computed in pfold with timegap or from
;                  the gtis if not given, and possibly returned if an
;                  undeclared variable is passed through the gap keyword
;             raterror: a given error vector for the countrate. If not
;                       given and we are not working on events,
;                       raterror is computed using Poissonian
;                       statistics.
;             nbins:    number of phase-bins to be used in creating the trial
;                       pulse (default=20)
;             sampling: how many periods per peak to use (default=10)   
;             linear  : if set, use a linear trial period array with
;                       step equal to the value of linear
;             trial_period : array contains the trial periods. If set,
;                            pstart is the smallest trial period and
;                            pstop the longest. 
;            tolerance: parameter defining the lower limit for the gap
;                       length; the reference is the time difference
;                       between the first and second entry in the time
;                       array; tolerance defines the maximum allowed relative
;                       deviation from this reference bin length; 
;                       default: 1e-8; this parameter is passed to timegap
;                       (see timegap.pro for further
;                       explanation). Only used if no gti is given. 
;             pdot,pddot: period derivative and 2nd period derivative
;                       (NB: these are NOT stepped, so they make only
;                       sense, e.g., when epfold is called from a
;                       two-parametric period search)
;             time0:    epoch for folding
;
; KEYWORD PARAMETERS:
;             fitchi   : fit a Gauss distribution to the chi^2 and use
;                        the center of the gauss to determine the
;                        period, instead of using the maximum of the
;                        chi^2 distribution. This is mainly needed by
;                        eperror.pro. In case this keyword is set
;                        chierg will be 3 dimensional array and the
;                        fit result is stored in
;                        chierg[2,*]. maxchierg[1] still contains the
;                        maximum chi^2 and NOT the chi^2 of the fitted
;                        period !!
;             chatty   : if set, be chatty
;             fchatty  : chatty keyword for pfold
;             mjd      : if set, time is assumed in days, but the
;                        period will be in seconds anyway (apart from
;                        mjd being set, period is in the same units as time)
;   
; OUTPUTS:
;             chierg   : 2D-array, chierg[0,*] contains the trial
;                        period, and chierg[1,*] the respective chi^2
;                        value. 
;
; OPTIONAL OUTPUTS:
;             maxchierg: 2 dim. array, maxchierg[0] contains the
;                        period of maximum chi^2, and maxchierg[1] the
;                        respective maximum chi^2
;             persig   : uncertainty of that period, from triangular
;                        approximation, not to be interpreted in a
;                        statistical sense!
;             gap      : see above (optional inputs)
;
;
; COMMON BLOCKS:
;             none
;
;
; SIDE EFFECTS:
;             none
;
;
; RESTRICTIONS:
;
;             the input lightcurve has to be given in count rates (not
;             in photon numbers). To prevent convergence problems of
;             the fit (only if keyword /fitchi is set), the maximum of
;             the chi^2 distribution has to be in between pstart and
;             pstop, and the period interval has to be resonable
;             large.
;
;
; PROCEDURE:
;             This subroutine performs a period search using epoch
;             folding. For each trial period, a pulse profile is
;             generated using pfold. This profile is then tested for
;             constancy using a chi^2 test. The maximum chi^2,
;             i.e., the maximum deviation, is the most possible
;             period. This is done for all periods between pstart and
;             pstop using a grid-search operating on the minimum
;             possible period that can still be detected (given by
;             p^2/t). 
;
;             Caveat: The significance of the found periods SHOULD NOT
;             be tested using the chi^2 value obtained from this
;             routine, since the values are only asymptotically chi^2
;             distributed (see Davies, 1990, Schwarzenberg-Czerny,
;             1989, Larsson, 1996)
; 
;             Read (and understand) the references before using this
;             routine!
;
;             References:
;                 Davies, S.R., 1990, MNRAS 244, 93
;                 Larsson, S., 1996, A&AS 117, 197
;                 Leahy, D.A., Darbro, W., Elsner, R.F., et al., 1983,
;                    ApJ 266, 160-170
;                 Schwarzenberg-Czerny, A., 1989, MNRAS 241, 153
;
; EXAMPLE:
;
;
; 
   
   ;; ******* Set default values *********
   
   time=intime

   ;; set default value sampling per peak
   IF (n_elements(sampl) EQ 0) THEN sampl =10.              
   
   ;; set default value number of bins
   IF (n_elements(nbins) EQ 0) THEN nbins=15

   ;; set default tolerance:
   IF (n_elements(nbins) EQ 0) THEN tolerance=1e-8 

   ;; set default epoch
   IF (n_elements(time0) EQ 0) THEN time0=time[0]


   dim=n_elements(time) 

   IF (keyword_set(chatty)) THEN BEGIN 
       print,' '
       print,' >> Starting data analysis <<'
   ENDIF 

   ;; ******* DATA ANALYSIS *******
   IF ( n_elements(inrate) EQ 0) THEN BEGIN 

       ;; event data: no gti file known -> warning
       IF n_elements(gti) EQ 0 THEN BEGIN 
           message,'eventdata: GTI not known, using time of first',/informational
           message,'  and last event. This is most probably NOT',/informational
           message,'  what you want!',/informational
           gti=fltarr(2,1)
           gti[0,0]=time[0]
           gti[1,0]=time[n_elements(time)-1]
       ENDIF 

       ;; duration for events 
       ;; NB: this is total DURATION, not total EXPOSURE!
       IF size(gti,/n_dimensions) EQ 2 THEN ngti=(size(gti))[2] ELSE ngti=1
       tges=gti[1,ngti-1]-gti[0,0]

       IF (keyword_set(mjd)) THEN tges=tges*86400d0


   ENDIF ELSE BEGIN 

       rate=inrate
   
       ;; calculate error or use given error-vector
       IF (n_elements(inerror) EQ 0) THEN BEGIN 
           ;;JW are we sure we want to compute the error this way?
           ;;JW if we have negative rates, then the data have been 
           ;;JW background subtracted and the following formula does
           ;;JW not apply...
           error=sqrt(abs(rate))
       END ELSE BEGIN 
           error=inerror
       ENDELSE  

       mittel = total(rate)/float(dim) ;; calculate medium value
       rate = rate-mittel              ;; renormalize to mean countrate
       var = variance(rate)            ;; sample variance

       ;; Time resolution
       delta_t = temporary(shift(time,-1))-time
       delta_t[n_elements(time)-1]=delta_t[0] ;; somewhat arbitrary
       ;; calculate smallest delta_t and set it to bintime
       bintime=min(delta_t)
       ;; Max. Timescale = pdim*dtstep = TGes
       tges = time[dim-1]-time[0]+ bintime ; 

       IF (keyword_set(mjd)) THEN BEGIN 
           bintime=bintime*86400d0
           tges=tges*86400d0
       ENDIF 

       
   ENDELSE
  
   IF (keyword_set(chatty)) THEN BEGIN 
       IF ( n_elements(inrate) EQ 0) THEN BEGIN 
           print,' You work with event data '
           print,' total observation time ',tges
       ENDIF ELSE BEGIN
           rms = sqrt(total(rate^2.)/(dim-1.))
           print,' mean / rms             ',mittel,rms ; medium countrate
           print,' total observation time ',tges
           print,' effective observation  ',dim*bintime
           print,' duty cycle             ',100.*(dim*bintime)/tges,' %'
           print
       ENDELSE 
   ENDIF 
   
   ;; ******* CHI^2-statistics *******
   
    ;; Determine final dimension of result
   IF n_elements(trial_period) EQ 0 THEN BEGIN 
       IF (keyword_set(chatty)) THEN BEGIN 
           print,'Determining number of trial periods'
           print,'   pstart',pstart
           print,'   pstop ',pstop
       ENDIF 

       ergdim=0L
       p = pstart
       IF n_elements(linear) EQ 0 THEN BEGIN 
           WHILE (p LE pstop) DO BEGIN
               ergdim = ergdim + 1L
               p=p+p*p/(tges*sampl)
           ENDWHILE
       ENDIF ELSE BEGIN 
           ergdim = long((pstop - pstart)/double(linear)) +1L
       ENDELSE 
       IF (keyword_set(chatty)) THEN BEGIN
           print,'number of periods',ergdim
           print,'  ... done'
       ENDIF 
   ENDIF ELSE BEGIN 
       ;; sort the array in ascending order 
       trial_period = trial_period[sort(trial_period)]
       ;; dimension of the result
       ergdim = n_elements(trial_period)
       ;; determine the first and last trial period
       pstart = trial_period[0]
       pstop = trial_period[ergdim-1]
   ENDELSE

   IF (keyword_set(fitchi)) THEN BEGIN 
       ;; We need a 3 dimensional array to store the chi^2 fit 
       chierg=dblarr(3,ergdim)
   ENDIF ELSE BEGIN 
       chierg=dblarr(2,ergdim)
   ENDELSE

   IF (keyword_set(chatty)) THEN BEGIN 
       tstart=systime(1)
       tlast=systime(1)
   ENDIF 

   ptest = pstart
   i=0L
   WHILE (ptest LE pstop) DO BEGIN
       chierg[0,i]=ptest
       IF ( n_elements(inrate) EQ 0) THEN BEGIN 
           ;; Work with event data
           pfold, time, profile, gti=gti,period=ptest,nbins=nbins, $
             ttot=ttot,npts=npts,pdot=pdot,pddot=pddot,time0=time0,$
             chatty=fchatty,mjd=mjd
           ndx=where(finite(profile))
           IF (i EQ 0) THEN BEGIN 
               meanrate=mean(profile[ndx])
           ENDIF 
           IF ( ndx[0] NE -1 ) THEN BEGIN 
               ;; variance for a Poisson process is just the rate
               chierg[1,i]=total((profile[ndx]-meanrate)^2. / meanrate)
           ENDIF ELSE BEGIN 
               chierg[1,i]=!values.f_nan
           ENDELSE 
       ENDIF ELSE BEGIN 
           ;; compute profile for testperiod ptest.
           ;; if gti is given and gaps are not known, computes gap
           ;; if gap is known, it is used for determining the bad bins
           ;; of the lightcurve
           pfold,time,rate,profile,period=ptest,nbins=nbins, $
             proferr=proferr,raterr=error,npts=npts, $
             dt=delta_t,tolerance=tolerance,gti=gti,gap=gap, $
             pdot=pdot,pddot=pddot,time0=time0,chatty=fchatty,mjd=mjd
       
           ;;JW  this is  Davies, eq.4 / Larsson, eq.1
           ndx=where(finite(profile))
           IF (ndx[0] NE -1) THEN BEGIN 
               chierg[1,i]=total(npts*profile[ndx]^2.)/var
           ENDIF ELSE BEGIN
               chierg[1,i]=!values.f_nan
           ENDELSE 
       ENDELSE
       ;; Next test period
       IF n_elements(trial_period) EQ 0 THEN BEGIN 
           IF n_elements(linear) EQ 0 THEN BEGIN 
               ptest = ptest + ptest*ptest/(tges*sampl)
           ENDIF 
           IF n_elements(linear) NE 0 THEN BEGIN 
               ptest = ptest + double(linear)      
           ENDIF   
       ENDIF ELSE BEGIN 
           IF i+1L LT ergdim THEN ptest = trial_period[i+1L] $
           ELSE ptest = pstop + 1.
       ENDELSE  
       i=i+1L
       
       ;; Current status
       IF (keyword_set(chatty)) THEN BEGIN 
           IF (systime(1)-tlast GE 10) THEN BEGIN
               tlast=systime(1)
               perc=i*100./ergdim
               titer=(tlast-tstart)/double(i)
               remain=(ergdim-i)*titer
               print,strtrim(string(i*100./ergdim),2)+'% done'
               print,'   '+strtrim(string(titer),2)+'s per iteration'
               IF (remain LE 600) THEN BEGIN 
                   print,'   '+strtrim(string(format='(F9.1)',remain),2)+'s remaining'
               ENDIF ELSE BEGIN 
                   IF (remain LE 7200) THEN BEGIN 
                       print,'   '+strtrim(string(format='(F7.1)',remain/60.),2)+'min remaining'
                   ENDIF ELSE BEGIN 
                       print,'   '+strtrim(string(format='(F7.1)',remain/3600.),2)+'h remaining'
                   ENDELSE 
               ENDELSE 
           ENDIF 
       ENDIF 
   ENDWHILE

   ndx=where(finite(chierg[1,*]))
   IF (ndx[0] EQ -1) THEN BEGIN 
       ;;JW I do not see how we could fold data and not get a single
       ;;JW chi^2 value, so we should never get here....
       message,'This should not happen'
   ENDIF 

   eer=n_elements(ndx)
   chimean=total(chierg[1,ndx])/eer
   chisig=sqrt(total((chierg[1,ndx]-chimean)^2.)/(eer-1))
   maxchi=max(chierg[1,ndx],maxindex)

   IF (keyword_set(fitchi)) THEN BEGIN 
       ;; create estimates using currently known values:

       ;; save us some indirect indexing via ndx
       per=chierg[0,ndx]
       ccc=chierg[1,ndx]

       period = per[maxindex]   ; first estimate of period
       sy     = sqrt(ccc)       ; error bar (well, not really...)
       lin    = 0.d             ; constant

       p0=dblarr(4)             ; Coefficients for fitting function
       p0[0] = lin              ; mean constant 
       p0[1] = period           ; first estimate of period
       p0[2] = 0.1*period*period/tges ; sigma
       p0[3] = max(ccc)         ; area

       ;; Fit function looks like
       expr = 'P[0] + GAUSS1(x, P[1:3])' 
       ;; Get final cofficients
       p = mpfitexpr(expr, per,ccc, sy,p0,/quiet)  

       ;; compute fitted function::
       x = per
       ok = execute('chierg[2,*]='+expr)

       ;; Period:
       period=p[1]
   ENDIF ELSE BEGIN 
       period=chierg[0,ndx[maxindex]]
   ENDELSE 
   
   maxchierg = [period,maxchi]
 
   persig=period^2./tges
   
   IF (keyword_set(chatty)) THEN BEGIN 
       print,'Total time: '+strtrim(string(systime(1)-tstart),2)+'s'
       print,' Max. Chisquare : '+strtrim(string(maxchi),2)+ $
         ' ( P= '+strtrim(string(period),2)+')'
   ENDIF 
END
