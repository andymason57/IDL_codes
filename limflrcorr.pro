function limflrcorr, rate, time, dt, segs=segs, gpsz=gpsz, flrlim=flrlim
 
; -----------------------------------------------------------------------

; Produce the Flare GTIs given the SN corrected limit and using the
; in-band flare background timeseries 
;   INPUTS:
;          rate   - array of count rates
;          time   - time array to match the rates
;          dt     - time resolution of the time series
;          flrlim - S/N corrected limit to find flares
;          gpsz   - minimum size of good time intervals to be returned
;  
;   OUTPUTS: 
;          segs   - array of start/stop indices for the GTIs
;          GTI    - array of start/stop times for the GTIs


; ------------------------------------------------------------------------

 
;Remove NaN from both rate and time arrays
  bad_points = where(finite(rate,/NAN), badcount, COMPLEMENT=good_points, NCOMPLEMENT=good_count)
  if good_count gt 0 then rate = rate[good_points] 
  
  mask = where(rate le flrlim)
  time = time[mask]
  


; Figure out where there are non-consecutive time steps
  trk = make_array(n_elements(time)-1, /index)
  utrk = trk+1
  
  tdiffs = time[utrk]-time[trk]
  


; Define these as the gaps in the time series
  ; Add in a limit for when IDL is stupid
  lim = 1e-6
   gps = where(tdiffs gt dt*(1.0+lim))
;   print, "gps: ", n_elements(gps)
  
; ------------------------------------------------------------------------

; Create a GTI array only for segments of good time which are large enough
  if gps[0] ne -1 then begin
      
     segs = make_array(n_elements(gps)+1, 2)
     segs[*,0] = [0,gps+1]
     segs[*,1] = [gps, n_elements(time)-1]
     
     gti = time[segs]
     
; Find the differences in time for the Good Time intervals
     diffgti = gti[*,1]-gti[*,0]
     
    zum = where(diffgti gt gpsz)

; Filter the GTI based on this mask

     gti = gti[zum,*]
     segs = segs[zum,*]
endif else gti = [min(time), max(time)]

; -----------------------------------------------------------------------

return, gti

end
