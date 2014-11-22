function limflrcorr, rate, time, dt, segs=segs, gpsz=gpsz, flrlim=flrlim
 
 
;err_log = '/mnt/4tbdata/error_log.txt' 
;CATCH, err
;If err NE 0 THEN BEGIN
   ;cancel the error handler for this procedure 
;  CATCH, /CANCEL 
  ; write out detectionID and error code/message of calling procedure 
;  openw, snow, err_log, /get_lun
;  printf, snow, 'OBSID : ', obsID_path
;  printf, snow, 'Detection ID: ', detID
;  printf, snow, 'Error index: ', err
;  printf, snow, 'Error message: ', !ERROR_STATE.MSG
;  print, 'limflrcorr runtime error handler!!!!!!!!!!!!!!!!!!!!'
;  free_lun, snow
;ENDIF  

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

  ;print, "RATE :", rate
  ;print, "FLRLIM :", flrlim

;Remove NaN from both rate and time arrays
  bad_points = where(finite(rate,/NAN), badcount, COMPLEMENT=good_points, NCOMPLEMENT=good_count)
  if good_count gt 0 then rate = rate[good_points] 
  
  mask = where(rate le flrlim)
  print, "mask: ", mask
  time = time[mask]
  
  print, "TIME : ", time 

; Figure out where there are non-consecutive time steps

  trk = make_array(n_elements(time)-1, /index)
  utrk = trk+1
  
  tdiffs = time[utrk]-time[trk]
  
  ;print, "tdiffs: ", tdiffs
  

; Define these as the gaps in the time series

  ; Add in a limit for when IDL is stupid
  lim = 1e-6
   gps = where(tdiffs gt dt*(1.0+lim))
  
; ------------------------------------------------------------------------

; Create a GTI array only for segments of good time which are large enough
  if gps[0] ne -1 then begin
      
     segs = make_array(n_elements(gps)+1, 2)
     segs[*,0] = [0,gps+1]
     segs[*,1] = [gps, n_elements(time)-1]
     
     gti = time[segs]

; Find the differences in time for the Good Time intervals
     diffgti = gti[*,1]-gti[*,0]
     
;     print, "diffgti: ", diffgti

    zum = where(diffgti gt gpsz)

; altered by me from gt to GE (greater than to greater than or equal)
;    zum = where(diffgti GE gpsz) 
;     print, "zum: ", zum    
 
; Filter the GTI based on this mask

     gti = gti[zum,*]
     segs = segs[zum,*]
endif else gti = [min(time), max(time)]

; -----------------------------------------------------------------------

return, gti

end
