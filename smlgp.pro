Function smlgp, rate, time, dt, psdt, altm=altm, cgti=cgti


; ----------------------------------------------------------------------

; Function to correct for small gaps in the light curve

; ----------------------------------------------------------------------

; Segment where the count rate is zero
  zum = where (rate eq 0)
  t_r = time[zum]
  rt_r = rate[zum]

  corrt = rate & n = 3 & altm = time
  
  print, "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'

  if n_elements(t_r) gt 1 then begin
      segs = ts_segment(t_r, dx=dt, minseg=floor(1000./dt), nseg=nsgs)

; ----------------------------------------------------------------------

; Correct these bins with the interpolated

      for i = 0, nsgs-1 do begin
          bg = where(time eq t_r[segs[i,0]])-1
          ed = where(time eq t_r[segs[i,1]])+1
          if bg eq -1 then begin
              bg = 0 
              n = n-1 
          endif
          if ed eq n_elements(corrt) then begin
              ed = n_elements(corrt)-1 
              n = n-1
          endif

          sta = rate[bg]
          sto = rate[ed]
          nvls = segs[i,1] - segs[i,0] + n
          nwvls = interpol([sta,sto], nvls)
          nwvls = poidev(nwvls)
          corrt[bg:ed] = mean(rt_r); nwvls
      endfor
  endif
  
; ------------------------------------------------------------------------

; Correct for small gaps within the GTI (Telemetry drop outs)

  segs = ts_segment(time, dx=dt)

  ngts = n_elements(segs[*,0])

  gti = time[segs]

  if ngts gt 1 then begin

  zum = make_array(ngts-1, /index)
  
  diffs = gti[zum+1,0] - gti[zum, 1]

  twil = where (diffs lt psdt/5., ct)

  if ct gt 0 then begin

      for i = 0, ct-1 do begin
          bg = where(time eq gti[twil[i]+1,0])-1
          ed = where(time eq gti[twil[i],1])+1
          if bg eq -1 then begin
              bg = 0 
              n = n-1 
          endif
          if ed eq n_elements(corrt) then begin
              ed = n_elements(corrt)-1 
              n = n-1
          endif
          sta = rate[bg]
          sto = rate[ed]
          nvls = diffs[twil[i]] + n
          nwvls = interpol([sta,sto], nvls)
          nwvls = poidev(nwvls)
          nels = n_elements(corrt)
          corrt = [corrt[0:sta-1], nwvls, corrt[sto+1:nels-1]]
          nwt = (make_array(nvls, /index)*dt) + time[bg-1]
          altm = [altm[0:sta-1], nwt, altm[sto+1:nels-1]]
          gti[twil[i]+1, 0] = -1
          gti[twil[i], 1] = -1
      endfor

      smcor = where(gti[*,0] ne -1, ct)
      cgti = fltarr(ct, 2)
      cgti[*,0] = gti[smcor,0]
      smcor = where(gti[*,1] ne -1, ct)
      cgti[*,1] = gti[smcor,1]

  endif
 
  endif else cgti = gti

  return, corrt

end

      


          